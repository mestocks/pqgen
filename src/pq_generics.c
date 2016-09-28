#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_args.h>
#include <rwk_htable.h>
#include "pq_generics.h"
#include "pq_sfstats.h"

void pq_init_parameters(struct pq_parameters *params)
{
  rwk_create_hash(&params->arghash, 128);
  params->KCOLS = NULL;
  params->update = pq_update_parameters;
}

void pq_free_parameters(struct pq_parameters *params)
{
  free(params->KCOLS);
  rwk_free_hash(&params->arghash);
}

void pq_update_parameters(struct pq_parameters *params, unsigned int argc, char **argv)
{  
  char *chr_ptr;
  char *end_ptr;
  char *kcol_ptr;
  char *fcol_ptr;
  
  rwk_multi_insert_hash(&params->arghash, argc, argv);
  chr_ptr = (char *)rwk_lookup_hash(&params->arghash, "-c");
  end_ptr = (char *)rwk_lookup_hash(&params->arghash, "-p");
  kcol_ptr = (char *)rwk_lookup_hash(&params->arghash, "-k");
  fcol_ptr = (char *)rwk_lookup_hash(&params->arghash, "-f");
  
  params->CHROM = atoi(chr_ptr) - 1;
  params->POS = atoi(end_ptr) - 1;
  params->FCOL = atoi(fcol_ptr) - 1;
  
  if (params->KCOLS != NULL) {
    free(params->KCOLS);
  }

  params->nkargs = rwk_countkargs(kcol_ptr);
  params->KCOLS = malloc(params->nkargs * sizeof(int));
  rwk_kargsplit(params->KCOLS, kcol_ptr);
}

void pq_print_bed(struct Data *data)
{
  printf("%s\t%llu\t%llu\t%s\t%llu\t%s\n", data->chr, *data->start, *data->end, data->name, *data->nvsites, data->value);
}

void pq_print_csv(struct Data *data)
{
  printf("%s,%llu,%llu,%s,%llu,%s\n", data->chr, *data->start, *data->end, data->name, *data->nvsites, data->value);
}


void pq_swfree(struct SWrap *wrap)
{
  int i;
  for (i = 0; i < wrap->nouts; i++) {
    free(wrap->outs[i]);
  }
  for (i = 0; i < wrap->nvalues; i++) {
    free(wrap->values[i]);
  }
  free(wrap->outs);
  free(wrap->values);
}


void nref_nalt(void **counts, char **array, struct pq_parameters *params)
{
  int i;
  char *curr;
  int nuniq_alleles;
  int allele_index;
  int allele_counts[512];
  for (i = 0; i < 512; i++) {
    allele_counts[i] = 0;
  }

  *(long long int *)counts[0] = 0;
  *(long long int *)counts[1] = 0;
  
  for (i = 0; i < params->nkargs; i++) {
    curr = array[params->KCOLS[i]];
    if (strlen(curr) != 3 || curr[1] != '/') {
      *(long long int *)counts[0] = 0;
      *(long long int *)counts[1] = 0;
      return;
    } else if (strlen(curr) == 1) {
      // does not work (returns 0). this is because wrap->nsam is set to params->nkargs
      allele_index = (int)curr[0];
      allele_counts[allele_index]++;
    } else if (strlen(curr) == 3) {
      allele_index = (int)curr[0];
      allele_counts[allele_index]++;
      allele_index = (int)curr[2];
      allele_counts[allele_index]++;
    }
  }  
  
  nuniq_alleles = 0;
  for (i = 0; i < 512; i++) {
    if (allele_counts[i] != 0) {
      if (nuniq_alleles > 1) {
	*(long long int *)counts[0] = 0;
	*(long long int *)counts[1] = 0;
	break;
      } else {
	if ((char)i != '.') {
	  *(long long int *)counts[nuniq_alleles] = allele_counts[i];
	  nuniq_alleles++;
	}
      }
    }
  }
}

void pq_swupdate_theta(struct SWrap *wrap, char **array, struct pq_parameters *params)
{
  int s, nref, nalt, nminor;
  s = nref = nalt = nminor= 0;

  nref_nalt(&wrap->values[3], array, params);
  nref = *(long long int *)wrap->values[3];
  nalt = *(long long int *)wrap->values[4];
  
  if (nref + nalt == wrap->nsam) {
    if (nref > 0 && nalt > 0) {
      s = 1;
    } else {
      s = 0;
    }
    *(long long int *)wrap->values[1] += s;
    if (nalt < nref) {
      nminor = nalt;
    } else {
      nminor = nref;
    }
    *(long long int *)wrap->values[2] += pqPairwiseDiffs(wrap->nsam, nminor);
    *(long long int *)wrap->values[0] += 1;
  }
}

void pq_swwrite_theta(struct SWrap *wrap, struct pq_parameters *params)
{
  long long int s, pisum, nvsites;
  double tw, combs, pi, tajd;
  nvsites = *(long long int *)wrap->values[0];
  s =  *(long long int *)wrap->values[1];
  pisum =  *(long long int *)wrap->values[2];
  tw = pqWattersonsTheta(wrap->nsam, s);
  combs = pqPairwiseCombs(wrap->nsam);
  pi = pqTajimasTheta(combs, pisum);
  tajd = pqTajimasD(wrap->nsam, s, tw, pi);
  sprintf(wrap->outs[0], "%d", wrap->nsam);
  sprintf(wrap->outs[1], "%lli", nvsites);
  sprintf(wrap->outs[2], "%lli", s);
  sprintf(wrap->outs[3], "%f", tw);
  sprintf(wrap->outs[4], "%f", pi);
  sprintf(wrap->outs[5], "%f", tajd);
}

void pq_swclear_theta(struct SWrap *wrap)
{
  int i;
  for (i = 0; i < wrap->nvalues; i++) {
    *(long long int *)wrap->values[i] = 0;
  }
}


void pq_swinit_theta(struct SWrap *wrap, int nsam)
{
  int i;
  wrap->nsam = nsam;
  wrap->nouts = 6;
  wrap->nvalues = 5;
  wrap->outs = calloc(wrap->nouts, sizeof(char *));
  wrap->values = calloc(wrap->nvalues, sizeof(void *));
  for (i = 0; i < wrap->nouts; i++) {
    wrap->outs[i] = calloc(128, sizeof(char));
  }
  for (i = 0; i < wrap->nvalues; i++) {
    wrap->values[i] = malloc(sizeof(long long int));
    *(long long int *)wrap->values[i] = 0;
  }
  wrap->update = pq_swupdate_theta;
  wrap->write = pq_swwrite_theta;
  wrap->clear = pq_swclear_theta;
}


void pq_swupdate_div(struct SWrap *wrap, char **array, struct pq_parameters *params)
{

  char *ref;
  char *alt;

  ref = array[params->KCOLS[0]];
  alt = array[params->KCOLS[1]];

  if (strcmp(ref, ".") != 0 && strcmp(alt, ".") != 0) {
    if (strcmp(ref, alt) != 0) {
      *(long long int *)wrap->values[0]+=1;
    }
    *(long long int *)wrap->values[1]+=1;
  }
}

void pq_swwrite_div(struct SWrap *wrap, struct pq_parameters *params)
{
  long long int diffs, nvsites;
  diffs =  *(long long int *)wrap->values[0];
  nvsites =  *(long long int *)wrap->values[1];
  sprintf(wrap->outs[0], "%lli", nvsites);
  sprintf(wrap->outs[1], "%f", (double)diffs / (double)nvsites);
}

void pq_swclear_div(struct SWrap *wrap)
{
  int i;
  for (i = 0; i < wrap->nvalues; i++) {
    *(long long int *)wrap->values[i] = 0;
  }
}

void pq_swinit_div(struct SWrap *wrap, int nsam)
{
  // nvsites, diffs
  int i;
  wrap->nsam = nsam;
  wrap->nouts = 2;
  wrap->nvalues = 2;
  wrap->outs = calloc(wrap->nouts, sizeof(char *));
  wrap->values = calloc(wrap->nvalues, sizeof(void *));
  for (i = 0; i < wrap->nouts; i++) {
    wrap->outs[i] = calloc(128, sizeof(char));
  }
  for (i = 0; i < wrap->nvalues; i++) {
    wrap->values[i] = malloc(sizeof(long long int));
    *(long long int *)wrap->values[i] = 0;
  }
  wrap->update = pq_swupdate_div;
  wrap->write = pq_swwrite_div;
  wrap->clear = pq_swclear_div;
}
