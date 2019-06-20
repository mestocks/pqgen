#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pq_htable.h"
#include "pq_genetics.h"
#include "pq_generics.h"
#include "pq_args.h"
// values = [nvcodons, nsyn, nnon, ds, dn, ps, pn]
// outs   = [nsam, nvcodons, nsyn, nnon, ds, dn, ps, pn]


void pq_swupdate_pnds(struct StatObject *stats, char **array)
{
  int k;
  int diff_syn;
  int diff_nsyn;
  void *aptr;
  void *ref_aa;
  void *ref_sptr;
  double syn_muts;
  double nsyn_muts;

  unsigned int bad_codons;
  bad_codons = 0;
  for (k = 0; k < NKARGS; k++) {
    bad_codons += 1 - pq_alldna(array[KCOLS[k]]);
  }

  if (bad_codons == 0) {
    pq_dna_upper(array[KCOLS[0]]);
    ref_aa = pq_lookup_hash(&CODON_TO_AMINO, array[KCOLS[0]]);
    ref_sptr = pq_lookup_hash(&CODON_TO_NSYN, array[KCOLS[0]]);
    syn_muts = *(double *)ref_sptr;
    nsyn_muts = 9.0 - syn_muts;
    
    *(long long int *)stats->values[0] += 1;
    *(double *)stats->values[1] += syn_muts / 3.0;
    *(double *)stats->values[2] += nsyn_muts / 3.0;
    
    diff_syn = 0;
    diff_nsyn = 0;
    for (k = 1; k < NKARGS; k++) {
      pq_dna_upper(array[KCOLS[k]]);
      aptr = pq_lookup_hash(&CODON_TO_AMINO, array[KCOLS[k]]);
      if (strcmp(array[KCOLS[0]], array[KCOLS[k]]) == 0) {
	
      } else {
	if (strcmp(aptr, ref_aa) != 0) {
	  diff_nsyn++;
	} else {
	  diff_syn++;
	}
      }
    }
    
    if (stats->nsam == diff_syn) {
      *(long long int *)stats->values[3] += 1;
    } else if (diff_syn > 0) {
      *(long long int *)stats->values[5] += 1;
    }
    if (stats->nsam == diff_nsyn) {
      *(long long int *)stats->values[4] += 1;
    } else if (diff_nsyn > 0) {
      *(long long int *)stats->values[6] += 1;
    }
  }
}


void pq_swwrite_pnds(struct StatObject *stats)
{
  sprintf(stats->outs[0], "%d", stats->nsam);
  sprintf(stats->outs[1], "%lli", *(long long int *)stats->values[0]);
  sprintf(stats->outs[2], "%f", *(double *)stats->values[1]);
  sprintf(stats->outs[3], "%f", *(double *)stats->values[2]);
  sprintf(stats->outs[4], "%lli", *(long long int *)stats->values[3]);
  sprintf(stats->outs[5], "%lli", *(long long int *)stats->values[4]);
  sprintf(stats->outs[6], "%lli", *(long long int *)stats->values[5]);
  sprintf(stats->outs[7], "%lli", *(long long int *)stats->values[6]);
}


void pq_swclear_pnds(struct StatObject *stats)
{
  int i;
  for (i = 0; i < stats->nvalues; i++) {
    *(long long int *)stats->values[i] = 0;
  }
}


void PQ_PNDS_INIT(struct StatObject *stats, int nsam)
{
  int i;
  stats->nsam = nsam - 1;
  stats->nvalues = 7;
  stats->nouts = 8;
  stats->outs = calloc(stats->nouts, sizeof(char *));
  stats->values = calloc(stats->nvalues, sizeof(void *));
  for (i = 0; i < stats->nouts; i++) {
    stats->outs[i] = calloc(128, sizeof(char));
  }
  stats->values[0] = malloc(sizeof(long long int));
  *(long long int *)stats->values[0] = 0;
  stats->values[1] = malloc(sizeof(double));
  *(double *)stats->values[1] = 0.0;
  stats->values[2] = malloc(sizeof(double));
  *(double *)stats->values[2] = 0.0;
  for (i = 3; i < stats->nvalues; i++) {
    stats->values[i] = malloc(sizeof(long long int));
    *(long long int *)stats->values[i] = 0;
  }
  stats->update = pq_swupdate_pnds;
  stats->write = pq_swwrite_pnds;
  stats->clear = pq_swclear_pnds;
}

