#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_htable.h>
#include <pq_genetics.h>
#include <pq_generics.h>
#include <pq_args.h>
// values = [nvcodons, nsyn, nnon, ds, dn, ps, pn]
// outs   = [nsam, nvcodons, nsyn, nnon, ds, dn, ps, pn]


void pq_swupdate_pnds(struct SWrap *wrap, char **array)
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
    ref_aa = rwk_lookup_hash(&CODON_TO_AMINO, array[KCOLS[0]]);
    ref_sptr = rwk_lookup_hash(&CODON_TO_NSYN, array[KCOLS[0]]);
    syn_muts = *(double *)ref_sptr;
    nsyn_muts = 9.0 - syn_muts;
    
    *(long long int *)wrap->values[0] += 1;
    *(double *)wrap->values[1] += syn_muts / 3.0;
    *(double *)wrap->values[2] += nsyn_muts / 3.0;
    
    diff_syn = 0;
    diff_nsyn = 0;
    for (k = 1; k < NKARGS; k++) {
      pq_dna_upper(array[KCOLS[k]]);
      aptr = rwk_lookup_hash(&CODON_TO_AMINO, array[KCOLS[k]]);
      if (strcmp(array[KCOLS[0]], array[KCOLS[k]]) == 0) {
	
      } else {
	if (strcmp(aptr, ref_aa) != 0) {
	  diff_nsyn++;
	} else {
	  diff_syn++;
	}
      }
    }
    
    if (wrap->nsam == diff_syn) {
      *(long long int *)wrap->values[3] += 1;
    } else if (diff_syn > 0) {
      *(long long int *)wrap->values[5] += 1;
    }
    if (wrap->nsam == diff_nsyn) {
      *(long long int *)wrap->values[4] += 1;
    } else if (diff_nsyn > 0) {
      *(long long int *)wrap->values[6] += 1;
    }
  }
}


void pq_swwrite_pnds(struct SWrap *wrap)
{
  sprintf(wrap->outs[0], "%d", wrap->nsam);
  sprintf(wrap->outs[1], "%lli", *(long long int *)wrap->values[0]);
  sprintf(wrap->outs[2], "%f", *(double *)wrap->values[1]);
  sprintf(wrap->outs[3], "%f", *(double *)wrap->values[2]);
  sprintf(wrap->outs[4], "%lli", *(long long int *)wrap->values[3]);
  sprintf(wrap->outs[5], "%lli", *(long long int *)wrap->values[4]);
  sprintf(wrap->outs[6], "%lli", *(long long int *)wrap->values[5]);
  sprintf(wrap->outs[7], "%lli", *(long long int *)wrap->values[6]);
}


void pq_swclear_pnds(struct SWrap *wrap)
{
  int i;
  for (i = 0; i < wrap->nvalues; i++) {
    *(long long int *)wrap->values[i] = 0;
  }
}


void PQ_PNDS_INIT(struct SWrap *wrap, int nsam)
{
  int i;
  wrap->nsam = nsam - 1;
  wrap->nvalues = 7;
  wrap->nouts = 8;
  wrap->outs = calloc(wrap->nouts, sizeof(char *));
  wrap->values = calloc(wrap->nvalues, sizeof(void *));
  for (i = 0; i < wrap->nouts; i++) {
    wrap->outs[i] = calloc(128, sizeof(char));
  }
  wrap->values[0] = malloc(sizeof(long long int));
  *(long long int *)wrap->values[0] = 0;
  wrap->values[1] = malloc(sizeof(double));
  *(double *)wrap->values[1] = 0.0;
  wrap->values[2] = malloc(sizeof(double));
  *(double *)wrap->values[2] = 0.0;
  for (i = 3; i < wrap->nvalues; i++) {
    wrap->values[i] = malloc(sizeof(long long int));
    *(long long int *)wrap->values[i] = 0;
  }
  wrap->update = pq_swupdate_pnds;
  wrap->write = pq_swwrite_pnds;
  wrap->clear = pq_swclear_pnds;
}

