#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_htable.h>
#include <pq_generics.h>
#include <pq_args.h>

void pq_swupdate_div(struct SWrap *wrap, char **array)
{

  char *ref;
  char *alt;

  ref = array[KCOLS[0]];
  alt = array[KCOLS[1]];

  if (strcmp(ref, ".") != 0 && strcmp(alt, ".") != 0) {
    if (strcmp(ref, alt) != 0) {
      *(long long int *)wrap->values[0]+=1;
    }
    *(long long int *)wrap->values[1]+=1;
  }
}

void pq_swwrite_div(struct SWrap *wrap)
{
  long long int diffs, nvsites;
  diffs =  *(long long int *)wrap->values[0];
  nvsites =  *(long long int *)wrap->values[1];
  sprintf(wrap->outs[0], "%lli", nvsites);
  if (strcmp((char *)rwk_lookup_hash(&ARGHASH, "-b"), "0") == 0) {
    sprintf(wrap->outs[1], "%f", (double)diffs);
  } else {
    sprintf(wrap->outs[1], "%f", (double)diffs / (double)nvsites);
  }
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
