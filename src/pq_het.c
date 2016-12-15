#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_htable.h>

#include <pq_args.h>
#include <pq_generics.h>
#include <pq_genetics.h>


void pq_swupdate_het(struct SWrap *wrap, char **array)
{
  int nhet, nref, nalt;
  double het;
  
  pq_gtstat(&wrap->values[2], array);
  nhet = *(int *)wrap->values[2];
  nref = *(int *)wrap->values[3];
  nalt = *(int *)wrap->values[4];

  if (nref + nalt > 0) {
    het = (double)nhet / (((double)nref + (double)nalt) / 2);
    *(double *)wrap->values[1] += het;
    *(long long int *)wrap->values[0] += 1;
  }
}


void pq_swwrite_het(struct SWrap *wrap)
{
  double het;
  sprintf(wrap->outs[0], "%d", wrap->nsam);
  sprintf(wrap->outs[1], "%lli", *(long long int *)wrap->values[0]);

  if (strcmp((char *)rwk_lookup_hash(&ARGHASH, "-b"), "0") == 0) {
    het = *(double *)wrap->values[1];
  } else {
    het = *(double *)wrap->values[1] / *(long long int *)wrap->values[0];
  }
  sprintf(wrap->outs[2], "%f", het);
}


void pq_swclear_het(struct SWrap *wrap)
{
  *(long long int *)wrap->values[0] = 0;
  *(double *)wrap->values[1] = 0.0;
}


// nvalues = [nvsites, het, nhet, nref, nalt]
// nouts   = [nsam, nvsites, uhet]

void pq_het_init(struct SWrap *wrap, int nsam)
{
  int i;
  wrap->nsam = nsam;
  wrap->nouts = 3;
  wrap->nvalues = 5;
  wrap->outs = calloc(wrap->nouts, sizeof(char *));
  wrap->values = calloc(wrap->nvalues, sizeof(void *));
  for (i = 0; i < wrap->nouts; i++) {
    wrap->outs[i] = calloc(128, sizeof(char));
  }
  
  wrap->values[0] = malloc(sizeof(long long int));
  *(long long int *)wrap->values[0] = 0;

  wrap->values[1] = malloc(sizeof(double));
  *(double *)wrap->values[1] = 0.0;
  
  for (i = 2; i < wrap->nvalues; i++) {
    wrap->values[i] = malloc(sizeof(int));
    *(int *)wrap->values[i] = 0;
  }
  
  wrap->update = pq_swupdate_het;
  wrap->write = pq_swwrite_het;
  wrap->clear = pq_swclear_het;
}

