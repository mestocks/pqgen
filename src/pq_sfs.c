#include <stdio.h>
#include <stdlib.h>

#include <pq_args.h>
#include <pq_generics.h>
#include <pq_theta.h>


void pq_swupdate_sfs(struct SWrap *wrap, char **array)
{
  int nref, nalt;
  nref = nalt = 0;

  nref_nalt(&wrap->values[1], array);
  nref = *(long long int *)wrap->values[1];
  nalt = *(long long int *)wrap->values[2];
  
  if (nref + nalt == wrap->nsam) {
    if (nref < nalt) {
      *(long long int *)wrap->values[3 + nref] += 1;
    } else {
      *(long long int *)wrap->values[3 + nalt] += 1;
    }
    *(long long int *)wrap->values[0] += 1;
  }
}

void pq_swwrite_sfs(struct SWrap *wrap)
{
  int i;
  sprintf(wrap->outs[0], "%d", wrap->nsam);
  sprintf(wrap->outs[1], "%lli", *(long long int *)wrap->values[0]);
  for (i = 2; i < wrap->nouts; i++) {
    sprintf(wrap->outs[i], "%lli", *(long long int *)wrap->values[i + 1]);
  }
}

void pq_swclear_sfs(struct SWrap *wrap)
{
  int i;
  for (i = 0; i < wrap->nvalues; i++) {
    *(long long int *)wrap->values[i] = 0;
  }
}

// nvalues = [nvsites, nref, nalt, sfs1, sfs2, ...]
// nouts   = [nsam, nvsites, sfs1, sfs2, ...]
// sfs: number of occurances [0, 1, ..., n + 1]

void PQ_SFS_INIT(struct SWrap *wrap, int nsam)
{
  int i;
  wrap->nsam = nsam;
  wrap->nouts = 2 + ((nsam / 2) + 1);
  wrap->nvalues = 3 + ((nsam / 2) + 1);
  wrap->outs = calloc(wrap->nouts, sizeof(char *));
  wrap->values = calloc(wrap->nvalues, sizeof(void *));
  for (i = 0; i < wrap->nouts; i++) {
    wrap->outs[i] = calloc(128, sizeof(char));
  }
  for (i = 0; i < wrap->nvalues; i++) {
    wrap->values[i] = malloc(sizeof(long long int));
    *(long long int *)wrap->values[i] = 0;
  }
  wrap->update = pq_swupdate_sfs;
  wrap->write = pq_swwrite_sfs;
  wrap->clear = pq_swclear_sfs;
}
