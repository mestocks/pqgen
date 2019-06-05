#include <stdio.h>
#include <stdlib.h>

#include "pq_args.h"
#include "pq_generics.h"
#include "pq_genetics.h"


void pq_swupdate_sfs(struct SWrap *wrap, char **array)
{
  int nref, nalt;

  pq_gtstat(&wrap->values[1], array);
  nref = *(int *)wrap->values[2];
  nalt = *(int *)wrap->values[3];
  
  if (nref + nalt == wrap->nsam) {
    if (nref < nalt) {
      *(long long int *)wrap->values[4 + nref] += 1;
    } else {
      *(long long int *)wrap->values[4 + nalt] += 1;
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
    sprintf(wrap->outs[i], "%lli", *(long long int *)wrap->values[i + 2]);
  }
}

void pq_swclear_sfs(struct SWrap *wrap)
{
  int i;
  *(long long int *)wrap->values[0] = 0;
  for (i = 1; i < 4; i++) {
    *(int *)wrap->values[i] = 0;
  }
  for (i = 4; i < wrap->nvalues; i++) {
    *(long long int *)wrap->values[i] = 0;
  }
}

// nvalues = [nvsites, nhet, nref, nalt, sfs1, sfs2, ...]
// nouts   = [nsam, nvsites, sfs1, sfs2, ...]
// sfs: number of occurances [0, 1, ..., n + 1]

void pq_sfs_init(struct SWrap *wrap, int nsam)
{
  int i;
  wrap->nsam = nsam;
  wrap->nouts = 2 + ((nsam / 2) + 1);
  wrap->nvalues = 4 + ((nsam / 2) + 1);
  wrap->outs = calloc(wrap->nouts, sizeof(char *));
  wrap->values = calloc(wrap->nvalues, sizeof(void *));
  for (i = 0; i < wrap->nouts; i++) {
    wrap->outs[i] = calloc(128, sizeof(char));
  }
  
  wrap->values[0] = malloc(sizeof(long long int));
  *(long long int *)wrap->values[0] = 0;
  for (i = 1; i < 4; i++) {
    wrap->values[i] = malloc(sizeof(int));
    *(int *)wrap->values[i] = 0;
  }
  for (i = 4; i < wrap->nvalues; i++) {
    wrap->values[i] = malloc(sizeof(long long int));
    *(long long int *)wrap->values[i] = 0;
  }
  wrap->update = pq_swupdate_sfs;
  wrap->write = pq_swwrite_sfs;
  wrap->clear = pq_swclear_sfs;
}
