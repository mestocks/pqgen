#include <stdio.h>
#include <stdlib.h>

#include "pq_args.h"
#include "pq_generics.h"
#include "pq_genetics.h"


void update_sfs(struct StatObject *stats, char **array)
{
  int nref, nalt;

  pq_gtstat(&stats->values[1], array);
  nref = *(int *)stats->values[2];
  nalt = *(int *)stats->values[3];
  
  if (nref + nalt == stats->nsam) {
    if (nref < nalt) {
      *(long long int *)stats->values[4 + nref] += 1;
    } else {
      *(long long int *)stats->values[4 + nalt] += 1;
    }
    *(long long int *)stats->values[0] += 1;
  }
}

void write_sfs(struct StatObject *stats)
{
  int i;
  sprintf(stats->outs[0], "%d", stats->nsam);
  sprintf(stats->outs[1], "%lli", *(long long int *)stats->values[0]);
  for (i = 2; i < stats->nouts; i++) {
    sprintf(stats->outs[i], "%lli", *(long long int *)stats->values[i + 2]);
  }
}

void clear_sfs(struct StatObject *stats)
{
  int i;
  *(long long int *)stats->values[0] = 0;
  for (i = 1; i < 4; i++) {
    *(int *)stats->values[i] = 0;
  }
  for (i = 4; i < stats->nvalues; i++) {
    *(long long int *)stats->values[i] = 0;
  }
}

// nvalues = [nvsites, nhet, nref, nalt, sfs1, sfs2, ...]
// nouts   = [nsam, nvsites, sfs1, sfs2, ...]
// sfs: number of occurances [0, 1, ..., n + 1]

void init_sfs(struct StatObject *stats, int nsam)
{
  int i;
  stats->nsam = nsam;
  stats->nouts = 2 + ((nsam / 2) + 1);
  stats->nvalues = 4 + ((nsam / 2) + 1);
  stats->outs = calloc(stats->nouts, sizeof(char *));
  stats->values = calloc(stats->nvalues, sizeof(void *));
  for (i = 0; i < stats->nouts; i++) {
    stats->outs[i] = calloc(128, sizeof(char));
  }
  
  stats->values[0] = malloc(sizeof(long long int));
  *(long long int *)stats->values[0] = 0;
  for (i = 1; i < 4; i++) {
    stats->values[i] = malloc(sizeof(int));
    *(int *)stats->values[i] = 0;
  }
  for (i = 4; i < stats->nvalues; i++) {
    stats->values[i] = malloc(sizeof(long long int));
    *(long long int *)stats->values[i] = 0;
  }
  stats->update = update_sfs;
  stats->write = write_sfs;
  stats->clear = clear_sfs;
}
