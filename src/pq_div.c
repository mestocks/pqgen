#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pq_htable.h"
#include "pq_generics.h"
#include "pq_args.h"

void pq_swupdate_div(struct StatObject *stats, char **array)
{

  char *ref;
  char *alt;

  ref = array[KCOLS[0]];
  alt = array[KCOLS[1]];

  if (strcmp(ref, ".") != 0 && strcmp(alt, ".") != 0) {
    if (strcmp(ref, alt) != 0) {
      *(long long int *)stats->values[0]+=1;
    }
    *(long long int *)stats->values[1]+=1;
  }
}

void pq_swwrite_div(struct StatObject *stats)
{
  long long int diffs, nvsites;
  diffs =  *(long long int *)stats->values[0];
  nvsites =  *(long long int *)stats->values[1];
  sprintf(stats->outs[0], "%lli", nvsites);
  if (strcmp((char *)pq_lookup_hash(&ARGHASH, "-b"), "0") == 0) {
    sprintf(stats->outs[1], "%f", (double)diffs);
  } else {
    sprintf(stats->outs[1], "%f", (double)diffs / (double)nvsites);
  }
}

void pq_swclear_div(struct StatObject *stats)
{
  int i;
  for (i = 0; i < stats->nvalues; i++) {
    *(long long int *)stats->values[i] = 0;
  }
}

void PQ_DIV_INIT(struct StatObject *stats, int nsam)
{
  // nvsites, diffs
  int i;
  stats->nsam = nsam;
  stats->nouts = 2;
  stats->nvalues = 2;
  stats->outs = calloc(stats->nouts, sizeof(char *));
  stats->values = calloc(stats->nvalues, sizeof(void *));
  for (i = 0; i < stats->nouts; i++) {
    stats->outs[i] = calloc(128, sizeof(char));
  }
  for (i = 0; i < stats->nvalues; i++) {
    stats->values[i] = malloc(sizeof(long long int));
    *(long long int *)stats->values[i] = 0;
  }
  stats->update = pq_swupdate_div;
  stats->write = pq_swwrite_div;
  stats->clear = pq_swclear_div;
}
