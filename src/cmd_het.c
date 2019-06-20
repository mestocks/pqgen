#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pq_htable.h"

#include "pq_args.h"
#include "pq_generics.h"
#include "pq_genetics.h"


void update_het(struct StatObject *stats, char **array)
{
  int nhet, nref, nalt;
  double het;
  
  pq_gtstat(&stats->values[2], array);
  nhet = *(int *)stats->values[2];
  nref = *(int *)stats->values[3];
  nalt = *(int *)stats->values[4];

  if (nref + nalt > 0) {
    het = (double)nhet / (((double)nref + (double)nalt) / 2);
    *(double *)stats->values[1] += het;
    *(long long int *)stats->values[0] += 1;
  }
}


void write_het(struct StatObject *stats)
{
  double het;
  sprintf(stats->outs[0], "%d", stats->nsam);
  sprintf(stats->outs[1], "%lli", *(long long int *)stats->values[0]);

  if (strcmp((char *)pq_lookup_hash(&ARGHASH, "-b"), "0") == 0) {
    het = *(double *)stats->values[1];
  } else {
    het = *(double *)stats->values[1] / *(long long int *)stats->values[0];
  }
  sprintf(stats->outs[2], "%f", het);
}


void clear_het(struct StatObject *stats)
{
  *(long long int *)stats->values[0] = 0;
  *(double *)stats->values[1] = 0.0;
}


// nvalues = [nvsites, het, nhet, nref, nalt]
// nouts   = [nsam, nvsites, uhet]

void init_het(struct StatObject *stats, int nsam)
{
  int i;
  stats->nsam = nsam;
  stats->nouts = 3;
  stats->nvalues = 5;
  stats->outs = calloc(stats->nouts, sizeof(char *));
  stats->values = calloc(stats->nvalues, sizeof(void *));
  for (i = 0; i < stats->nouts; i++) {
    stats->outs[i] = calloc(128, sizeof(char));
  }
  
  stats->values[0] = malloc(sizeof(long long int));
  *(long long int *)stats->values[0] = 0;

  stats->values[1] = malloc(sizeof(double));
  *(double *)stats->values[1] = 0.0;
  
  for (i = 2; i < stats->nvalues; i++) {
    stats->values[i] = malloc(sizeof(int));
    *(int *)stats->values[i] = 0;
  }
  
  stats->update = update_het;
  stats->write = write_het;
  stats->clear = clear_het;
}

