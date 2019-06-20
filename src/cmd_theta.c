#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pq_htable.h"

#include "pq_args.h"
#include "pq_generics.h"
#include "pq_genetics.h"
#include "pq_sfstats.h"


void update_theta(struct StatObject *stats, char **array)
{
  int refalt[2];

  refalt[0] = refalt[1] = 0;
  count_alleles_from_genotypes(refalt, array);
  //  gt_to_loopalt(refalt, array);

  // Missing data?
  if (refalt[0] + refalt[1] == stats->nsam) {
    if (refalt[1] > 0 && refalt[0] > 0) {
      *(long long int *)stats->values[1] += 1;
    }
    *(long long int *)stats->values[0] += 1;
    *(long long int *)stats->values[2] += PairwiseDiffs(stats->nsam, refalt[0]);
  }
}


void write_theta(struct StatObject *stats)
{
  long long int s, pisum, nvsites;
  double tw, combs, pi, tajd;

  nvsites = *(long long int *)stats->values[0];
  s =  *(long long int *)stats->values[1];
  pisum =  *(long long int *)stats->values[2];
  tw = WattersonsTheta(stats->nsam, s);
  combs = PairwiseCombs(stats->nsam);
  pi = TajimasTheta(combs, pisum);
  tajd = TajimasD(stats->nsam, s, tw, pi);

  sprintf(stats->outs[0], "%d", stats->nsam);
  sprintf(stats->outs[1], "%lli", nvsites);
  sprintf(stats->outs[2], "%lli", s);
  if (strcmp((char *)pq_lookup_hash(&ARGHASH, "-b"), "0") == 0) {
    sprintf(stats->outs[3], "%f", tw);
    sprintf(stats->outs[4], "%f", pi);
  } else {
    sprintf(stats->outs[3], "%f", tw / nvsites);
    sprintf(stats->outs[4], "%f", pi / nvsites);
  }
  sprintf(stats->outs[5], "%f", tajd);
}


void clear_theta(struct StatObject *stats)
{
  int i;
  for (i = 0; i < 3; i++) {
    *(long long int *)stats->values[i] = 0;
  }
    for (i = 3; i < stats->nvalues; i++) {
    *(int *)stats->values[i] = 0;
  }
}


// values = [nvsites, s, pisum, nhet, nref, nalt]
// out    = [nsam, nvsites, s, tw, pi, tajd]

void init_theta(struct StatObject *stats, int nsam)
{
  int i;
  stats->nsam = nsam;
  stats->nouts = 6;
  stats->nvalues = 6;
  stats->outs = calloc(stats->nouts, sizeof(char *));
  stats->values = calloc(stats->nvalues, sizeof(void *));
  for (i = 0; i < stats->nouts; i++) {
    stats->outs[i] = calloc(128, sizeof(char));
  }
  
  for (i = 0; i < 3; i++) {
    stats->values[i] = malloc(sizeof(long long int));
    *(long long int *)stats->values[i] = 0;
  }
  for (i = 3; i < stats->nvalues; i++) {
    stats->values[i] = malloc(sizeof(int));
    *(int *)stats->values[i] = 0;
  }
  
  stats->update = update_theta;
  stats->write = write_theta;
  stats->clear = clear_theta;
}

