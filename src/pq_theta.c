#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pq_htable.h>

#include <pq_args.h>
#include <pq_generics.h>
#include <pq_genetics.h>
#include <pq_sfstats.h>


void pq_swupdate_theta(struct SWrap *wrap, char **array)
{
  int nref, nalt;

  pq_gtstat(&wrap->values[3], array);
  nref = *(int *)wrap->values[4];
  nalt = *(int *)wrap->values[5];
    
  if (nref + nalt == wrap->nsam) {
    if (nalt > 0 && nref > 0) {
      *(long long int *)wrap->values[1] += 1;
    }
    *(long long int *)wrap->values[2] += pqPairwiseDiffs(wrap->nsam, nref);
    //*(long long int *)wrap->values[2] += nref * (wrap->nsam - nref);
    *(long long int *)wrap->values[0] += 1;
  }
}


void pq_swwrite_theta(struct SWrap *wrap)
{
  long long int s, pisum, nvsites;
  double tw, combs, pi, tajd;

  nvsites = *(long long int *)wrap->values[0];
  s =  *(long long int *)wrap->values[1];
  pisum =  *(long long int *)wrap->values[2];
  tw = pqWattersonsTheta(wrap->nsam, s);
  combs = pqPairwiseCombs(wrap->nsam);
  pi = pqTajimasTheta(combs, pisum);
  tajd = pqTajimasD(wrap->nsam, s, tw, pi);

  sprintf(wrap->outs[0], "%d", wrap->nsam);
  sprintf(wrap->outs[1], "%lli", nvsites);
  sprintf(wrap->outs[2], "%lli", s);
  if (strcmp((char *)pq_lookup_hash(&ARGHASH, "-b"), "0") == 0) {
    sprintf(wrap->outs[3], "%f", tw);
    sprintf(wrap->outs[4], "%f", pi);
  } else {
    sprintf(wrap->outs[3], "%f", tw / nvsites);
    sprintf(wrap->outs[4], "%f", pi / nvsites);
  }
  sprintf(wrap->outs[5], "%f", tajd);
}


void pq_swclear_theta(struct SWrap *wrap)
{
  int i;
  for (i = 0; i < 3; i++) {
    *(long long int *)wrap->values[i] = 0;
  }
    for (i = 3; i < wrap->nvalues; i++) {
    *(int *)wrap->values[i] = 0;
  }
}


// values = [nvsites, s, pisum, nhet, nref, nalt]
// out    = [nsam, nvsites, s, tw, pi, tajd]

void pq_theta_init(struct SWrap *wrap, int nsam)
{
  int i;
  wrap->nsam = nsam;
  wrap->nouts = 6;
  wrap->nvalues = 6;
  wrap->outs = calloc(wrap->nouts, sizeof(char *));
  wrap->values = calloc(wrap->nvalues, sizeof(void *));
  for (i = 0; i < wrap->nouts; i++) {
    wrap->outs[i] = calloc(128, sizeof(char));
  }
  
  for (i = 0; i < 3; i++) {
    wrap->values[i] = malloc(sizeof(long long int));
    *(long long int *)wrap->values[i] = 0;
  }
  for (i = 3; i < wrap->nvalues; i++) {
    wrap->values[i] = malloc(sizeof(int));
    *(int *)wrap->values[i] = 0;
  }
  
  wrap->update = pq_swupdate_theta;
  wrap->write = pq_swwrite_theta;
  wrap->clear = pq_swclear_theta;
}

