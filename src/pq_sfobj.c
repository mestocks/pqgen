#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "pq_sfobj.h"
#include "pq_sfstats.h"

void pqThetaWAdd(struct pqThetaW *ptr, int s) {
  ptr->s = ptr->s + s;
  ptr->nsites++;
}

double pqThetaWEval(struct pqThetaW *ptr) {
  double tw = pqWattersonsTheta(ptr->nsam, ptr->s);
  return tw;
}

void pqThetaWInit(struct pqThetaW *ptr, int nsam) {
  ptr->s = 0;
  ptr->nsites = 0;
  ptr->nsam = nsam;
  ptr->add = pqThetaWAdd;
  ptr->eval = pqThetaWEval;
  ptr->reset = pqThetaWReset;
}

void pqThetaWReset(struct pqThetaW *ptr) {
  ptr->s = 0;
  ptr->nsites = 0;
}


// Tajima's theta

void pqThetaPiAdd(struct pqThetaPi *ptr, int nminor) {
  ptr->pisum += pqPairwiseDiffs(ptr->nsam, nminor);
  ptr->nsites++;
}

double pqThetaPiEval(struct pqThetaPi *ptr) {
  // n choose 2: (n - (k - i)) / i
  // i = 1: n - (2 - 1) / 1 = n - 1
  // i = 2: n - (2 - 2) / 2 = n / 2
  double combs = pqPairwiseCombs(ptr->nsam);
  return pqTajimasTheta(combs, ptr->pisum);
}

void pqThetaPiInit(struct pqThetaPi *ptr, int nsam) {
  ptr->pisum = 0;
  ptr->nsites = 0;
  ptr->nsam = nsam;
  ptr->add = pqThetaPiAdd;
  ptr->eval = pqThetaPiEval;
  ptr->reset = pqThetaPiReset;
}

void pqThetaPiReset(struct pqThetaPi *ptr) {
  ptr->pisum = 0;
  ptr->nsites = 0;
}

// Folded site frequency spectrum

int pqSizeofFoldedSFS(int nsam) {
  if (nsam % 2 == 0) {
    return nsam / 2;
  } else {
    return (nsam - 1) / 2;
  }
}

//int pqFoldAllele(int nsam, int nminor) {
  
  /*
    n b x y
    3 1 1 0 [0, 1, 2, 3] -> [0]
    3 1 2 0
    4 2 1 0 [0, 1, 2, 3, 4] -> [0, 1]
    4 2 2 1
    4 2 3 0
    5       [0, 1, 2, 3, 4, 5] -> [0, 1]
   */
//}

//void pqFoldedSFS(struct pqFoldedSFS *ptr, int nsam) {
  
//}
