#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double pqWattersonsTheta(int nsam, int s) {
  int i;
  double a1, tw;

  if (s == 0) {
    tw = 0.0; }
  else {
    a1 = 0.0;
    i = nsam - 1;
    while (i >= 1) {
      a1 += 1.0 / i;
      i--;
    }
    tw = s / a1;
  }
  return tw;
}

long long int pqPairwiseDiffs(int nsam, int nminor) {
  long long int diffs = nminor * (nsam - nminor);
  return diffs;
}

double pqPairwiseCombs(int nsam) {
  // n choose 2: (n - (k - i)) / i
  // i = 1: n - (2 - 1) / 1 = n - 1
  // i = 2: n - (2 - 2) / 2 = n / 2
  double combs = (nsam - 1) * (nsam / 2);
  return combs;
}

double pqTajimasTheta(double combs, long long int diffs) {
  return diffs / (double)combs;
}
