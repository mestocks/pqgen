#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double WattersonsTheta(int nsam, int s) {
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

long long int PairwiseDiffs(int nsam, int nminor) {
  long long int diffs = nminor * (nsam - nminor);
  return diffs;
}

double PairwiseCombs(int nsam) {
  // n choose 2: (n - (k - i)) / i
  // i = 1: n - (2 - 1) / 1 = n - 1
  // i = 2: n - (2 - 2) / 2 = n / 2
  double combs = (nsam - 1) * (nsam / 2);
  return combs;
}

double TajimasTheta(double combs, long long int diffs) {
  return diffs / (double)combs;
}

double TajimasD(int nsam, long long int s, double tw, double pi) {
  double a1, a2, b1, b2, c1, c2, e1, e2, V;
  
  a1 = 0;
  a2 = 0;
  for (int i = 1; i < nsam; i++) {
    a1 += 1.0 / i;
    a2 += 1.0 / (i * i);
  }
  
  b1 = (nsam + 1.0) / (3.0 * (nsam - 1.0));
  b2 = 2.0 * (nsam * nsam + nsam + 3.0) / (9.0 * nsam * (nsam - 1.0));
  c1 = b1 - 1.0 / a1;
  c2 = b2 - (nsam + 2.0) / (a1 * nsam) + a2 / (a1 * a1);
  e1 = c1 / a1;
  e2 = c2 / (a1 * a1 + a2);
  V = e1 * s + e2 * s * (s - 1.0);
  
  return (pi - tw) / sqrt(V);
}
