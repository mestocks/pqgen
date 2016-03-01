#ifndef pq_sfstats_h__
#define pq_sfstats_h__

double pqWattersonsTheta(int nsam, int s);
long long int pqPairwiseDiffs(int nsam, int nminor);
double pqPairwiseCombs(int nsam);
double pqTajimasTheta(double combs, long long int diffs);

#endif
