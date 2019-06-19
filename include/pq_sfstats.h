#ifndef pq_sfstats_h__
#define pq_sfstats_h__

double WattersonsTheta(int nsam, int s);
long long int PairwiseDiffs(int nsam, int nminor);
double PairwiseCombs(int nsam);
double TajimasTheta(double combs, long long int diffs);
double TajimasD(int nsam, long long int s, double tw, double pi);

#endif
