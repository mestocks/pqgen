#ifndef pq_theta_h__
#define pq_theta_h__

#define PQ_THETA_FRMT 2
#define PQ_THETA_NAME "theta"
#define PQ_THETA_DESC "Calculate site frequency based statistics"

extern void nref_nalt(void **counts, char **array);
extern void PQ_THETA_INIT(struct SWrap *wrap, int nsam);

#endif
