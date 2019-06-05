#ifndef cmd_theta_h__
#define cmd_theta_h__

#define PQ_THETA_FRMT 2
#define PQ_THETA_NAME "theta"
#define PQ_THETA_DESC "Calculate site frequency based statistics"
#define PQ_THETA_DEFS "-f 1 -c 1 -p 3 -k 5-%d -b 1 -d '\t'"

extern void pq_theta_init(struct SWrap *wrap, int nsam);

#endif
