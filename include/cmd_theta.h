#ifndef cmd_theta_h__
#define cmd_theta_h__

#define THETA_FRMT 2
#define THETA_NAME "theta"
#define THETA_DESC "Calculate site frequency based statistics"
#define THETA_DEFS "-f 1 -c 1 -p 3 -k 5-%d -b 1 -d '\t'"

extern void init_theta(struct SWrap *wrap, int nsam);

#endif
