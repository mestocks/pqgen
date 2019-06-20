#ifndef pq_div_h__
#define pq_div_h__

#define PQ_DIV_FRMT 1
#define PQ_DIV_NAME "div"
#define PQ_DIV_DESC "Calculate divergence based statistics"
#define PQ_DIV_DEFS "-f 1 -c 1 -p 3 -k 5,6 -b 1 -d '\t'"

extern void PQ_DIV_INIT(struct StatObject *stats, int nsam);

#endif
