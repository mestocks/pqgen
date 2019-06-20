#ifndef pq_pnds_h__
#define pq_pnds_h__

#define PQ_PNDS_FRMT 1
#define PQ_PNDS_NAME "pnds"
#define PQ_PNDS_DESC "Count the number of silent and replacement substitutions and polymorphisms"
#define PQ_PNDS_DEFS "-f 1 -c 1 -p 3 -k 7-%d -d '\t'"

extern void PQ_PNDS_INIT(struct StatObject *stats, int nsam);

#endif
