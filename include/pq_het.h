#ifndef pq_het_h__
#define pq_het_h__

#define PQ_HET_FRMT 1
#define PQ_HET_NAME "het"
#define PQ_HET_DESC "Calculate heterozygosity"
#define PQ_HET_DEFS "-f 1 -c 1 -p 3 -k 5-%d -b 1 -d '\t'"

extern void pq_het_init(struct SWrap *wrap, int nsam);

#endif
