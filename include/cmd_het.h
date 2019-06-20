#ifndef cmd_het_h__
#define cmd_het_h__

#define HET_FRMT 1
#define HET_NAME "het"
#define HET_DESC "Calculate heterozygosity"
#define HET_DEFS "-f 1 -c 1 -p 3 -k 5-%d -b 1 -d '\t'"

extern void init_het(struct SWrap *wrap, int nsam);

#endif
