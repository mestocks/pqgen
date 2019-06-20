#ifndef cmd_sfs_h__
#define cmd_sfs_h__

#define PQ_SFS_FRMT 2
#define PQ_SFS_NAME "sfs"
#define PQ_SFS_DESC "Output the folded site frequency spectrum"
#define PQ_SFS_DEFS "-f 1 -c 1 -p 3 -k 5-%d -d '\t'"

extern void pq_sfs_init(struct SWrap *wrap, int nsam);

#endif
