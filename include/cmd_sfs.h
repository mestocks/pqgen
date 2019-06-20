#ifndef cmd_sfs_h__
#define cmd_sfs_h__

#define SFS_FRMT 2
#define SFS_NAME "sfs"
#define SFS_DESC "Output the folded site frequency spectrum"
#define SFS_DEFS "-f 1 -c 1 -p 3 -k 5-%d -d '\t'"

extern void init_sfs(struct SWrap *wrap, int nsam);

#endif
