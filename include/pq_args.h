#ifndef pq_htable_h__
#define pq_htable_h__
#endif

#ifndef pq_args_h__
#define pq_args_h__

extern int *KCOLS;
extern int NKARGS;
extern unsigned int CHROM;
extern unsigned int POS;
extern unsigned int FCOL;
extern struct HashTable ARGHASH;

extern void pq_update_args(unsigned int argc, char **argv);
extern void pq_init_args();
extern void pq_free_args();

#endif
