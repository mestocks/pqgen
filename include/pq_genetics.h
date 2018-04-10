#ifndef pq_htable_h__
#define pq_htable_h__
#endif

#ifndef pq_genetics_h__
#define pq_genetics_h__

extern void pq_complement(char *codon);
extern void pq_dna_upper(char *codon);
extern void pq_reverse(char *codon);
extern unsigned int pq_alldna(const char *dna);

extern void pq_gtstat(void **info, char **array);

extern struct HashTable CODON_TO_NSYN;
extern struct HashTable CODON_TO_AMINO;

extern void file2charHash(char *fname, int hsize);
extern void file2doubleHash(char *fname, int hsize);

#endif
