#ifndef pq_genetics_h__
#define pq_genetics_h__

extern unsigned int pq_alldna(const char *dna);
extern void pq_complement(char *codon);
extern void pq_dna_upper(char *codon);
extern void pq_reverse(char *codon);

#endif
