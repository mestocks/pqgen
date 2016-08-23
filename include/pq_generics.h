#ifndef pq_generics_h__
#define pq_generics_h__

struct pq_parameters {
  unsigned int CHROM;
  unsigned int POS;
  unsigned int FCOL;
  int *KCOLS;
  int nkargs;
  struct rwkHashTable arghash;
  void (*update)(struct pq_parameters *ptr, unsigned int argc, char **argv);
};

extern void pq_init_parameters(struct pq_parameters *params);
extern void pq_free_parameters(struct pq_parameters *params);
extern void pq_update_parameters(struct pq_parameters *params, unsigned int argc, char **argv);

struct Data {
  char *chr;
  unsigned long long int *start;
  unsigned long long int *end;
  char *name;
  unsigned long long int *nvsites;
  char *value;
};

extern void pq_print_bed(struct Data *data);
extern void pq_print_csv(struct Data *data);

#endif
