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


struct SWrap {
  int nsam;
  int nouts;
  int nvalues;
  char **outs;
  void **values;
  void (*update)(struct SWrap *, char **);
  void (*write)(struct SWrap *);
  void (*clear)(struct SWrap *);
};

extern void pq_swfree(struct SWrap *wrap);


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
