#ifndef pq_generics_h__
#define pq_generics_h__

extern int genotype_format(char *str);

struct pq_parameters {
  unsigned int CHROM;
  unsigned int POS;
  unsigned int FCOL;
  int *KCOLS;
  int nkargs;
  struct rwkHashTable arghash;
  void (*update)(struct pq_parameters *ptr, unsigned int argc, char **argv);
};

extern int kcols_to_segs(char **array, int *kcol, int nkargs);

extern void pq_init_parameters(struct pq_parameters *params);
extern void pq_free_parameters(struct pq_parameters *params);
extern void pq_update_parameters(struct pq_parameters *params, unsigned int argc, char **argv);


struct SWrap {
  int nsam;
  int nouts;
  int nvalues;
  char **outs;
  void **values;
  void (*update)(struct SWrap *, char **, struct pq_parameters *);
  void (*write)(struct SWrap *, struct pq_parameters *);
  void (*clear)(struct SWrap *);
};

extern void tgt_checkit(void **counts, char **array, struct pq_parameters *params);

extern void pq_swfree(struct SWrap *wrap);

extern void pq_swupdate_theta(struct SWrap *wrap, char **array, struct pq_parameters *params);
extern void pq_swwrite_theta(struct SWrap *wrap, struct pq_parameters *params);
extern void pq_swclear_theta(struct SWrap *wrap);
extern void pq_swinit_theta(struct SWrap *wrap, int nsam);

extern void pq_swupdate_div(struct SWrap *wrap, char **array, struct pq_parameters *params);
extern void pq_swwrite_div(struct SWrap *wrap, struct pq_parameters *params);
extern void pq_swclear_div(struct SWrap *wrap);
extern void pq_swinit_div(struct SWrap *wrap, int nsam);


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
