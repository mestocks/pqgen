#ifndef pq_generics_h__
#define pq_generics_h__


struct StatObject {
  int nsam;
  int nouts;
  int nvalues;
  char **outs;
  void **values;
  void (*update)(struct StatObject *, char **);
  void (*write)(struct StatObject *);
  void (*clear)(struct StatObject *);
};

extern void pq_swfree(struct StatObject *stats);


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
