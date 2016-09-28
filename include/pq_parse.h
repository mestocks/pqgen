#ifndef pq_parse_h__
#define pq_parse_h__

struct GenericRow {
  char **array;
  int row_index;
  int ncols;
  int chr_col;
  int pos_col;
  int fac_col;
  char *(*chrom)(struct GenericRow *);
  char *(*factor)(struct GenericRow *);
  long long int (*pos)(struct GenericRow *);
  void (*update)(struct GenericRow *, char *, const char *);
};

extern void init_row(struct GenericRow *row, int ncols, int chr_col, int pos_col, int fac_col);
extern void free_row(struct GenericRow *row);

#endif
