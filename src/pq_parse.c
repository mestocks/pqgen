#include <stdlib.h>
#include <pq_parse.h>
#include <rwk_parse.h>


char *fetch_chrom(struct GenericRow *row)
{
  return row->array[row->chr_col];
}

char *fetch_fcol(struct GenericRow *row)
{
  return row->array[row->fac_col];
}

long long int fetch_pos(struct GenericRow *row)
{
  return atoll(row->array[row->pos_col]);
}

void update_row(struct GenericRow *row, char *buffer, const char *delim)
{
  rwk_str2array(row->array, buffer, row->ncols, delim);
}

void init_row(struct GenericRow *row, int ncols, int chr_col, int pos_col, int fac_col)
{
  row->ncols = ncols;
  row->array = calloc(row->ncols, sizeof (char*));
  row->chr_col = chr_col;
  row->pos_col = pos_col;
  row->fac_col = fac_col;
  row->pos = fetch_pos;
  row->chrom = fetch_chrom;
  row->factor = fetch_fcol;
  row->update = update_row;
}

void free_row(struct GenericRow *row)
{
  free(row->array);
}
