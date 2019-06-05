#include <stdlib.h>
#include <string.h>

#include "pq_parse.h"

int pq_countcols(const char *buffer, const char *delim)
{
  int count;
  const char *tmp;
  count = 0;
  tmp = buffer;
  while (*tmp) {
    if (*delim == *tmp) {
      count++;
    }
    tmp++;
  }
  return count + 1;
}

/*
 *
 * pq_k2array
 * 
 * Convert a textual representation of a range of integers into an 
 * array of integers.
 *
 * int *array - a pointer to an array of integers of size nk. This 
 * will be populated with the range denoted by text.
 *
 * char *text - a pointer to a NULL terminated array of characters 
 * denoting the range of integers. Each integer is comma separated, 
 * but a range of integers can be denoted by a dash (and can be in 
 * descending order). For example, "1,2,3" "1-3" "1,4,7" "2,3-6,9,10" 
 * "1,9-2". 
 *
 * int nk - the number of integers given in the integer range. This 
 * can be given by pq_nkrange().
 *
 * Examples:
 *
 * "1,2,3"       -> [1, 2, 3]
 * "1-3"         -> [1, 2, 3]
 * "1,4,7"       -> [1, 4, 7]
 * "2,3-6,9,10"  -> [2, 3, 4, 5, 6, 9, 10]
 * "1,9-2"       -> [1, 9, 8, 7, 6, 5, 4, 3, 2]
 */

void pq_k2array(int *array, char *text, int nk)
{
  int a, i, j, c;
  int strsize;
  int cval, pval;
  char *tmp, *curr, *prev;

  strsize = strlen(text);
  curr = calloc(strsize + 1, sizeof *curr);
  prev = calloc(strsize + 1, sizeof *curr);

  tmp = text;
  prev[0] = '\0';
  a = c = i = j = 0;
  while (*tmp) {
    if (*tmp == ',') {
      curr[i] = prev[j] = '\0';
      cval = atoi(curr);
      pval = atoi(prev);
      c = cval >= pval ? 1: -1;
      while (pval != cval) {
	array[a] = pval;
	a++;
	pval = pval + c;
      }
      array[a] = pval;
      i = j = 0;
      a++;
      tmp++;
      curr[0] = prev[0] = '\0';
    } else if (*tmp == '-') {
      prev[j] = curr[0] = '\0';
      i = 0;
      j++;
      tmp++;
    } else {
      curr[i] = prev[j] = *tmp;
      i++;
      j++;
      tmp++;
    }
  }
  curr[i] = prev[j] = '\0';
  cval = atoi(curr);
  pval = atoi(prev);
  c = cval >= pval ? 1: -1;
  while (pval != cval) {
    array[a] = pval;
    a++;
    pval = pval + c;
  }
  array[a] = pval;
  
  free(curr);
  free(prev);
}


int pq_nkrange(char *text)
{
  int a, i, j, c;
  int strsize;
  int cval, pval;
  char *tmp, *curr, *prev;

  if (*text == '\0') {
    return 0;
  }

  strsize = strlen(text);
  curr = calloc(strsize + 1, sizeof *curr);
  prev = calloc(strsize + 1, sizeof *prev);

  tmp = text;
  prev[0] = '\0';
  a = c = i = j = 0;
  while (*tmp) {
    if (*tmp == ',') {
      curr[i] = prev[j] = '\0';
      cval = atoi(curr);
      pval = atoi(prev);
      c = cval >= pval ? 1: -1;
      while (pval != cval) {
	a++;
	pval = pval + c;
      }
      i = j = 0;
      a++;
      tmp++;
      curr[0] = prev[0] = '\0';
    } else if (*tmp == '-') {
      prev[j] = curr[0] = '\0';
      i = 0;
      j++;
      tmp++;
    } else {
      curr[i] = prev[j] = *tmp;
      i++;
      j++;
      tmp++;
    }
  }
  curr[i] = prev[j] = '\0';
  cval = atoi(curr);
  pval = atoi(prev);
  c = cval >= pval ? 1: -1;
  while (pval != cval) {
    a++;
    pval = pval + c;
  }
  a++;

  free(curr);
  free(prev);
  
  return a;
}

/*
 * pq_str2array - replace delim with '/0' and insert (char *) into array
 *
 *               +-    +-        +-
 * buffer: ["str1\0str2\0... strn\0"] 
 *           ^     ^         ^
 * array:  [ ptr*, ptr*, ... ptr* ]
 *
 * Note that buffer is not const (delim is replaced with '/0'). No memory is 
 * allocated (char pointers from buffer are used in array). This also 
 * means that the pointers stored in array are freed when buffer is freed.
 *
 * if (ndel + 1 > ncols) returns -1
 * if (ndel <= ncols) NULL pointer added to array indices ndel to ncols - 1 & returns 0
 *
 * Example usage:
 *
 * char **array = calloc(ncols, sizeof (char *));
 * if (pq_str2array(array, buffer, ncols, &delim) == -1) {
 *   free(array);
 *   // free buffer if necessary
 *   exit(1);
 * } else {
 *   // do something with array
 * }
 * free(array);
 * // free buffer if necessary
 */
int pq_str2array(char **array, char *buffer, const int ncols, const char *delim)
{
  int coln;
  char *tmp;
  char newline;

  coln = 0;
  tmp = buffer;
  newline = '\n';
  
  array[coln] = tmp;
  while (*tmp && newline != *tmp) {
    
    if (*delim == *tmp) {
      *tmp = '\0';
      coln++;
      tmp++;
      if (coln == ncols) {
	return -1;
      } 
      array[coln] = tmp;
    } else {
      tmp++;
    }
  }
  *tmp = '\0';
  
  while (coln < ncols - 1) {
    coln++;
    array[coln] = NULL;
  }
  
  return 0;
}



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
  pq_str2array(row->array, buffer, row->ncols, delim);
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
