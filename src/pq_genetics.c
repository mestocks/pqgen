#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_htable.h>
#include <rwk_parse.h>
#include <pq_genetics.h>

unsigned int pq_alldna(const char *dna) {
  const char *tmp;
  unsigned int alldna;

  tmp = dna;
  alldna = 1;
  while (*tmp && alldna == 1) {
    switch (*tmp) {
    case 'A':
      break;
    case 'T':
      break;
    case 'G':
      break;
    case 'C':
      break;
    case 'a':
      break;
    case 't':
      break;
    case 'g':
      break;
    case 'c':
      break;
    default: alldna = 0;
    }
    tmp++;
  }
  return alldna;
}

void pq_complement(char *codon) {
  int bp;
  for (bp = 0; bp < 3; bp++) {
    switch (codon[bp]) {
    case 'A': codon[bp] = 'T'; break;
    case 'T': codon[bp] = 'A'; break;
    case 'G': codon[bp] = 'C'; break;
    case 'C': codon[bp] = 'G'; break;
    case 'a': codon[bp] = 't'; break;
    case 't': codon[bp] = 'a'; break;
    case 'g': codon[bp] = 'c'; break;
    case 'c': codon[bp] = 'g'; break;
    default: codon[bp] = codon[bp];
    }
  }
}

void pq_dna_upper(char *codon) {
  int bp;
  for (bp = 0; bp < 3; bp++) {
    switch (codon[bp]) {
    case 'A': codon[bp] = 'A'; break;
    case 'T': codon[bp] = 'T'; break;
    case 'G': codon[bp] = 'G'; break;
    case 'C': codon[bp] = 'C'; break;
    case 'a': codon[bp] = 'A'; break;
    case 't': codon[bp] = 'T'; break;
    case 'g': codon[bp] = 'G'; break;
    case 'c': codon[bp] = 'C'; break;
    default: codon[bp] = codon[bp];
    }
  }
}

void pq_reverse(char *codon) {
  char one  = codon[0];
  char three = codon[2];
  codon[0] = three;
  codon[2] = one;
}


struct rwkHashTable CODON_TO_NSYN;
struct rwkHashTable CODON_TO_AMINO;

void file2doubleHash(char *fname, int hsize) {

  FILE *fp;
  char *ptr;
  double *dptr;
  char **array;
  char buffer[2048];
  const char delim = ' ';
  
  fp = fopen(fname, "r");
  rwk_create_hash(&CODON_TO_NSYN, hsize);
  array = calloc(2, sizeof (char*));
  
  while (fgets(buffer, sizeof(buffer), fp)) {
    if (rwk_str2array(array, buffer, 2, &delim) == -1) {
      free(array);
      rwk_free_hash(&CODON_TO_NSYN);
      fclose(fp);
      exit(1);
    }
    ptr = malloc(2048 * sizeof (char));
    strcpy(ptr, array[0]);
    dptr = malloc(sizeof (double));
    *dptr = atof(array[1]);
    rwk_insert_hash(&CODON_TO_NSYN, ptr, dptr);
  }
  free(array);
  fclose(fp);
}

void file2charHash(char *fname, int hsize) {

  FILE *fp;
  char *ptr;
  char *cptr;
  char **array;
  char buffer[2048];
  const char delim = ' ';
  
  fp = fopen(fname, "r");
  rwk_create_hash(&CODON_TO_AMINO, hsize);
  array = calloc(2, sizeof (char*));
  
  while (fgets(buffer, sizeof(buffer), fp)) {
    if (rwk_str2array(array, buffer, 2, &delim) == -1) {
      free(array);
      rwk_free_hash(&CODON_TO_AMINO);
      fclose(fp);
      exit(1);
    }
    ptr = malloc(2048 * sizeof (char));
    strcpy(ptr, array[0]);
    cptr = malloc(2048 * sizeof (char));
    strcpy(cptr, array[1]);
    rwk_insert_hash(&CODON_TO_AMINO, ptr, cptr);
  }
  free(array);
  fclose(fp);
}
