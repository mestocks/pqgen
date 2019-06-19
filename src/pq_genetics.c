#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pq_parse.h"

#include "pq_htable.h"

#include "pq_args.h"
#include "pq_genetics.h"

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


// ~/~ ~/~
// ~|~ ~|~

// nhets

//    0 1 2
// 0 [     ]
// 1 [     ]
// 2 [     ]

// [00, 01, 02, 10, 11, 12, 20, 21, 22]
// [ 0,  1,  2,  3,  4,  5,  6,  7,  8]

// [x, y] -> (x * n) + y

void pq_gtstat(void **info, char **array)
{
  int i;
  int j;
  int a1;
  int a2;
  int nhet;
  int nref;
  int nalt;
  int nperm;
  char *curr;
  int *genos;
  int cgenos;
  int ngenos;
  int *alleles;
  int maxgenos;

  /* 
   * When maxgenos = 10 the assumption is therefore that
   * there are no alleles > 9 (i.e. length(~/~) == 3 &&
   * ~ can be treated as a single char).
   */

  maxgenos = 10;
  nperm = maxgenos * maxgenos;
  genos = calloc(nperm, sizeof(int));
  alleles = calloc(maxgenos, sizeof(int));

  /*
   * Add genotypes to array. 
   *
   * The c standard requires that the character set for 
   * digits is in numerical order so their numerical value 
   * can be obtained simply by subtracting '0'.
  */
  
  ngenos = 0;
  for (i = 0; i < NKARGS; i++) {
    curr = array[KCOLS[i]];
    if (curr[0] != '.') {
      a1 = curr[0] - '0';
      a2 = curr[2] - '0';
      genos[(a1 * maxgenos) + a2]++;
      alleles[a1]++;
      alleles[a2]++;
      ngenos++;
    }
  }
  
  nhet = 0;
  cgenos = 0;

  for (i = 0; i < maxgenos && cgenos != ngenos; i++) {
    for (j = 0; j < maxgenos && cgenos != ngenos; j++) {
      cgenos += genos[(i * maxgenos) + j];
      if (i != j) {
	nhet += genos[(i * maxgenos) + j];
      }
    }
  }

  nref = alleles[0];
  nalt = (2 * ngenos) - nref;

  *(int *)info[0] = nhet;
  *(int *)info[1] = nref;
  *(int *)info[2] = nalt;

  free(alleles);
  free(genos);
}

/* 
 * Assumes ref is 0 and anything non-zero is alt. Therefore,
 * different non-zero alleles are considered the same and will
 * not appear as a segregating site if the ref allele is not present.
 * For example:
 *    0/0 0/0   [4, 0]   no segregating site
 *    0/0 0/1   [3. 1]   segregating site
 *    0/0 0/2   [3. 1]   segregating site
 *    1/1 1/1   [0, 4]   no segregating site
 *    1/2 1/2   [0, 4]   no segregating site
 */
void count_refalt_from_genotypes(int *refalt, char **array)
{
  char *curr;
  int i, a1, a2;
  
  for (i = 0; i < NKARGS; i++) {
    curr = array[KCOLS[i]];
    if (curr[0] != '.') {
      a1 = curr[0] - '0';
      a2 = curr[2] - '0';
      refalt[a1 == 0]++;
      refalt[a2 == 0]++;
    }
  }
}

/*
 * Assumes there are only 2 alleles and uses the first allele that 
 * is encountered as the ref allele. Therefore, if there are more 
 * than 2 alleles then all non-ref alleles will be lumped together.
 * For example:
 *    0/0 0/0   [4, 0]
 *    0/0 0/1   [3, 1]
 *    1/2 1/2   [2, 2]
 *    1/2 2/3   [1, 3]
 */
void count_alleles_from_genotypes(int *refalt, char **array)
{
  char *curr;
  int i, a1, a2, is_first, ref;

  ref = is_first = 1;
  for (i = 0; i < NKARGS; i++) {
    curr = array[KCOLS[i]];
    if (curr[0] != '.') {
      a1 = curr[0] - '0';
      a2 = curr[2] - '0';
      if (is_first) {
	ref = a1;
	is_first = 0;
      }
      refalt[a1 == ref]++;
      refalt[a2 == ref]++;
    }
  }
}

void gt_to_loopalt(int *refalt, char **array)
{
  char *curr;
  int i, a1, a2, is_first, ref;

  size_t n;

  n = 0;
  while (array[KCOLS[n]][0] == '.' && n < NKARGS) {
    n++;
  }
  
  if (n < NKARGS) {
    ref = array[KCOLS[n]][0] - '0';
    for (i = n; i < NKARGS; i++) {
      curr = array[KCOLS[i]];
      if (curr[0] != '.') {
	a1 = curr[0] - '0';
	a2 = curr[2] - '0';
	refalt[a1 == ref]++;
	refalt[a2 == ref]++;
      }
    }
  }
}

struct HashTable CODON_TO_NSYN;
struct HashTable CODON_TO_AMINO;

void file2doubleHash(char *fname, int hsize) {

  FILE *fp;
  char *ptr;
  double *dptr;
  char **array;
  char buffer[2048];
  const char delim = ' ';
  
  fp = fopen(fname, "r");
  pq_create_hash(&CODON_TO_NSYN, hsize);
  array = calloc(2, sizeof (char*));
  
  while (fgets(buffer, sizeof(buffer), fp)) {
    if (pq_str2array(array, buffer, 2, &delim) == -1) {
      free(array);
      pq_free_hash(&CODON_TO_NSYN);
      fclose(fp);
      exit(1);
    }
    ptr = malloc(2048 * sizeof (char));
    strcpy(ptr, array[0]);
    dptr = malloc(sizeof (double));
    *dptr = atof(array[1]);
    pq_insert_hash(&CODON_TO_NSYN, ptr, dptr);
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
  pq_create_hash(&CODON_TO_AMINO, hsize);
  array = calloc(2, sizeof (char*));
  
  while (fgets(buffer, sizeof(buffer), fp)) {
    if (pq_str2array(array, buffer, 2, &delim) == -1) {
      free(array);
      pq_free_hash(&CODON_TO_AMINO);
      fclose(fp);
      exit(1);
    }
    ptr = malloc(2048 * sizeof (char));
    strcpy(ptr, array[0]);
    cptr = malloc(2048 * sizeof (char));
    strcpy(cptr, array[1]);
    pq_insert_hash(&CODON_TO_AMINO, ptr, cptr);
  }
  free(array);
  fclose(fp);
}
