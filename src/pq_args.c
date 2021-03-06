#include <stdlib.h>
#include <stdio.h>

#include "pq_parse.h"
#include "pq_htable.h"

int *KCOLS;
int NKARGS;
unsigned int CHROM;
unsigned int POS;
unsigned int FCOL;
struct HashTable ARGHASH;

void pq_update_args(unsigned int argc, char **argv)
{
  int k;
  char *chr_ptr;
  char *end_ptr;
  char *kcol_ptr;
  char *fcol_ptr;
  
  pq_multi_insert_hash(&ARGHASH, argc, argv);
  chr_ptr = (char *)pq_lookup_hash(&ARGHASH, "-c");
  end_ptr = (char *)pq_lookup_hash(&ARGHASH, "-p");
  kcol_ptr = (char *)pq_lookup_hash(&ARGHASH, "-k");
  fcol_ptr = (char *)pq_lookup_hash(&ARGHASH, "-f");
  
  CHROM = atoi(chr_ptr) - 1;
  POS = atoi(end_ptr) - 1;
  FCOL = atoi(fcol_ptr) - 1;
  
  if (KCOLS != NULL) {
    free(KCOLS);
  }

  NKARGS = pq_nkrange(kcol_ptr);
  KCOLS = malloc(NKARGS * sizeof(int));
  pq_k2array(KCOLS, kcol_ptr, NKARGS);
  for (k = 0; k < NKARGS; k++) {
    KCOLS[k]--;
  }
}

void pq_init_args()
{
  pq_create_hash(&ARGHASH, 128);
  KCOLS = NULL;
}

void pq_free_args()
{
  free(KCOLS);
  pq_free_hash(&ARGHASH);
}
