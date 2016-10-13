#include <stdlib.h>

#include <rwk_args.h>
#include <rwk_htable.h>

int *KCOLS;
int NKARGS;
unsigned int CHROM;
unsigned int POS;
unsigned int FCOL;
struct rwkHashTable ARGHASH;

void pq_update_args(unsigned int argc, char **argv)
{  
  char *chr_ptr;
  char *end_ptr;
  char *kcol_ptr;
  char *fcol_ptr;
  
  rwk_multi_insert_hash(&ARGHASH, argc, argv);
  chr_ptr = (char *)rwk_lookup_hash(&ARGHASH, "-c");
  end_ptr = (char *)rwk_lookup_hash(&ARGHASH, "-p");
  kcol_ptr = (char *)rwk_lookup_hash(&ARGHASH, "-k");
  fcol_ptr = (char *)rwk_lookup_hash(&ARGHASH, "-f");
  
  CHROM = atoi(chr_ptr) - 1;
  POS = atoi(end_ptr) - 1;
  FCOL = atoi(fcol_ptr) - 1;
  
  if (KCOLS != NULL) {
    free(KCOLS);
  }

  NKARGS = rwk_countkargs(kcol_ptr);
  KCOLS = malloc(NKARGS * sizeof(int));
  rwk_kargsplit(KCOLS, kcol_ptr);
}

void pq_init_args()
{
  rwk_create_hash(&ARGHASH, 128);
  KCOLS = NULL;
}

void pq_free_args()
{
  free(KCOLS);
  rwk_free_hash(&ARGHASH);
}
