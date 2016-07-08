/*
... | codon2nonsyn

Bed format output:
chr    codon.start    codon.end    name    score    +|-    codon.1    codon.2 ... codon.n    

Example input:
chr1      676151  676154  transcript1       0.65   +       ATG     ATG
chr1      676154  676157  transcript1       0.65   +       TCG     TCG
chr1      676157  676160  transcript1       0.65   +       ACG     ACA
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_parse.h>
#include <rwk_htable.h>
#include <pq_genetics.h>

void file2hash(char *fname, struct rwkSmallHashTable *hash) {

  FILE *fp;
  int hsize;
  int lwidth;
  char **array;
  char buffer[lwidth];
  const char delim = ' ';
  
  lwidth = 128;
  hsize = 6000000;
  fp = fopen(fname, "r");
  rwk_create_hash(hash, hsize);
  array = calloc(2, sizeof (char*));
  
  while (fgets(buffer, sizeof(buffer), fp)) {
    rwk_strsplit(array, buffer, &delim);
    rwk_insert_hash(hash, array[0], atof(array[1]));
  }
}

int main(int argc, char **argv) {

  struct rwkSmallHashTable myHash;
  file2hash(argv[1], &myHash);
  
  /*
  
  int ismp;
  int irow;
  int nnuc;
  int ncols;
  int start;
  int lwidth = 2048;
  const char delim = '\t';
  char old_delim = '\t';
  char buffer[lwidth];

  char **array;
  char **codons;

  irow = 0;
  start = 0;
  
  while (fgets(buffer, sizeof(buffer), stdin)) {

    if (irow == 0) {
      ncols = rwk_countcols(buffer, &delim);
      array = calloc(ncols, sizeof (char*));
      nnuc = ncols - 6;
      codons = calloc(nnuc, sizeof (char*));
      for (ismp = 0; ismp < nnuc; ismp++) {
	codons[ismp] = calloc(4, sizeof (char));
	codons[ismp][3] = '\0';
      }
    }
    rwk_strsplit(array, buffer, &delim);

    for (ismp = 0; ismp < nnuc; ismp++) {
      codons[ismp][start] = *array[ismp+6];
    }
    
    if (start == 2) {
      printf("%s\t%d\t%s\t%s\t%s\t%s",
	     array[0], atoi(array[1]) - 2, array[2],
	     array[3], array[4], array[5]);
      
      for (ismp = 0; ismp < nnuc; ismp++) {
	if ('-' == *array[5]) {
	  pq_reverse(codons[ismp]);
	  pq_complement(codons[ismp]);
	}
	printf("\t%s", codons[ismp]);
      }
      printf("\n");
      start = 0;
    } else {
      start++;
    }
        
    irow++;
  }

  for (ismp = 0; ismp < nnuc; ismp++) {
    free(codons[ismp]);
  }
  free(codons);
  free(array);

*/

  rwk_free_hash(&myHash);
  
  return 0; }
