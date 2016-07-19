/*
... | dna2codon

Bed format input:
chr    start    end    name    score    +|-    nuc.1    nuc.2 ... nuc.n

Bed format output:
chr    codon.start    codon.end    name    score    +|-    codon.1    codon.2 ... codon.n    

Example input:
chr1      676151  676152  transcript1       0.65   +       A       A
chr1      676152  676153  transcript1       0.65   +       T       T
chr1      676153  676154  transcript1       0.65   +       G       G
chr1      676154  676155  transcript1       0.65   +       T       T
chr1      676155  676156  transcript1       0.65   +       C       C
chr1      676156  676157  transcript1       0.65   +       G       G
chr1      676157  676158  transcript1       0.65   +       A       A
chr1      676158  676159  transcript1       0.65   +       C       C
chr1      676159  676160  transcript1       0.65   +       G       A

Example output:
chr1      676151  676154  transcript1       0.65   +       ATG     ATG
chr1      676154  676157  transcript1       0.65   +       TCG     TCG
chr1      676157  676160  transcript1       0.65   +       ACG     ACA
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_parse.h>
#include <pq_genetics.h>

int main(int argc, char **argv) {

  char usage[] = "usage: pq-dna2codon [--help]\n";
  char options[] = "Input:\nchr    start    end    name    score    +|-    nuc.1    nuc.2 ... nuc.n\n\nOutput:\nchr    codon.start    codon.end    name    score    +|-    codon.1    codon.2 ... codon.n\n\nExample input:\nchr1      676151  676152  transcript1       0.65   +       A       A\nchr1      676152  676153  transcript1       0.65   +       T       T\nchr1      676153  676154  transcript1       0.65   +       G       G\nchr1      676154  676155  transcript1       0.65   +       T       T\nchr1      676155  676156  transcript1       0.65   +       C       C\nchr1      676156  676157  transcript1       0.65   +       G       G\nchr1      676157  676158  transcript1       0.65   +       A       A\nchr1      676158  676159  transcript1       0.65   +       C       C\nchr1      676159  676160  transcript1       0.65   +       G       A\n\nExample output:\nchr1      676151  676154  transcript1       0.65   +       ATG     ATG\nchr1      676154  676157  transcript1       0.65   +       TCG     TCG\nchr1      676157  676160  transcript1       0.65   +       ACG     ACA\n";
  
  if (strcmp(argv[1], "--help") == 0) {
    printf("%s\n", usage);
    printf("%s\n", options);
    exit(0);
  }
  
  int ismp;
  int irow;
  int nnuc;
  int ncols;
  int start;
  int lwidth = 2048;
  const char delim = '\t';
  char buffer[lwidth];

  char **array;
  char **codons;

  irow = 0;
  start = 0;
  
  while (fgets(buffer, sizeof(buffer), stdin)) {

    if (irow == 0) {
      ncols = rwk_countcols(buffer, &delim);
      //ncols = rwkCountCols(buffer, old_delim);
      array = calloc(ncols, sizeof (char*));
      nnuc = ncols - 6;
      codons = calloc(nnuc, sizeof (char*));
      for (ismp = 0; ismp < nnuc; ismp++) {
	codons[ismp] = calloc(4, sizeof (char));
	codons[ismp][3] = '\0';
      }
    }
    //rwkStrtoArray(array, buffer, &old_delim);
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
  
  return 0; }
