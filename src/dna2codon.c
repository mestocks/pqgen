/*

... | dna2codon

chr    pos-1    pos    nuc

 0       1       2      3        4       5         6
chr    pos-1    pos    name    score    +|-    A|T|G|C|N|-

chr    pos-1    pos    codon    
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_parse.h>

int main(int argc, char **argv) {

  int irow;
  int ncols;
  int start;
  int lcols = 1024;
  int lwidth = 2048;
  char delim = '\t';
  char buffer[lwidth];

  char **array;
  char codon[4];

  irow = 0;
  start = 0;
  codon[3] = '\0';
  
  while (fgets(buffer, sizeof(buffer), stdin)) {

    if (irow == 0) {
      ncols = rwkCountCols(buffer, delim);
      array = calloc(ncols, sizeof (char*));
    }
    rwkStrtoArray(array, buffer, &delim);

    codon[start] = *array[6];
    if (start == 2) {
      if ('+' == *array[5]) {
      printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\n",
	     array[0], atoi(array[1]) - 2, array[2],
	     array[3], array[4], array[5], codon);
      }
      else if ('-' == *array[5]) {
	printf("%s\t%d\t%s\t%s\t%s\t%s\t%c%c%c\n",
	       array[0], atoi(array[1]) - 2, array[2],
	       array[3], array[4], array[5], codon[2], codon[1], codon[0]);
      }
      start = 0;
    } else {
      start++;
    }
    irow++;
  }
  
  free(array);
  
  return 0; }
