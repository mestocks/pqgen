#include <stdio.h>
#include <stdlib.h>

#include "pq_genetics.h"

void pqComplement(char *codon) {
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

void pqReverse(char *codon) {
  char one  = codon[0];
  char three = codon[2];
  codon[0] = three;
  codon[2] = one;
}
