#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_htable.h>

#include <pq_args.h>
#include <pq_generics.h>

/*
  #/#  #/#  #/#  #/#

  .genotypes
  [ value
  [ count
  [ next
  
*/

struct genotype_counts {
  char *value;
  int count;
  struct genotype_counts *next;
};

void pq_swupdate_het(struct SWrap *wrap, char **array)
{
  int i;
  int n;
  int nhet;
  char *curr;
  char *a1;
  char *a2;
  double het;

  n = 0;
  het = 0.0;
  nhet = 0;
  for (i = 0; i < NKARGS; i++) {
    curr = array[KCOLS[i]];
    a1 = &curr[0];
    if (curr[0] != '.') {
      if (curr[1] == '/' && curr[3] == '\0') {
	a2 = &curr[2];
	if (*a1 != *a2) {
	  nhet++;
	}
	n++;
      } else {
	// ignored
	// indel. e.g. AAT/AAAAT, AAT/A or A/AAAAGC
      }
    }
  }
  het = (double)nhet / (double)n;
  if (n > 0) {
    *(long long int *)wrap->values[0] += 1;
    *(double *)wrap->values[1] += het;
  }
}

void pq_swwrite_het(struct SWrap *wrap)
{
  double het;
  sprintf(wrap->outs[0], "%d", wrap->nsam);
  sprintf(wrap->outs[1], "%lli", *(long long int *)wrap->values[0]);

  if (strcmp((char *)rwk_lookup_hash(&ARGHASH, "-b"), "0") == 0) {
    het = *(double *)wrap->values[1];
  } else {
    het = *(double *)wrap->values[1] / *(long long int *)wrap->values[0];
  }
  sprintf(wrap->outs[2], "%f", het);
}

void pq_swclear_het(struct SWrap *wrap)
{
  *(long long int *)wrap->values[0] = 0;
  *(double *)wrap->values[1] = 0.0;
}

// nvalues = [nvsites, het]
// nouts   = [nsam, nvsites, uhet]

void PQ_HET_INIT(struct SWrap *wrap, int nsam)
{
  int i;
  wrap->nsam = nsam;
  wrap->nouts = 3;
  wrap->nvalues = 2;
  wrap->outs = calloc(wrap->nouts, sizeof(char *));
  wrap->values = calloc(wrap->nvalues, sizeof(void *));
  for (i = 0; i < wrap->nouts; i++) {
    wrap->outs[i] = calloc(128, sizeof(char));
  }
  
  wrap->values[0] = malloc(sizeof(long long int));
  *(long long int *)wrap->values[0] = 0;
    wrap->values[1] = malloc(sizeof(double));
  *(double *)wrap->values[1] = 0.0;
  
  wrap->update = pq_swupdate_het;
  wrap->write = pq_swwrite_het;
  wrap->clear = pq_swclear_het;
}

