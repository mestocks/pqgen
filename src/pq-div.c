#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_args.h>
#include <rwk_parse.h>
#include <rwk_htable.h>
#include <pq_generics.h>

#define LCOL 1024
#define LWIDTH 2048

int main(int argc, char **argv)
{
  char usage[] = "usage: pq-div [--help] [OPTIONS]\n";
  char options[] = "OPTIONS\n\n  -f <int>\n    column number (1-indexed) of the factor over which\n    the stats should be calculated. The default is to output \n    stats per chromosome, but the fourth name column could \n    be used instead to calculate over some group of features. [1]\n";

  if (argc == 2 && strcmp(argv[1], "--help") == 0) {
    printf("%s\n", usage);
    printf("%s\n", options);
    exit(0);
  }
  
  int nargs;
  char **def_array;
  unsigned int REF;
  unsigned int ALT;
  struct pq_parameters params;
  char defaults[] = "-f 1 -c 1 -p 3 -k 5,6";
  nargs = rwk_countcols(defaults, " ");
  def_array = calloc(nargs, sizeof(char *));
  rwk_str2array(def_array, defaults, nargs, " ");

  pq_init_parameters(&params);
  params.update(&params, nargs, def_array);
  params.update(&params, argc-1, argv+1);

  REF = params.KCOLS[0];
  ALT = params.KCOLS[1];
  free(def_array);

  
  unsigned int i;
  int *pint;
  char *pchar;
  char dna[] = "ATGCatgc";
  struct rwkHashTable DnaHash;
  
  rwk_create_hash(&DnaHash, 20);
  for (i = 0; i < strlen(dna); i++) {
    pint = malloc(sizeof (int));
    *pint = 1;
    pchar = malloc(2 * sizeof (char));
    pchar[0] = dna[i];
    pchar[1] = '\0';
    rwk_insert_hash(&DnaHash, pchar, pint);
  }
  
  int ncols;
  char **array;
  char buffer[LWIDTH];
  const char delim = '\t';
  
  char chr[LCOL];
  char factor[LCOL];
  unsigned long long int startpos;
  unsigned long long int stoppos;
  unsigned long long int start_region;
  unsigned long long int stop_region;
  
  char value[LCOL];
  char missing[LCOL];
  struct Data data;
  
  strcpy(missing, "NA");

  double dist = 0.0;
  unsigned long long int nvsites = 0;
  unsigned long long int diff = 0;
  unsigned int startindex = 0;

  int row1 = 0;
  while (fgets(buffer, sizeof(buffer), stdin)) {

    if (row1 == 0) {
      ncols = rwk_countcols(buffer, &delim);
      array = calloc(ncols, sizeof (char*));
      row1 = 1;
    }
    
    if (rwk_str2array(array, buffer, ncols, &delim) == -1) {
      free(array);
      rwk_free_hash(&DnaHash);
      rwk_free_hash(&DnaHash);
      pq_free_parameters(&params);
      exit(1);
    }
    
    // seg fault if row contains < ncols columns
    stoppos = atoll(array[params.POS]);
    startpos = stoppos - 1;

    if (startindex == 0) {
      strcpy(chr, array[params.CHROM]);
      strcpy(factor, array[params.FCOL]);
      start_region = startpos;
      startindex = 1;
    }

    if (strcmp(array[params.FCOL], factor) != 0) {
      if (nvsites == 0) {
	sprintf(value, "%s", missing);
      } else {
	if (diff == 0) {
	  dist = 0.0;
	} else {
	  dist = (double)diff / (double)nvsites;
	}
	sprintf(value, "%f", dist);
      }

      data.chr = chr;
      data.start = &start_region;
      data.end = &stop_region;
      data.name = factor;
      data.nvsites = &nvsites;
      data.value = value;
      pq_print_bed(&data);

      diff = 0;
      dist = 0;
      nvsites = 0;
      strcpy(factor, array[params.FCOL]);
      start_region = startpos;
      stop_region = stoppos;
    }
        
    if (rwk_lookup_hash(&DnaHash, array[REF]) != NULL && rwk_lookup_hash(&DnaHash, array[ALT]) != NULL) {
      nvsites++;
      if (strcmp(array[REF], array[ALT]) != 0) {
	diff++;
      }
      stop_region = stoppos;
    }
    strcpy(chr, array[params.CHROM]);
  }
  
  if (nvsites == 0) {
    sprintf(value, "%s", missing);
  } else {
    if (diff == 0) {
      dist = 0.0;
    } else {
      dist = (double)diff / (double)nvsites;
    }
    sprintf(value, "%f", dist);
  }
  
  data.chr = chr;
  data.start = &start_region;
  data.end = &stop_region;
  data.name = factor;
  data.nvsites = &nvsites;
  data.value = value;
  pq_print_bed(&data);
  
  free(array);
  rwk_free_hash(&DnaHash);
  pq_free_parameters(&params);
  
  return 0;
}
