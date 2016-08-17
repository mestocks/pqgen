#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_parse.h>
#include <rwk_htable.h>

#define LCOL 1024
#define LWIDTH 2048

struct Data {
  char *chr;
  unsigned long long int *start;
  unsigned long long int *end;
  char *name;
  unsigned long long int *nvsites;
  char *value;
};

void print_bed(struct Data *data)
{
  printf("%s\t%llu\t%llu\t%s\t%llu\t%s\n", data->chr, *data->start, *data->end, data->name, *data->nvsites, data->value);
}

int main(int argc, char **argv)
{
  unsigned int CHR = 0;
  unsigned int START = 1;
  unsigned int END = 2;
  unsigned int REF = 4;
  unsigned int ALT = 5;
  unsigned int FCOL = 0;

  char usage[] = "usage: pq-div [--help] [OPTIONS]\n";
  char options[] = "OPTIONS\n\n  -f <int>\n    column number (1-indexed) of the factor over which\n    the stats should be calculated. The default is to output \n    stats per chromosome, but the fourth name column could \n    be used instead to calculate over some group of features. [1]\n";

  if (argc == 2 && strcmp(argv[1], "--help") == 0) {
    printf("%s\n", usage);
    printf("%s\n", options);
    exit(0);
  } else if (argc > 1 && strcmp(argv[1], "-f") == 0) {
    FCOL = atoi(argv[2]) - 1;
  }

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
  
  int ncols = 6;
  char buffer[LWIDTH];
  const char delim = '\t';
  
  char **array;
  array = calloc(ncols, sizeof (char*));
  
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
  
  while (fgets(buffer, sizeof(buffer), stdin)) {
    
    if (rwk_str2array(array, buffer, ncols, &delim) == -1) {
      free(array);
      rwk_free_hash(&DnaHash);
      exit(1);
    }
    
    // seg fault if row contains < ncols columns
    startpos = atoll(array[START]);
    stoppos = atoll(array[END]);

    if (startindex == 0) {
      strcpy(chr, array[CHR]);
      strcpy(factor, array[FCOL]);
      start_region = startpos;
      startindex = 1;
    }

    if (strcmp(array[FCOL], factor) != 0) {
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
      print_bed(&data);

      diff = 0;
      dist = 0;
      nvsites = 0;
      strcpy(factor, array[FCOL]);
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
    strcpy(chr, array[CHR]);
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
  print_bed(&data);
  
  free(array);
  rwk_free_hash(&DnaHash);
  
  return 0;
}
