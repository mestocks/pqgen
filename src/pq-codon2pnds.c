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

void file2doublehash(char *fname, struct rwkHashTable *hash, int hsize) {

  FILE *fp;
  int lwidth;
  int cwidth;
  char **array;
  
  lwidth = 128;
  cwidth = 128;
  char buffer[lwidth];
  const char delim = ' ';
  
  fp = fopen(fname, "r");
  rwk_create_hash(hash, hsize);
  array = calloc(2, sizeof (char*));
  
  while (fgets(buffer, sizeof(buffer), fp)) {
    rwk_strsplit(array, buffer, &delim);
    char *ptr = malloc(cwidth * sizeof (char));
    strcpy(ptr, array[0]);
    double *dptr = malloc(sizeof (double));
    *dptr = atof(array[1]);
    rwk_insert_hash(hash, ptr, dptr);
  }
}

void file2charhash(char *fname, struct rwkHashTable *hash, int hsize) {

  FILE *fp;
  int lwidth;
  int cwidth;
  char **array;

  lwidth = 128;
  cwidth = 128;
  char buffer[lwidth];
  const char delim = ' ';
  
  fp = fopen(fname, "r");
  rwk_create_hash(hash, hsize);
  array = calloc(2, sizeof (char*));
  while (fgets(buffer, sizeof(buffer), fp)) {
    rwk_strsplit(array, buffer, &delim);
    char *ptr = malloc(cwidth * sizeof (char));
    strcpy(ptr, array[0]);
    char *cptr = malloc(cwidth * sizeof (char));
    strcpy(cptr, array[1]);
    rwk_insert_hash(hash, ptr, cptr);
  }
}

int main(int argc, char **argv) {

  /*
    
   ... | pq-pnds -f <int> -a <file> -s <file>
   
   */

  // parse arguments
  
  int i;
  int fcol;
  char *fname_aa;
  char *fname_syn;

  int default_fcol = 0;
  char default_fname_aa[] = "share/codon2aa";
  char default_fname_syn[] = "share/codon2syn";

  fcol = default_fcol;
  fname_aa = default_fname_aa;
  fname_syn = default_fname_syn;
  
  if (argc != 1) {
    for (i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-f") == 0) {
	fcol = atoi(argv[i + 1]) - 1;
      } else if (strcmp(argv[i], "-a") == 0) {
	fname_aa = argv[i + 1];
      } else if (strcmp(argv[i], "-s") == 0) {
	fname_syn = argv[i + 1];
      }
    }
  }

  // fill hash tables
  
  int hsize;
  struct rwkHashTable aaHash;
  struct rwkHashTable synHash;

  hsize = 128;
  file2charhash(fname_aa, &aaHash, hsize);
  file2doublehash(fname_syn, &synHash, hsize);

  //

  int row1;
  int factor1;

  int nsam;
  int ncols;
  int cwidth;
  int lwidth;
  char **array;
  const char delim = '\t';

  cwidth = 2048;
  lwidth = 2048;

  char chr[cwidth];
  char factor[cwidth];
  char buffer[lwidth];

  void *aptr;
  void *ref_aa;
  void *ref_sptr;
  double syn_muts;
  double nsyn_muts;

  int dn;
  int ds;
  int pn;
  int ps;
  int ncodons;
  int nvcodons;
  int diff_syn;
  int diff_nsyn;
  double synonymous;
  double nonsynonymous;

  long long int startpos;
  long long int stoppos;
  long long int start_region;
  long long int stop_region;
  
  row1 = 0;
  factor1 = 0;
  ncodons = 0;
  nvcodons = 0;
  diff_syn = 0;
  diff_nsyn = 0;
  synonymous = 0.0;
  nonsynonymous = 0.0;
  
  while (fgets(buffer, sizeof(buffer), stdin)) {
    
    if (row1 == 0) {
      ncols = rwk_countcols(buffer, &delim);
      nsam = ncols - 7;
      array = calloc(ncols, sizeof (char*));
      row1 = 1;
    }
    rwk_strsplit(array, buffer, &delim);
    startpos = atoll(array[1]);
    stoppos = atoll(array[2]);
    
    if (factor1 == 0) {
      strcpy(chr, array[0]);
      strcpy(factor, array[fcol]);
      start_region = startpos;
      factor1 = 1;
    }

    if (strcmp(array[fcol], factor) != 0) {
      if (nsam == 1) {
	printf("%s\t%lld\t%lld\t%s\t%d\t%d\t%f\t%f\t%d\t%d\n", chr, start_region, stop_region, factor, ncodons, nvcodons, synonymous, nonsynonymous, ds, dn);
      } else if (nsam > 1) {
	printf("%s\t%lld\t%lld\t%s\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%d\n", chr, start_region, stop_region, factor, ncodons, nvcodons, synonymous, nonsynonymous, ds, dn, ps, pn);
      }
      ds = 0;
      dn = 0;
      pn = 0;
      ps = 0;
      ncodons = 0;
      nvcodons = 0;
      synonymous = 0.0;
      nonsynonymous = 0.0;
      strcpy(factor, array[fcol]);
      start_region = startpos;
      stop_region = stoppos;
    }

    unsigned int ok_codon;
    unsigned int bad_codons = 0;
    for (i = 6; i < ncols; i++) {
      ok_codon = pq_alldna(array[i]);
      bad_codons += 1 - ok_codon;
    }
    
    if (bad_codons == 0) {
      pq_dna_upper(array[6]);
      ref_aa = rwk_lookup_hash(&aaHash, array[6]);
      ref_sptr = rwk_lookup_hash(&synHash, array[6]);
      syn_muts = *(double *)ref_sptr;
      nsyn_muts = 9.0 - syn_muts;
      
      synonymous += syn_muts / 3.0;
      nonsynonymous += nsyn_muts / 3.0;
      diff_syn = 0;
      diff_nsyn = 0;
      for (i = 7; i < ncols; i++) {
	pq_dna_upper(array[i]);
	aptr = rwk_lookup_hash(&aaHash, array[i]);
	if (strcmp(array[6], array[i]) == 0) {
	} else {
	  if (strcmp(aptr, ref_aa) != 0) {
	    diff_nsyn++;
	  } else {
	    diff_syn++;
	  }
	}
      }
      if (nsam == diff_syn) {
	ds++;
      } else if (diff_syn > 0) {
	ps++;
      }
      if (nsam == diff_nsyn) {
	dn++;
      } else if (diff_nsyn > 0) {
	pn++;
      }
      nvcodons++;
    }
    strcpy(chr, array[0]);
    stop_region = stoppos;
    ncodons++;
  }
  if (nsam == 1) {
    printf("%s\t%lld\t%lld\t%s\t%d\t%d\t%f\t%f\t%d\t%d\n", chr, start_region, stop_region, factor, ncodons, nvcodons, synonymous, nonsynonymous, ds, dn);
  } else if (nsam > 1) {
    printf("%s\t%lld\t%lld\t%s\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%d\n", chr, start_region, stop_region, factor, ncodons, nvcodons, synonymous, nonsynonymous, ds, dn, ps, pn);
  }
  
  free(array);
  rwk_free_hash(&aaHash);
  rwk_free_hash(&synHash);
  
  return 0; }
