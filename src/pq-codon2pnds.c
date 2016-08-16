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

#define CWIDTH 2048
#define LWIDTH 2048

/* 
   ATG M
   ATG 2
   
   
*/

void file2doublehash(char *fname, struct rwkHashTable *hash, int hsize) {

  FILE *fp;
  char *ptr;
  double *dptr;
  char **array;
  char buffer[LWIDTH];
  const char delim = ' ';
  
  fp = fopen(fname, "r");
  rwk_create_hash(hash, hsize);
  array = calloc(2, sizeof (char*));
  
  while (fgets(buffer, sizeof(buffer), fp)) {
    if (rwk_str2array(array, buffer, 2, &delim) == -1) {
      free(array);
      rwk_free_hash(hash);
      fclose(fp);
      exit(1);
    }
    ptr = malloc(CWIDTH * sizeof (char));
    strcpy(ptr, array[0]);
    dptr = malloc(sizeof (double));
    *dptr = atof(array[1]);
    rwk_insert_hash(hash, ptr, dptr);
  }
  free(array);
  fclose(fp);
}

/*
         ---> [ <char> <char> ... ]
 *ptr 
         ---> [ <double> ]

 */

void file2charhash(char *fname, struct rwkHashTable *hash, int hsize) {

  FILE *fp;
  char *ptr;
  char *cptr;
  char **array;
  char buffer[LWIDTH];
  const char delim = ' ';
  
  fp = fopen(fname, "r");
  rwk_create_hash(hash, hsize);
  array = calloc(2, sizeof (char*));
  
  while (fgets(buffer, sizeof(buffer), fp)) {
    if (rwk_str2array(array, buffer, 2, &delim) == -1) {
      free(array);
      rwk_free_hash(hash);
      fclose(fp);
      exit(1);
    }
    ptr = malloc(CWIDTH * sizeof (char));
    strcpy(ptr, array[0]);
    cptr = malloc(CWIDTH * sizeof (char));
    strcpy(cptr, array[1]);
    rwk_insert_hash(hash, ptr, cptr);
  }
  free(array);
  fclose(fp);
}

int main(int argc, char **argv) {

  /*
    
   ... | pq-pnds -f <int> -a <file> -s <file>
   
   */

  // parse arguments
  
  int i;
  int fcol;
  char *home;
  char *fname_aa;
  char *fname_syn;

  char default_fname_aa[] = ".config/pqgen/codon2aa";
  char default_fname_syn[] = ".config/pqgen/codon2syn";
  
  home = getenv("HOME");
  
  //char *default_fullpath_aa = malloc(strlen(home) + 2 + strlen(default_fname_aa));
  //char *default_fullpath_syn = malloc(strlen(home) + 2 + strlen(default_fname_syn));

  char default_fullpath_syn[1028];
  char default_fullpath_aa[1028];
  
  sprintf(default_fullpath_aa, "%s/%s", home, default_fname_aa);
  sprintf(default_fullpath_syn, "%s/%s", home, default_fname_syn);
  
  int default_fcol = 0;
  
  fcol = default_fcol;
  fname_aa = default_fullpath_aa;
  fname_syn = default_fullpath_syn;

  char usage[] = "usage: pq-codon2pnds [--help] [OPTIONS]\n";
  char options[] = "OPTIONS\n\n  -a <STR>\n     space delimited file where each line gives a '<codon> <amino_acid>'\n     key-value pair determining which amino acid a codon translates to. ['~/.config/pqgen/codon2aa']\n\n  -f <int>\n     column number (1-indexed) of the factor over which \n     the stats should be calculated. The default is to output \n     stats per chromosome, but the fourth name column could \n    be used instead to calculate over some group of features. [1]\n\n  -s <STR>\n     space delimited file where each line gives a '<codon> <syn_sites>'\n     key-value pair determining how many possible synonymous mutations there are\n     for each codon. ['~/.config/pqgen/codon2syn']\n\nInput:\nchr    codon.start    codon.end    name    score    +|-    codon.REF    codon.1 ... codon.n\n\nOutput:\nchr    codon.start    codon.end    name    score    +|-    ncodons    nvcodons    nsites.syn    nsites.nsyn    fix.syn    fix.nsyn    [poly.syn    poly.nsyn]\n\nExample input:\nchr1      676151  676154  transcript1       0.65   +       ATG     ATG\nchr1      676154  676157  transcript1       0.65   +       TCG     TCG\nchr1      676157  676160  transcript1       0.65   +       ACG     ACA\n\nExample output:\nchr1      676151  676154  transcript1       0.65   +       3	3	4.333333	4.666667	1	0\n";
  
  if (argc == 1 || strcmp(argv[1], "--help") == 0) {
    printf("%s\n", usage);
    printf("%s\n", options);
    exit(0);
  } else if (argc > 1) {
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
  
  int row1;
  int factor1;

  int nsam;
  int ncols;
  char **array;
  const char delim = '\t';

  char chr[CWIDTH];
  char factor[CWIDTH];
  char buffer[LWIDTH];

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

  ds = 0;
  dn = 0;
  pn = 0;
  ps = 0;
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
    if (rwk_str2array(array, buffer, ncols, &delim) == -1) {
      free(array);
      rwk_free_hash(&aaHash);
      rwk_free_hash(&synHash);
      exit(1);
    }
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

    unsigned int bad_codons;
    bad_codons = 0;
    for (i = 6; i < ncols; i++) {
      bad_codons += 1 - pq_alldna(array[i]);
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
