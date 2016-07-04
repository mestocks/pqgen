/*

  ... | pq-theta <nsam> [OPTIONS]

  <nsam> is a integer representing the total number of alleles 
  in the population. pq-stats expects stdin in bed format, with 
  the 5th and 6th columns representing the number of ref and 
  alt alleles respectively:
  
    chrom    pos-1    pos    name    nref    nalt

  OPTIONS

  -b <bool>
    takes values 0 or 1, indicating whether theta values should be 
    give per base pair or summed over the entire region. [1]
  
  -f <int>
    column number (1-indexed) of the factor over which 
    the stats should be calculated. The default is to output 
    stats per chromosome, but the fourth name column could 
    be used instead to calculate over some group of features. [1]

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pq_sfobj.h>
#include <rwk_parse.h>
#include <pq_sfstats.h>

int main(int argc, char **argv) {

  char help[] = "  ... | pqstats nsam [OPTIONS]";

  int i;
  int nsam;
  int fcol;
  int perbp;
  
  // <cmd> <nsam>
  fcol = 0;
  perbp = 1;
  if (argc == 1) {
    printf("%s\n", help);
  } else {
    nsam = atoi(argv[1]);
    for (i = 2; i < argc; i++) {
      if (strcmp(argv[i], "-b") == 0) {
	perbp = atoi(argv[i+1]);
      } else if (strcmp(argv[i], "-f") == 0) {
	fcol = atoi(argv[i+1]) - 1;
      }
    }
  }
  
  int ncols = 6;
  int lcols = 1024;
  int lwidth = 2048;
  char delim = '\t';
  char newline = '\n';
  char buffer[lwidth];

  char **array;
  array = calloc(ncols, sizeof (char*));
  
  char chr[lcols];
  char factor[lcols];
  long long int startpos, stoppos;
  long long int start_region, stop_region;

  int s, pi;
  int nref, nalt;
  double tw_val, pi_val;
  
  struct pqThetaW thetaW;
  struct pqThetaPi thetaPi;
  
  pqThetaWInit(&thetaW, nsam);
  pqThetaPiInit(&thetaPi, nsam);
  
  int coln;
  char *tmp;
  int startindex = 0;
  while (fgets(buffer, sizeof(buffer), stdin)) {
    rwkStrtoArray(array, buffer, &delim);
    
    startpos = atoll(array[1]);
    stoppos = atoll(array[2]);
    nref = atoi(array[4]);
    nalt = atoi(array[5]);

    if (startindex == 0) {
      strcpy(chr, array[0]);
      strcpy(factor, array[fcol]);
      start_region = startpos;
      startindex = 1;
    }

    if (strcmp(array[fcol], factor) != 0) {
      tw_val = thetaW.eval(&thetaW);
      pi_val = thetaPi.eval(&thetaPi);

      if (perbp == 1) {
	printf("%s\t%lld\t%lld\t%s\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",
	       chr, start_region, stop_region, factor,
	       thetaW.nsam, thetaW.nsites, thetaW.s,
	       tw_val / thetaW.nsites, pi_val / thetaPi.nsites,
	       pqTajimasD(thetaW.nsam, thetaW.s, tw_val, pi_val));
      } else {
	printf("%s\t%lld\t%lld\t%s\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",
	       chr, start_region, stop_region, factor,
	       thetaW.nsam, thetaW.nsites, thetaW.s,
	       tw_val, pi_val,
	       pqTajimasD(thetaW.nsam, thetaW.s, tw_val, pi_val));
      }
      
      thetaW.reset(&thetaW);
      thetaPi.reset(&thetaPi);
      strcpy(factor, array[fcol]);
      start_region = startpos;
      stop_region = stoppos;
    }
    
    if (nref > 0 && nalt > 0) {
      s = 1;
      pi = nref;
    } else {
      s = 0;
      pi = 0;
    }

    if (nref + nalt == nsam) {
      thetaW.add(&thetaW, s);
      thetaPi.add(&thetaPi, pi);
      stop_region = stoppos;
    }
    strcpy(chr, array[0]);
  }

  tw_val = thetaW.eval(&thetaW);
  pi_val = thetaPi.eval(&thetaPi);
  if (perbp == 1) {
    printf("%s\t%lld\t%lld\t%s\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",
	   chr, start_region, stop_region, factor,
	   thetaW.nsam, thetaW.nsites, thetaW.s,
	   tw_val / thetaW.nsites, pi_val / thetaPi.nsites,
	   pqTajimasD(thetaW.nsam, thetaW.s, tw_val, pi_val));
  } else {
    printf("%s\t%lld\t%lld\t%s\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",
	   chr, start_region, stop_region, factor,
	   thetaW.nsam, thetaW.nsites, thetaW.s,
	   tw_val, pi_val,
	   pqTajimasD(thetaW.nsam, thetaW.s, tw_val, pi_val));
  }
    
  free(array);
  
  return 0; }
