#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pq_args.h>
#include <pq_genetics.h>

#include <rwk_args.h>
#include <rwk_parse.h>
#include <rwk_htable.h>

#include <pq_parse.h>
#include <pq_generics.h>

#include <pq_div.h>
#include <pq_het.h>
#include <pq_pnds.h>
#include <pq_sfs.h>
#include <pq_theta.h>

#include <pq_doc.h>
#include <pq_version.h>
#include <rwk_version.h>

#define LCOL 1024
#define LWIDTH 5120

void display_help()
{
  printf("usage: %s\n\n", PQ_USAGE);
  printf("   %s        %s\n", PQ_DIV_NAME, PQ_DIV_DESC);
  printf("   %s        %s\n", PQ_HET_NAME, PQ_HET_DESC);
  printf("   %s       %s\n", PQ_PNDS_NAME, PQ_PNDS_DESC);
  printf("   %s        %s\n", PQ_SFS_NAME, PQ_SFS_DESC);
  printf("   %s      %s\n", PQ_THETA_NAME, PQ_THETA_DESC);
}

void display_version()
{
  printf("pqgen version %s (librawk version %s)\n", PQ_VERSION, RWK_VERSION);
}

int main(int argc, char **argv)
{
  if (argc == 1) {
    display_help();
    exit(0);
  }
  if (argc > 1 && strcmp(argv[1], "--help") == 0) {
    display_help();
    exit(0);
  }
  if (argc > 1 && strcmp(argv[1], "--version") == 0) {
    display_version();
    exit(0);
  }
  
  int frm_multi;
  struct SWrap wrap;
  void (*swrap_init)(struct SWrap *, int);

  /*
    swrap_init function;
    frm_multi (ploidy);
    
   */
   
  if (strcmp(argv[1], "theta") == 0) {
    frm_multi = 2;
    swrap_init = pq_swinit_theta;
  } else if (strcmp(argv[1], "sfs") == 0) {
    frm_multi = 2;
    swrap_init = pq_swinit_sfs;
  } else if (strcmp(argv[1], "het") == 0) {
    frm_multi = 1;
    swrap_init = pq_swinit_het;
  } else if (strcmp(argv[1], "div") == 0) {
    frm_multi = 1;
    swrap_init = pq_swinit_div;
  } else if (strcmp(argv[1], "pnds") == 0) {
    char *home;
    char *fname_aa;
    char *fname_syn;
    char default_fname_aa[] = ".config/pqgen/codon2aa";
    char default_fname_syn[] = ".config/pqgen/codon2syn";
    home = getenv("HOME");
    char default_fullpath_aa[1028];
    char default_fullpath_syn[1028];
    sprintf(default_fullpath_aa, "%s/%s", home, default_fname_aa);
    sprintf(default_fullpath_syn, "%s/%s", home, default_fname_syn);
    fname_aa = default_fullpath_aa;
    fname_syn = default_fullpath_syn;
    file2charHash(fname_aa, 128);
    file2doubleHash(fname_syn, 128);
    frm_multi = 1;
    swrap_init = pq_swinit_pnds;
  } else {
    printf("pqgen: '%s' is not a command. See 'pqgen --help'.\n", argv[1]);
    exit(0);
  }
  
  int nargs;
  char **def_array;
  char defaults[128];
  
  int ncols;
  int nalleles;
  char buffer[LWIDTH];
  const char delim = '\t';

  char chr[LCOL];
  char factor[LCOL];
  unsigned long long int startpos;
  unsigned long long int stoppos;
  unsigned long long int start_region;
  unsigned long long int stop_region;

  int i;
  int row_index;
  int startindex;

  row_index = 0;
  startindex = 0;
  
  struct GenericRow row;
  
  while (fgets(buffer, sizeof(buffer), stdin)) {
    if (row_index == 0) {
      ncols = rwk_countcols(buffer, &delim);

      if (strcmp(argv[1], "theta") == 0) {
	sprintf(defaults, "-f 1 -c 1 -p 3 -k 5-%d -b 1", ncols);
      } else if (strcmp(argv[1], "het") == 0) {
	sprintf(defaults, "-f 1 -c 1 -p 3 -k 5-%d", ncols);
      } else if (strcmp(argv[1], "sfs") == 0) {
	sprintf(defaults, "-f 1 -c 1 -p 3 -k 5-%d", ncols);
      } else if (strcmp(argv[1], "div") == 0) {
	sprintf(defaults, "-f 1 -c 1 -p 3 -k 5,6 -b 1");
      } else if (strcmp(argv[1], "pnds") == 0) {
	sprintf(defaults, "-f 1 -c 1 -p 3 -k 7-%d", ncols);
      }
      
      nargs = rwk_countcols(defaults, " ");
      def_array = calloc(nargs, sizeof(char *));
      rwk_str2array(def_array, defaults, nargs, " ");
      
      pq_init_args();
      pq_update_args(nargs, def_array);
      pq_update_args(argc-2, argv+2);
      
      free(def_array);

      init_row(&row, ncols, CHROM, POS, FCOL);

      nalleles = frm_multi * NKARGS;
      swrap_init(&wrap, nalleles);
      
      row_index = 1;
    }
    
    row.update(&row, buffer, &delim);

    stoppos = row.pos(&row);
    startpos = stoppos - 1;

    if (startindex == 0) {
      strcpy(chr, row.chrom(&row));
      strcpy(factor, row.factor(&row));
      start_region = startpos;
      startindex = 1;
    }

    if (strcmp(row.factor(&row), factor) != 0) {
      wrap.write(&wrap);
      printf("%s\t%llu\t%llu\t%s", chr, start_region, stop_region, factor);
      for (i = 0; i < wrap.nouts; i++) {
	printf("\t%s", (char *)wrap.outs[i]);
      }  
      printf("\n");
      wrap.clear(&wrap);      
      
      strcpy(factor, row.factor(&row));
      start_region = startpos;
      stop_region = stoppos;   
    }

    wrap.update(&wrap, row.array);
    
    stop_region = stoppos;
    strcpy(chr, row.chrom(&row));
  }

  wrap.write(&wrap);
  printf("%s\t%llu\t%llu\t%s", chr, start_region, stop_region, factor);
  for (i = 0; i < wrap.nouts; i++) {
    printf("\t%s", (char *)wrap.outs[i]);
  }
  printf("\n");
  wrap.clear(&wrap);

  rwk_free_hash(&CODON_TO_NSYN);
  rwk_free_hash(&CODON_TO_AMINO);
  pq_free_args();
  pq_swfree(&wrap);
  free_row(&row);
  
  return 0;
}
