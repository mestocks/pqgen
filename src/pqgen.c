#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_args.h>
#include <rwk_parse.h>
#include <rwk_htable.h>

#include <pq_args.h>
#include <pq_parse.h>
#include <pq_limits.h>
#include <pq_generics.h>
#include <pq_genetics.h>

#include <pq_div.h>
#include <pq_het.h>
#include <pq_sfs.h>
#include <pq_pnds.h>
#include <pq_theta.h>

#include <pq_help.h>
#include <pq_version.h>

struct pq_command {
  int frmt;
  const char *name;
  const char *desc;
  void (*init)(struct SWrap *, int);
};

struct pq_command CMD[] =
  {
    {.name = PQ_DIV_NAME,
     .desc = PQ_DIV_DESC,
     .frmt = PQ_DIV_FRMT,
     .init = PQ_DIV_INIT},
    
    {.name = PQ_HET_NAME,
     .desc = PQ_HET_DESC,
     .frmt = PQ_HET_FRMT,
     .init = PQ_HET_INIT},
    
    {.name = PQ_PNDS_NAME,
     .desc = PQ_PNDS_DESC,
     .frmt = PQ_PNDS_FRMT,
     .init = PQ_PNDS_INIT},
    
    {.name = PQ_SFS_NAME,
     .desc = PQ_SFS_DESC,
     .frmt = PQ_SFS_FRMT,
     .init = PQ_SFS_INIT},
    
    {.name = PQ_THETA_NAME,
     .desc = PQ_THETA_DESC,
     .frmt = PQ_THETA_FRMT,
     .init = PQ_THETA_INIT}
  };

void display_help()
{
  fprintf(stderr, "usage: %s\n\n", PQ_USAGE);

  int i;
  int ncmds;
  int tsize;
  tsize = 11;
  ncmds = sizeof(CMD) / sizeof(CMD[0]);
  for (i = 0; i < ncmds; i++) {
    char gap[] = "           ";
    gap[tsize - strlen(CMD[i].name)] = '\0';
    fprintf(stderr, "   %s%s%s\n", CMD[i].name, gap, CMD[i].desc);
  }
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
    pq_display_version();
    exit(0);
  }

  int i;
  int ncmds;
  int frm_multi;
  struct SWrap Stat;
  void (*swrap_init)(struct SWrap *, int);
  
  i = 0;
  ncmds = sizeof(CMD) / sizeof(CMD[0]);
  while (i < ncmds) {
    if (strcmp(argv[1], CMD[i].name) == 0) {
      frm_multi = CMD[i].frmt;
      swrap_init = CMD[i].init;
      break;
    }
    i++;
  }
  
  if (i == ncmds) {
    fprintf(stderr, "pqgen: '%s' is not a command. See 'pqgen --help'.\n", argv[1]);
    exit(0);
  }

  if (strcmp(argv[1], "pnds") == 0) {
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
  }
  
  int nargs;
  char **def_array;
  char defaults[128];
  
  int ncols;
  int nalleles;
  char buffer[PQ_LWIDTH];
  const char delim = '\t';

  char chr[PQ_LCOL];
  char factor[PQ_LCOL];
  unsigned long long int startpos;
  unsigned long long int stoppos;
  unsigned long long int start_region;
  unsigned long long int stop_region;

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
	sprintf(defaults, "-f 1 -c 1 -p 3 -k 5-%d -b 1", ncols);
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
      swrap_init(&Stat, nalleles);
      
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
      Stat.write(&Stat);
      printf("%s\t%llu\t%llu\t%s", chr, start_region, stop_region, factor);
      for (i = 0; i < Stat.nouts; i++) {
	printf("\t%s", (char *)Stat.outs[i]);
      }  
      printf("\n");
      Stat.clear(&Stat);      
      
      strcpy(factor, row.factor(&row));
      start_region = startpos;
      stop_region = stoppos;   
    }

    Stat.update(&Stat, row.array);
    
    stop_region = stoppos;
    strcpy(chr, row.chrom(&row));
  }

  Stat.write(&Stat);
  printf("%s\t%llu\t%llu\t%s", chr, start_region, stop_region, factor);
  for (i = 0; i < Stat.nouts; i++) {
    printf("\t%s", (char *)Stat.outs[i]);
  }
  printf("\n");
  Stat.clear(&Stat);

  rwk_free_hash(&CODON_TO_NSYN);
  rwk_free_hash(&CODON_TO_AMINO);
  pq_free_args();
  pq_swfree(&Stat);
  free_row(&row);
  
  return 0;
}
