#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "pq_htable.h"

#include "pq_args.h"
#include "pq_parse.h"
#include "pq_limits.h"
#include "pq_generics.h"
#include "pq_genetics.h"
#include "pq_pnds.h"
#include "pq_version.h"


int main(int argc, char **argv)
{
  if (argc > 1 && strcmp(argv[1], "--help") == 0) {
    printf("Help!\n");
    exit(0);
  }
  
  if (argc > 1 && strcmp(argv[1], "--version") == 0) {
    pq_display_version();
    exit(0);
  }

  int i;
  int frm_multi;
  const char *cmd_defaults;
  struct StatObject Stat;
  
  void (*stat_init)(struct StatObject *, int);

  frm_multi = PQ_PNDS_FRMT;
  stat_init = PQ_PNDS_INIT;
  cmd_defaults = PQ_PNDS_DEFS;
  
  char *fname_aa;
  char *fname_syn;
  char default_fname_aa[] = "codon2aa";
  char default_fname_syn[] = "codon2syn";
  
  char default_fullpath_aa[1028];
  char default_fullpath_syn[1028];
  sprintf(default_fullpath_aa, "%s/%s", DATADIR, default_fname_aa);
  sprintf(default_fullpath_syn, "%s/%s", DATADIR, default_fname_syn);
  
  fname_aa = default_fullpath_aa;
  fname_syn = default_fullpath_syn;
  file2charHash(fname_aa, 128);
  file2doubleHash(fname_syn, 128);
  
  FILE *fp;
  
  int nargs;
  char **def_array;
  char defaults[128];
  
  int ncols;
  int nalleles;
  char buffer[PQ_LWIDTH];
  char delim = '\t';

  char chr[PQ_LCOL];
  char factor[PQ_LCOL];
  unsigned long long int startpos;

  unsigned long long int stoppos;
  unsigned long long int start_region;
  unsigned long long int stop_region;
  
  struct GenericRow row;

  int argc_wo_file;
  
  fp = stdin;
  argc_wo_file = argc;
  if (argc > 1 && access(argv[argc - 1], F_OK) != -1) {
    fp = fopen(argv[argc - 1], "r");
    argc_wo_file = argc - 1;
  }

  fgets(buffer, sizeof(buffer), fp);

  // add defaults but with dummy -k variable
  sprintf(defaults, cmd_defaults, 6);
  
  nargs = pq_countcols(defaults, " ");
  def_array = calloc(nargs, sizeof(char *));
  pq_str2array(def_array, defaults, nargs, " ");
  
  pq_init_args();
  pq_update_args(nargs, def_array);
  pq_update_args(argc_wo_file-1, argv+1);
  
  // space delimeters must be escaped on the command line (i.e. -d '\ ')
  delim = ((char *)pq_lookup_hash(&ARGHASH, "-d"))[1];
  ncols = pq_countcols(buffer, &delim);
  sprintf(defaults, cmd_defaults, ncols);
  pq_str2array(def_array, defaults, nargs, " ");
  
  pq_update_args(nargs, def_array);
  pq_update_args(argc_wo_file-1, argv+1);
  
  free(def_array);
  
  init_row(&row, ncols, CHROM, POS, FCOL);
  
  nalleles = frm_multi * NKARGS;
  stat_init(&Stat, nalleles);

  row.update(&row, buffer, &delim);
  strcpy(chr, row.chrom(&row));
  strcpy(factor, row.factor(&row));
  start_region = row.pos(&row) - 1;
  
  goto POST_INIT;
  
  while (fgets(buffer, sizeof(buffer), fp)) {
    
    row.update(&row, buffer, &delim);

  POST_INIT:
    
    stoppos = row.pos(&row);
    startpos = stoppos - 1;
    
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

  pq_free_hash(&CODON_TO_NSYN);
  pq_free_hash(&CODON_TO_AMINO);
  pq_free_args();
  pq_swfree(&Stat);
  free_row(&row);
  fclose(fp);
  
  return 0;
}
