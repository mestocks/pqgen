#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <pq_htable.h>

#include <pq_args.h>
#include <pq_parse.h>
#include <pq_limits.h>
#include <pq_generics.h>
#include <pq_genetics.h>

#include <pq_het.h>
#include <pq_sfs.h>
#include <pq_theta.h>

#include <pq_help.h>
#include <pq_version.h>

struct pq_command {
  int frmt;
  const char *name;
  const char *desc;
  const char *defs;
  void (*init)(struct SWrap *, int);
};

struct pq_command CMD[] =
  {
    {.name = PQ_HET_NAME,
     .desc = PQ_HET_DESC,
     .frmt = PQ_HET_FRMT,
     .init = pq_het_init,
     .defs = PQ_HET_DEFS},

    {.name = PQ_SFS_NAME,
     .desc = PQ_SFS_DESC,
     .frmt = PQ_SFS_FRMT,
     .init = pq_sfs_init,
     .defs = PQ_SFS_DEFS},
    
    {.name = PQ_THETA_NAME,
     .desc = PQ_THETA_DESC,
     .frmt = PQ_THETA_FRMT,
     .init = pq_theta_init,
     .defs = PQ_THETA_DEFS}
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
  const char *cmd_defaults;
  struct SWrap Stat;
  
  void (*swrap_init)(struct SWrap *, int);
  
  i = 0;
  ncmds = sizeof(CMD) / sizeof(CMD[0]);
  while (i < ncmds) {
    if (strcmp(argv[1], CMD[i].name) == 0) {
      frm_multi = CMD[i].frmt;
      swrap_init = CMD[i].init;
      cmd_defaults = CMD[i].defs;
      break;
    }
    i++;
  }
  
  if (i == ncmds) {
    fprintf(stderr, "pqgen: '%s' is not a command. See 'pqgen --help'.\n", argv[1]);
    exit(0);
  }

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
  if (argc > 2 && access(argv[argc - 1], F_OK) != -1) {
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
  pq_update_args(argc_wo_file-2, argv+2);
  
  // space delimeters must be escaped on the command line (i.e. -d '\ ')
  delim = ((char *)pq_lookup_hash(&ARGHASH, "-d"))[1];
  ncols = pq_countcols(buffer, &delim);
  sprintf(defaults, cmd_defaults, ncols);
  pq_str2array(def_array, defaults, nargs, " ");
  
  pq_update_args(nargs, def_array);
  pq_update_args(argc_wo_file-2, argv+2);
  
  free(def_array);
  
  init_row(&row, ncols, CHROM, POS, FCOL);
  
  nalleles = frm_multi * NKARGS;
  swrap_init(&Stat, nalleles);

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

  pq_free_args();
  pq_swfree(&Stat);
  free_row(&row);
  fclose(fp);
  
  return 0;
}
