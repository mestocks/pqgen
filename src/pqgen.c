#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <pq_version.h>
#include <rwk_version.h>

int main(int argc, char **argv) {

  char *path;
  char *home;
  char usage[] = "usage: pqgen [--version] [--help] <command> [<args>]\n";
  char commands[] = "  theta         Calculate site frequency based stats\n  dna2codon     Convert nucleotide sequences into codons\n  codon2pnds    Counts synonymous and non-synonymous sites from codons\n";
  char prefix[] = "/.local/bin/pq-";
  
  if (argc == 1 || (strcmp(argv[1], "--help") == 0)) {
    printf("%s\n", usage);
    printf("%s\n", commands);
    exit(0);
    
  } else if (strcmp(argv[1], "--version") == 0) {
    printf("pq-genetics version %s (linked to librawk version %s)\n", PQ_VERSION, RWK_VERSION);
    exit(0);
    
  } else {
    home = getenv("HOME");
    path = malloc(strlen(home) + strlen(prefix) + strlen(argv[1]) + 1);
    sprintf(path, "%s%s%s", home, prefix, argv[1]);

    if (access(path, F_OK) != -1) {
      execv(path, argv + 1);
      
    } else {
      printf("%s command does not exist\n\n", argv[1]);
      printf("%s\n", usage);
      printf("%s\n", commands);
      exit(0);
    }
  }
    
  return 0; }
