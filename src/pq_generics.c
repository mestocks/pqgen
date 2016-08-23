#include <stdio.h>
#include <stdlib.h>

#include <rwk_args.h>
#include <rwk_htable.h>
#include "pq_generics.h"

void pq_init_parameters(struct pq_parameters *params)
{
  rwk_create_hash(&params->arghash, 128);
  params->KCOLS = NULL;
  params->update = pq_update_parameters;
}

void pq_free_parameters(struct pq_parameters *params)
{
  free(params->KCOLS);
  rwk_free_hash(&params->arghash);
}

void pq_update_parameters(struct pq_parameters *params, unsigned int argc, char **argv)
{  
  char *chr_ptr;
  char *end_ptr;
  char *kcol_ptr;
  char *fcol_ptr;
  
  rwk_multi_insert_hash(&params->arghash, argc, argv);
  chr_ptr = (char *)rwk_lookup_hash(&params->arghash, "-c");
  end_ptr = (char *)rwk_lookup_hash(&params->arghash, "-p");
  kcol_ptr = (char *)rwk_lookup_hash(&params->arghash, "-k");
  fcol_ptr = (char *)rwk_lookup_hash(&params->arghash, "-f");
  
  params->CHROM = atoi(chr_ptr) - 1;
  params->POS = atoi(end_ptr) - 1;
  params->FCOL = atoi(fcol_ptr) - 1;
  
  if (params->KCOLS != NULL) {
    free(params->KCOLS);
  }

  params->nkargs = rwk_countkargs(kcol_ptr);
  params->KCOLS = malloc(params->nkargs * sizeof(int));
  rwk_kargsplit(params->KCOLS, kcol_ptr);
}

void pq_print_bed(struct Data *data)
{
  printf("%s\t%llu\t%llu\t%s\t%llu\t%s\n", data->chr, *data->start, *data->end, data->name, *data->nvsites, data->value);
}

void pq_print_csv(struct Data *data)
{
  printf("%s,%llu,%llu,%s,%llu,%s\n", data->chr, *data->start, *data->end, data->name, *data->nvsites, data->value);
}
