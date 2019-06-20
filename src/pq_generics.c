#include <stdio.h>
#include <stdlib.h>

#include "pq_generics.h"


void pq_print_bed(struct Data *data)
{
  printf("%s\t%llu\t%llu\t%s\t%llu\t%s\n", data->chr, *data->start, *data->end, data->name, *data->nvsites, data->value);
}

void pq_print_csv(struct Data *data)
{
  printf("%s,%llu,%llu,%s,%llu,%s\n", data->chr, *data->start, *data->end, data->name, *data->nvsites, data->value);
}


void pq_swfree(struct StatObject *stats)
{
  int i;
  for (i = 0; i < stats->nouts; i++) {
    free(stats->outs[i]);
  }
  for (i = 0; i < stats->nvalues; i++) {
    free(stats->values[i]);
  }
  free(stats->outs);
  free(stats->values);
}


