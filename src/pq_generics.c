#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_args.h>
#include <rwk_htable.h>
#include <pq_generics.h>


void pq_print_bed(struct Data *data)
{
  printf("%s\t%llu\t%llu\t%s\t%llu\t%s\n", data->chr, *data->start, *data->end, data->name, *data->nvsites, data->value);
}

void pq_print_csv(struct Data *data)
{
  printf("%s,%llu,%llu,%s,%llu,%s\n", data->chr, *data->start, *data->end, data->name, *data->nvsites, data->value);
}


void pq_swfree(struct SWrap *wrap)
{
  int i;
  for (i = 0; i < wrap->nouts; i++) {
    free(wrap->outs[i]);
  }
  for (i = 0; i < wrap->nvalues; i++) {
    free(wrap->values[i]);
  }
  free(wrap->outs);
  free(wrap->values);
}


