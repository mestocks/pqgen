#include <stdio.h>

#include "pq_version.h"

void pq_display_version()
{
  fprintf(stderr, "pqgen version %s\n", PQ_VERSION);
}
