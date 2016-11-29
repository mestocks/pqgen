#include <stdio.h>

#include <pq_version.h>
#include <rwk_version.h>

void pq_display_version()
{
  fprintf(stderr, "pqgen version %s (librawk version %s)\n", PQ_VERSION, RWK_VERSION);
}
