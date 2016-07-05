
#include <stdio.h>

#include <pq_version.h>
#include <rwk_version.h>

int main(int argc, char **argv) {

  printf("pq-genetics-v%s\nusing librawk-v%s\n", PQ_VERSION, RWK_VERSION);
  
  return 0; }
