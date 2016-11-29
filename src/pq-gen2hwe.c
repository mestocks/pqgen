#include <stdio.h>
#include <stdlib.h>

#include <rwk_parse.h>

#define COLOR "\x1B[33m"
#define NORMAL "\x1B[0m"

double SNPHWE_pValue(int obs_hets, int obs_hom1, int obs_hom2)
			    //	double *het_probs)
{
  int obs_homc = (obs_hom1 < obs_hom2) ? obs_hom2 : obs_hom1;
  int obs_homr = (obs_hom1 < obs_hom2) ? obs_hom1 : obs_hom2;
  
  int rare_copies = 2 * obs_homr + obs_hets;
  int genotypes   = obs_hets + obs_homc + obs_homr;
  
  double *het_probs = malloc((rare_copies + 1) * sizeof(double));
  
  // start at midpoint
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
  
  // check to ensure that midpoint and rare alleles have same parity
  if ((rare_copies & 1) ^ (mid & 1)) {
    mid ++;
  }

  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;
  
  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
    het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets *
      (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
    sum += het_probs[curr_hets - 2];
    
    // 2 fewer heterozygotes for next iteration -> add one rare,
    // one common homozygote
    curr_homr++;
    curr_homc++;
  }
  
  
  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
    het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr *
      curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0));
    sum += het_probs[curr_hets + 2];
    
    // add 2 heterozygotes for next iteration -> subtract one rare,
    // one common homozygote
    curr_homr--;
    curr_homc--;
  }
  
  
  for (int i=0; i <= rare_copies; i++) {
    het_probs[i] /= sum;
  }
  
  double p_hwe = 0.0;
  // p-value calculation for p_hwe
  for (int i=0; i <= rare_copies; i++) {
    if (het_probs[i] > het_probs[obs_hets])
      continue;
    p_hwe += het_probs[i];
  }
  
  return (p_hwe > 1.0) ? 1.0 : p_hwe;
}

int old_main(int argc, char **argv)
{
  int ohets;
  int ohom1;
  int ohom2;
  double pval;

  // 1st and 4th columns from Table 1 in Wigginton et al., 2005
  for (ohets = 21; ohets <= 21; ohets += 2) {
    ohom1 = (21 - ohets) / 2;
    ohom2 = 100 - (ohets + ohom1);
    pval = SNPHWE_pValue(ohets, ohom1, ohom2);
    printf("%d\t%f\n", ohets, pval);
  }
  
  return 0;
}

#define LWIDTH 5120

int main(int argc, char **argv)
{
  int ohets;
  int ohom1;
  int ohom2;
  double pval;

  int i;
  int ncols;
  char **array;
  char buffer[LWIDTH];
  const char delim = '\t';

  int rown;
  rown = 0;
  
  while (fgets(buffer, sizeof(buffer), stdin)) {

    if (rown == 0) {
      ncols = rwk_countcols(buffer, &delim);
      array = calloc(ncols, sizeof (char *));
    }

    if (rwk_str2array(array, buffer, ncols, &delim) == -1) {
      free(array);
      exit(1);
    }

    ohom1 = atoi(array[4]);
    ohets = atoi(array[5]);
    ohom2 = atoi(array[6]);
    pval = SNPHWE_pValue(ohets, ohom1, ohom2);
    
    printf("%s", array[0]);
    for (i = 1; i < 7; i++) {
      printf("\t%s", array[i]);
    }
    printf("\t%f\n", pval);
  }
  
  return 0;
}
