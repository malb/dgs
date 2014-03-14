#include <math.h>
#include "dgs/dgs_gauss.h"

#define NTRIALS 100000
#define TOLERANCE 0.1

int test_uniform_boundaries_dp(double sigma, double c, size_t tau, dgs_disc_gauss_alg_t algorithm) {

  printf("σ: %6.2f, c: %6.2f. τ: %2ld, precision: double, algorithm: %d\n",sigma, c, tau, algorithm);

  dgs_disc_gauss_dp_t *self = dgs_disc_gauss_dp_init(sigma, c, tau, algorithm);
  
  long lower_bound = ((long)self->c) - ceil(self->sigma * self->tau);
  long upper_bound = ((long)self->c) + ceil(self->sigma * self->tau);
  
  for(size_t i=0; i<NTRIALS; i++) {
    long r = self->call(self);
    if(__DGS_UNLIKELY(r < lower_bound))
      dgs_die("r (%ld) < lower_bound (%ld)", r, lower_bound);
    else if(__DGS_UNLIKELY(r > upper_bound))
      dgs_die("r (%ld) > upper_bound (%ld)", r, upper_bound);
  }
  return 0;
}

int test_mean_dp(double sigma, double c, size_t tau, dgs_disc_gauss_alg_t algorithm) {

  printf("σ: %6.2f, c: %6.2f. τ: %2ld, precision: double, algorithm: %d\n",sigma, c, tau, algorithm);

  dgs_disc_gauss_dp_t *self = dgs_disc_gauss_dp_init(sigma, c, tau, algorithm);

  double mean = 0.0;
  
  for(size_t i=0; i<NTRIALS; i++) {
    long r = self->call(self);
    mean += r;
    
  }

  mean /=NTRIALS;

  if(fabs(mean - c) > TOLERANCE)
    dgs_die("expected mean %6.2f but got %6.2f",c, mean);
  
  return 0;
}


int main(int argc, char *argv[]) {
  printf("# testing [⌊c⌋-⌈στ⌉,…, ⌊c⌋+⌈στ⌉] boundaries #\n");
  test_uniform_boundaries_dp( 3.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_ONLINE);  
  test_uniform_boundaries_dp(10.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_ONLINE);  
  test_uniform_boundaries_dp( 3.3, 1.0, 1, DGS_DISC_GAUSS_UNIFORM_ONLINE);  
  test_uniform_boundaries_dp( 2.0, 1.5, 2, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  printf("\n");

  test_uniform_boundaries_dp( 3.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_TABLE);  
  test_uniform_boundaries_dp(10.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_TABLE);  
  test_uniform_boundaries_dp( 3.3, 1.0, 1, DGS_DISC_GAUSS_UNIFORM_TABLE);  
  test_uniform_boundaries_dp( 2.0, 1.5, 2, DGS_DISC_GAUSS_UNIFORM_TABLE);  
  printf("\n");

  test_uniform_boundaries_dp( 3.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);  
  test_uniform_boundaries_dp(10.0, 0.0, 2, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);  
  test_uniform_boundaries_dp( 3.3, 1.0, 1, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);  
  test_uniform_boundaries_dp( 2.0, 2.0, 2, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  printf("\n");

  printf("# testing c is center #\n");
  test_mean_dp( 3.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_ONLINE);  
  test_mean_dp(10.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_ONLINE);  
  test_mean_dp( 3.3, 1.0, 6, DGS_DISC_GAUSS_UNIFORM_ONLINE);  
  test_mean_dp( 2.0, 1.5, 6, DGS_DISC_GAUSS_UNIFORM_ONLINE);
  printf("\n");

  test_mean_dp( 3.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_TABLE);  
  test_mean_dp(10.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_TABLE);  
  test_mean_dp( 3.3, 1.0, 6, DGS_DISC_GAUSS_UNIFORM_TABLE);  
  test_mean_dp( 2.0, 1.5, 6, DGS_DISC_GAUSS_UNIFORM_TABLE);  
  printf("\n");

  test_mean_dp( 3.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);  
  test_mean_dp(10.0, 0.0, 6, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);  
  test_mean_dp( 3.3, 1.0, 6, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);  
  test_mean_dp( 2.0, 2.0, 6, DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  printf("\n");

  test_mean_dp( 3.0, 0.0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);  
  test_mean_dp(10.0, 0.0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);  
  test_mean_dp( 3.3, 1.0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);  
  test_mean_dp( 2.0, 2.0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);

  printf("\n");

  return 0;
}
