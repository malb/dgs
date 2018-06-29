#include <math.h>
#include "dgs/dgs_rround.h"

#define NTRIALS 1<<18
#define TOLERANCE 0.1
#define BOUND 2

int test_defaults_dp() {
  dgs_rround_dp_t *self;
  self = dgs_rround_dp_init(6, DGS_RROUND_DEFAULT);
  if (self->algorithm != DGS_RROUND_KARNEY)
    dgs_die("automatic choice of Karney's algorithm failed (%d)", self->algorithm);
  dgs_rround_dp_clear(self);

  printf("passed\n");
  return 0;
}

int test_uniform_boundaries_dp(double sigma, double c, size_t tau, dgs_rround_alg_t algorithm) {
  dgs_rround_dp_t *self = dgs_rround_dp_init(tau, algorithm);

    printf("σ: %6.2f, c: %6.2f. τ: %2ld, precision: double, algorithm: %d\n", sigma, c, self->tau, self->algorithm);


  long lower_bound = ((long)c) - ceil(sigma * self->tau);
  long upper_bound = ((long)c) + ceil(sigma * self->tau);

  for(size_t i=0; i<NTRIALS; i++) {
    long r = self->call(self, sigma, c);
    if(__DGS_UNLIKELY(r < lower_bound))
      dgs_die("r (%ld) < lower_bound (%ld)", r, lower_bound);
    else if(__DGS_UNLIKELY(r > upper_bound))
      dgs_die("r (%ld) > upper_bound (%ld)", r, upper_bound);
  }
  return 0;
}

/**
 We test if the proportional probabilities holds
*/
#define RHO(x) exp(-(x*x)/(2*sigma*sigma))

int test_ratios_dp(double sigma, size_t tau, dgs_rround_alg_t algorithm) {
  printf("σ: %6.2f, c:    0.0. τ: %2ld, precision: double, algorithm: %d\n", sigma, tau, algorithm);

  dgs_rround_dp_t *self = dgs_rround_dp_init(tau, algorithm);

  double ctr[2*BOUND+1];

  for(size_t i=0; i<NTRIALS; i++) {
    long r = self->call(self, sigma, 0);
    if (abs(r) <= BOUND)
      ctr[r+BOUND] += 1;
  }

  for(long i=-BOUND; i<=BOUND; i++) {
    double left  = ctr[BOUND+1]/ctr[BOUND+i];
    double right = RHO(0)/RHO(i);
    if (fabs(log(left/right)) >= 0.4)
      dgs_die("exp(-((-c)²)/(2σ²))/exp(-((%d-c)²)/(2σ²)) = %7.5f != %7.5f (%7.5f)", i, right, left, fabs(log(left/right)));
  }
  return 0;
}

int test_mean_dp(double sigma, double c, size_t tau, dgs_rround_alg_t algorithm) {

  printf("σ: %6.2f, c: %6.2f. τ: %2ld, precision: double, algorithm: %d\n",sigma, c, tau, algorithm);

  dgs_rround_dp_t *self = dgs_rround_dp_init(tau, algorithm);

  double mean = 0.0;

  for(size_t i=0; i<NTRIALS; i++) {
    long r = self->call(self, sigma, c);
    mean += r;

  }

  mean /=NTRIALS;

  if(fabs(mean - c) > TOLERANCE)
    dgs_die("expected mean %6.2f but got %6.2f",c, mean);

  return 0;
}


int main(int argc, char *argv[]) {
  printf("# testing defaults #\n");
  test_defaults_dp();
  printf("\n");

  printf("# testing proportional probabilities #\n");
  test_ratios_dp( 3.0, 6, DGS_RROUND_UNIFORM_ONLINE);
  test_ratios_dp( 2.0, 6, DGS_RROUND_UNIFORM_ONLINE);
  test_ratios_dp( 4.0, 3, DGS_RROUND_UNIFORM_ONLINE);
  test_ratios_dp(15.4, 3, DGS_RROUND_UNIFORM_ONLINE);
  printf("\n");

  test_ratios_dp( 3.0, 6, DGS_RROUND_KARNEY);
  test_ratios_dp( 2.0, 6, DGS_RROUND_KARNEY);
  test_ratios_dp( 4.0, 3, DGS_RROUND_KARNEY);
  test_ratios_dp(15.4, 3, DGS_RROUND_KARNEY);
  printf("\n");

  test_ratios_dp( 3.0, 6, DGS_RROUND_CONVOLUTION);
  test_ratios_dp( 2.0, 6, DGS_RROUND_CONVOLUTION);
  test_ratios_dp( 4.0, 3, DGS_RROUND_CONVOLUTION);
  test_ratios_dp(15.4, 3, DGS_RROUND_CONVOLUTION);
  printf("\n");

  printf("# testing [⌊c⌋-⌈στ⌉,…, ⌊c⌋+⌈στ⌉] boundaries #\n");
  test_uniform_boundaries_dp( 3.0, 0.0, 2, DGS_RROUND_UNIFORM_ONLINE);
  test_uniform_boundaries_dp(10.0, 0.0, 2, DGS_RROUND_UNIFORM_ONLINE);
  test_uniform_boundaries_dp( 3.3, 1.0, 1, DGS_RROUND_UNIFORM_ONLINE);
  test_uniform_boundaries_dp( 2.0, 1.5, 2, DGS_RROUND_UNIFORM_ONLINE);
  printf("\n");

  printf("# testing c is center #\n");
  test_mean_dp( 3.0, 0.0, 6, DGS_RROUND_UNIFORM_ONLINE);
  test_mean_dp(10.0, 0.0, 6, DGS_RROUND_UNIFORM_ONLINE);
  test_mean_dp( 3.3, 1.0, 6, DGS_RROUND_UNIFORM_ONLINE);
  test_mean_dp( 2.0, 1.5, 6, DGS_RROUND_UNIFORM_ONLINE);
  printf("\n");

  test_mean_dp( 3.0, 0.0, 6, DGS_RROUND_KARNEY);
  test_mean_dp(10.0, 0.0, 6, DGS_RROUND_KARNEY);
  test_mean_dp( 3.3, 1.0, 6, DGS_RROUND_KARNEY);
  test_mean_dp( 2.0, 1.5, 6, DGS_RROUND_KARNEY);
  printf("\n");

  test_mean_dp( 3.0, 0.0, 6, DGS_RROUND_CONVOLUTION);
  test_mean_dp(10.0, 0.0, 6, DGS_RROUND_CONVOLUTION);
  test_mean_dp( 3.3, 1.0, 6, DGS_RROUND_CONVOLUTION);
  test_mean_dp( 2.0, 1.5, 6, DGS_RROUND_CONVOLUTION);
  printf("\n");

  return 0;
}
