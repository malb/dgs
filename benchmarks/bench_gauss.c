#include <stdio.h>
#include <math.h>

#include "dgs.h"
#include "bench.h"

double run_dp(double sigma, double c, int tau, dgs_disc_gauss_alg_t alg, size_t ntrials, unsigned long long *t) {
  double variance = 0.0;
  gmp_randstate_t state;
  gmp_randinit_default(state);

  dgs_disc_gauss_dp_t *gen = dgs_disc_gauss_dp_init(sigma, c, tau, alg);

  *t =  walltime(0);
  for(size_t i=0; i<ntrials; i++) {
    long r = gen->call(gen);
    variance += ((double)r)*((double)r);
  }
  *t = walltime(*t);

  dgs_disc_gauss_dp_clear(gen);
  gmp_randclear(state);

  variance /= ntrials;
  return sqrt(variance);
}

double run_mp(double sigma_, double c_, int tau, dgs_disc_gauss_alg_t alg, size_t ntrials, unsigned long long *t) {
  mpfr_set_default_prec(80);
  mpfr_t sigma;
  mpfr_init_set_d(sigma, sigma_, MPFR_RNDN);
  gmp_randstate_t state;
  gmp_randinit_default(state);

  mpfr_t c;
  mpfr_init_set_d(c, c_, MPFR_RNDN);
  dgs_disc_gauss_mp_t *gen = dgs_disc_gauss_mp_init(sigma, c, tau, alg);

  double variance = 0.0;
  mpz_t r;
  mpz_init(r);
  
  *t =  walltime(0);
  for(size_t i=0; i<ntrials; i++) {
    gen->call(r, gen, state);
    variance += mpz_get_d(r)*mpz_get_d(r);
  }
  *t = walltime(*t);
  dgs_disc_gauss_mp_clear(gen);
  mpfr_clear(sigma);
  mpz_clear(r);
  mpfr_clear(c);

  gmp_randclear(state);
  
  variance /= ntrials;
  return sqrt(variance);
}

double run(double sigma, double c, int tau, int prec, dgs_disc_gauss_alg_t alg, size_t ntrials, unsigned long long *t) {
  if (prec == MP)
    return run_mp(sigma, c, tau, alg, ntrials, t);
  else 
    return run_dp(sigma, c, tau, alg, ntrials, t);
}


const char *alg_to_str(dgs_disc_gauss_alg_t alg) {
  switch(alg) {
  case DGS_DISC_GAUSS_UNIFORM_TABLE: return "uniform+table";
  case DGS_DISC_GAUSS_UNIFORM_ONLINE: return "uniform+online";
  case DGS_DISC_GAUSS_UNIFORM_LOGTABLE: return "uniform+logtable";
  case DGS_DISC_GAUSS_SIGMA2_LOGTABLE: return "sigma2+logtable";
  case DGS_DISC_GAUSS_ALIAS:           return "alias";
  default: return "unknown";
  }
}

int main(int argc, char *argv[]) {
  const double sigma2 = sqrt(1.0/(2.0*log(2.0)));
  
  unsigned long long t;

  cmdline_params_gauss_z_t params;
  params.sigma = 2 * sigma2;
  params.tau = 6;
  params.c = 0;
  params.ntrials = 10000000;
  params.algorithm = DGS_DISC_GAUSS_UNIFORM_TABLE;
  params.precision = MP;


  parse_gauss_z_cmdline(&params, argc, argv);

  if (params.algorithm == DGS_DISC_GAUSS_SIGMA2_LOGTABLE) {
    int k = round(params.sigma/sigma2);
    params.sigma = k*sigma2;
  }

  printf("%s :: σ: %.2f, c: %.2f. τ: %ld, precision: %d, algorithm: %d -- ", argv[0],
         params.sigma, params.c, params.tau, params.precision, params.algorithm);
  
  run(params.sigma, params.c, params.tau, params.precision, params.algorithm, params.ntrials, &t);
  double walltime = t/100000.0/params.ntrials*(1000.0*1000.0); // μs

  printf("wall time: %8.3f μs per call (rate: %8.3f per second)\n", walltime, 1000.0*1000.0/walltime);

  mpfr_free_cache();
  
  return 0;
}
