#include <stdio.h>
#include <math.h>

#include "dgs.h"
#include "bench.h"

double run_dp(double sigma, double c, int tau, dgs_rround_alg_t alg, size_t ntrials, unsigned long long *t) {
  double variance = 0.0;
  gmp_randstate_t state;
  gmp_randinit_default(state);

  dgs_rround_dp_t *gen = dgs_rround_dp_init(tau, alg);

  *t =  walltime(0);
  for(size_t i=0; i<ntrials; i++) {
    long r = gen->call(gen, sigma, c);
    variance += ((double)r)*((double)r);
  }
  *t = walltime(*t);

  dgs_rround_dp_clear(gen);
  gmp_randclear(state);

  variance /= ntrials;
  return sqrt(variance);
}

double run_mp(double sigma_, double c_, int tau, dgs_rround_alg_t alg, size_t ntrials, unsigned long long *t) {
  mpfr_set_default_prec(80);
  mpfr_t sigma;
  mpfr_init_set_d(sigma, sigma_, MPFR_RNDN);
  gmp_randstate_t state;
  gmp_randinit_default(state);

  mpfr_t c;
  mpfr_init_set_d(c, c_, MPFR_RNDN);
  dgs_rround_mp_t *gen = dgs_rround_mp_init(tau, alg, mpfr_get_default_prec());

  double variance = 0.0;
  mpz_t r;
  mpz_init(r);
  
  *t =  walltime(0);
  for(size_t i=0; i<ntrials; i++) {
    gen->call(r, gen, sigma, c, state);
    variance += mpz_get_d(r)*mpz_get_d(r);
  }
  *t = walltime(*t);
  dgs_rround_mp_clear(gen);
  mpfr_clear(sigma);
  mpz_clear(r);
  mpfr_clear(c);

  gmp_randclear(state);
  
  variance /= ntrials;
  return sqrt(variance);
}

double run(double sigma, double c, int tau, int prec, dgs_rround_alg_t alg, size_t ntrials, unsigned long long *t) {
  if (prec == MP)
    return run_mp(sigma, c, tau, alg, ntrials, t);
  else 
    return run_dp(sigma, c, tau, alg, ntrials, t);
}


const char *alg_to_str(dgs_disc_gauss_alg_t alg) {
  switch(alg) {
  case DGS_RROUND_UNIFORM_ONLINE: return "uniform+online";
  case DGS_RROUND_KARNEY: return "karney";
  default: return "unknown";
  }
}

int main(int argc, char *argv[]) {
  
  unsigned long long t;

  cmdline_params_rround_z_t params;
  params.sigma = 3.0;
  params.tau = 6;
  params.c = 0;
  params.ntrials = 10000000;
  params.algorithm = DGS_RROUND_DEFAULT;
  params.precision = MP;


  parse_rround_z_cmdline(&params, argc, argv);

  printf("%s :: σ: %.2f, c: %.2f. τ: %ld, precision: %d, algorithm: %d -- ", argv[0],
         params.sigma, params.c, params.tau, params.precision, params.algorithm);
  
  run(params.sigma, params.c, params.tau, params.precision, params.algorithm, params.ntrials, &t);
  double walltime = t/100000.0/params.ntrials*(1000.0*1000.0); // μs

  printf("wall time: %8.3f μs per call (rate: %8.3f per second)\n", walltime, 1000.0*1000.0/walltime);

  mpfr_free_cache();
  
  return 0;
}
