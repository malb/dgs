#include <unistd.h>
#include "bench.h"

unsigned long long walltime(unsigned long long t0) {
  static time_t base_sec;
  struct timeval tp;
  gettimeofday(&tp, NULL);
  if (__DGS_UNLIKELY(base_sec == 0))
    base_sec = tp.tv_sec;
  return (tp.tv_sec - base_sec) * 1000000 + tp.tv_usec - t0;
}

void print_gauss_z_help(const char *name) {
  //printf("%s:\n\n", name);
  printf("REQUIRED:\n");
  printf(" s -- Gaussian width parameter σ, double > 0\n");
  printf(" t -- cutoff parameter for sampling from uniform distribution,\n");
  printf("      values outside of [⌊c⌋-⌈στ⌉, ⌊c⌋+⌈στ⌉] are considered to have probability zero.\n");
  printf("      integer > 0\n");
  printf(" c -- center of Gaussian distribution. double\n");
  printf(" a -- algorithm:\n");
  printf("      %d -- sample from uniform distribution, call exp()\n", DGS_DISC_GAUSS_UNIFORM_ONLINE);
  printf("      %d -- sample from uniform distribution, exp() calls tabulated\n", DGS_DISC_GAUSS_UNIFORM_TABLE);
  printf("      %d -- sample from uniform distribution, exp() calls as Bernoulli oracles\n", DGS_DISC_GAUSS_UNIFORM_LOGTABLE);
  printf("      %d -- sample from k⋅σ2 distribution, exp() calls as Bernoulli oracles \n", DGS_DISC_GAUSS_SIGMA2_LOGTABLE);
  printf(" p -- precision: 0 for double precision, 1 for arbitrary precision\n");
  printf(" n -- number of trials > 0 (default: 100000)\n");
}

void parse_gauss_z_cmdline(cmdline_params_gauss_z_t *params, int argc, char *argv[]) {
  int c;
  while ((c = getopt(argc, argv, "s:t:c:a:p:hn:")) != -1)
    switch (c) {
    case 's':
      params->sigma = atof(optarg); break;
    case 't':
      params->tau = atol(optarg); break;
    case 'c':
      params->c = atof(optarg); break;
    case 'a':
      params->algorithm = atoi(optarg); break;
    case 'p':
      params->precision = atoi(optarg); break;
    case 'n':
      params->ntrials = atoi(optarg); break;
    case 'h':
      print_gauss_z_help(argv[0]);
      exit(0);
    default:
      print_gauss_z_help(argv[0]);
      abort();
    }
  if (params->sigma <= 0.0)
    dgs_die("σ > 0 required, but got σ = %f",params->sigma);
  if (params->tau <= 0)
    dgs_die("τ > 0 required, but got τ = %d",params->tau);
  if (params->precision != 0 && params->precision != 1)
    dgs_die("precision must be either 0 or 1, but got %d",params->precision);
}
