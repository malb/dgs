#include <time.h>
#include <sys/time.h>
#include "dgs/dgs_gauss.h"

#define DP 0
#define MP 1

unsigned long long walltime(unsigned long long t0);

typedef struct {
  double sigma;
  long tau;
  double c;
  dgs_disc_gauss_alg_t algorithm;
  int precision;
  size_t ntrials;
} cmdline_params_gauss_z_t;


void parse_gauss_z_cmdline(cmdline_params_gauss_z_t *params, int argc, char *argv[]);
void print_gauss_z_help(const char *name);
