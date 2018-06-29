/**
   Discrete Gaussian Rounding over the Integers.

   A discrete Gaussian distribution on the Integers is a distribution where the
   integer `x` is sampled with probability proportional to `exp(-(x-c)²/(2σ²))`.
   It is denoted by `D_{σ,c}` where `σ` is the width parameter (close to the
   standard deviation) and `c` is the center. This file contains algorithms
   suitable for situations, where for each set of parameters (σ,c) only a small
   number of samples are desired. I.e. no precomputation is done during
   initialization that depends on the parameters of the distribution. Drawing
   samples from such a distribution can be viewed as rounding `c` to a nearby
   integer, where 'nearby' is defined by the width parameter `σ`.

   AVAILABLE ALGORITHMS:

 - ``DGS_RROUND_UNIFORM_ONLINE`` - samples are drawn from a uniform
   distribution and accepted with probability proportional to
   `\exp(-(x-c)²/(2σ²))` where `\exp(-(x-c)²/(2σ²))` is computed in each
   invocation. Typically this is very slow.

 - ``DGS_RROUND_KARNEY`` - Use Karney's algorithm. This is better than
   uniform rejection sampling.

 - ``DGS_RROUND_CONVOLUTION`` - Use convolution to reduce to alias
   sampling.

  AVAILABLE PRECISIONS:

  - ``mp`` - multi-precision using MPFR, cf. ``dgs_gauss_mp.c``

  - ``dp`` - double precision using machine doubles, cf. ``dgs_gauss_dp.c``.

  For readers unfamiliar with the implemented algorithms it makes sense to start
  with ``dgs_gauss_dp.c`` which implements the same algorithms as
  ``dgs_gauss_mp.c`` should be easier to read.

  TYPICAL USAGE::

      dgs_rround_dp_t *D = dgs_rround_dp_init(<tau>, <algorithm>);
      D->call(D, <sigma>, <c>); // as often as needed
      dgs_rround_dp_clear(D);

   .. author:: Michael Walter

*/

#ifndef DGS_RROUND__H
#define DGS_RROUND__H

#include "dgs_gauss.h"

#define DGS_RROUND_SIGMA_LOG2_MAX 30
#define DGS_RROUND_SIGMA_MAX (1 << DGS_RROUND_SIGMA_LOG2_MAX)

/**
   Available Algorithms
*/

typedef enum {
  DGS_RROUND_DEFAULT        = 0x0,  //<pick algorithm
  DGS_RROUND_UNIFORM_ONLINE = 0x1,  //<call dgs_disc_gauss_mp_call_uniform_online
  DGS_RROUND_KARNEY         = 0x2,  //<call dgs_disc_gauss_mp_call_karney
  DGS_RROUND_CONVOLUTION    = 0x3,  //<call dgs_disc_gauss_mp_call_convolution
} dgs_rround_alg_t;

struct _dgs_rround_mp_t;

typedef struct _dgs_rround_dp_t {

  /**
     Cutoff `τ`, samples outside the range `(⌊c⌉-⌈στ⌉,...,⌊c⌉+⌈στ⌉)` are
     considered to have probability zero. This bound applies to algorithms
     which sample from the uniform distribution.
  */

  size_t tau;

  dgs_bern_uniform_t *B;
  dgs_bern_dp_t *B_half_exp;

  dgs_disc_gauss_dp_t **base_samplers;
  int log_base, digits, flips;
  double s_bar2, bm_sample;
  uint64_t mask, pool;

  dgs_rround_alg_t algorithm;  //<  which algorithm to use

  /**
   Return a ``long`` sampled from this rounder

   :param self: discrete Gaussian rounder.
   :param sigma: noise parameter.
   :param c: center.

  */

  long (*call)(struct _dgs_rround_dp_t *self, double sigma, double c);

} dgs_rround_dp_t;

/**
 Create a new double-precision discrete Gaussian rounder.

 :param tau: cutoff `τ`
 :param algorithm: algorithm to use.

*/

dgs_rround_dp_t *dgs_rround_dp_init(size_t tau, dgs_rround_alg_t algorithm);

/**
   Sample from ``dgs_rround_dp_t`` by rejection sampling using the uniform distribution

   :param self: Discrete Gaussian rounder
   :param sigma: noise parameter
   :param c: center

 */

long dgs_rround_dp_call_uniform_online(dgs_rround_dp_t *self, double sigma, double c);

/**
 Sample from ``dgs_rround_dp_t`` using Karney's algorithm. This is similar
 to the sigma2 method, but can sample from any c and σ.

 :param self: discrete Gaussian sampler
*/
long dgs_rround_dp_call_karney(dgs_rround_dp_t *self, double sigma, double c);

/**
 Sample from ``dgs_rround_dp_t`` using the convolution sampler

 :param self: discrete Gaussian sampler
*/
long dgs_rround_dp_call_convolution(dgs_rround_dp_t *self, double sigma, double c);

/**
   Free memory.

   :param self: discrete Gaussian rounder

 */
void dgs_rround_dp_clear(dgs_rround_dp_t *self);

/**
   Multi-precision Discrete Gaussians `D_{σ,c}`

   Return integer `x` with probability

   `ρ_{σ,c}(x) = exp(-(x-c)²/(2σ²))/exp(-(\ZZ-c)²/(2σ²))`

   where `exp(-(\ZZ-c)²/(2σ²)) ≈ \sum_{i=-τσ}^{τσ} exp(-(i-c)²/(2σ^²))` is the
   probability for all of the integers.

*/

typedef struct _dgs_rround_mp_t {

  dgs_bern_uniform_t *B;
  dgs_bern_mp_t *B_half_exp;

  /**
     Cutoff `τ`, samples outside the range `(⌊c⌉-⌈στ⌉,...,⌊c⌉+⌈στ⌉)` are
     considered to have probability zero. This bound applies to algorithms
     which sample from the uniform distribution.
  */
  size_t tau;

  dgs_rround_alg_t algorithm;  //<  which algorithm to use

  /**
   Return an ``mpz_t`` sampled from this sampler

   :param rop: target value.
   :param self: discrete Gaussian sampler.
   :param state: entropy pool.

  */

  void (*call)(mpz_t rop, struct _dgs_rround_mp_t *self, const mpfr_t sigma, const mpfr_t c, gmp_randstate_t state);

  /**
   * Temporary variables:
   */
  mpfr_t c_r;  //< `c_r := c % 1`
  mpz_t c_z;   //< c_z := c - (c_r)

  mpz_t sigma_z;

  /**
   We sample ``x`` with ``abs(x) < upper_bound`` in
   ``DGS_RROUND_UNIFORM_ONLINE``.
   */

  mpz_t upper_bound;

  /**
   We sample ``x`` with ``abs(x) <= upper_bound - 1`` in
   ``DGS_RROUND_UNIFORM_ONLINE``.
   */

  mpz_t upper_bound_minus_one;

  /**
     There are ``2*upper_bound -1`` elements in the range
     ``-upper_bound+1,...,upper_bound-1``.
   */

  mpz_t two_upper_bound_minus_one;

  /**
   Precomputed `-1/(2σ²)`.
  */

  mpfr_t f;

  mpz_t x;   //< space for temporary integer
  mpfr_t y;  // space for temporary rational number
  mpfr_t z;  // space for temporary rational number
  mpfr_t tmp;

  dgs_disc_gauss_mp_t *wide_sampler;
  dgs_disc_gauss_mp_t **base_samplers;
  int log_base, digits, flips, pool;
  mpfr_t s_bar2, bm_sample;
  uint64_t mask;
} dgs_rround_mp_t;

dgs_rround_mp_t *dgs_rround_mp_init(size_t tau, dgs_rround_alg_t algorithm, mpfr_prec_t prec);

/**
  Sample from ``dgs_rround_mp_t`` by rejection sampling using the uniform distribution.

  :param rop: return value
  :param self: discrete Gaussian rounder
  :param sigma: noise parameter
  :param c: center
  :param state: state

 */

void dgs_rround_mp_call_uniform_online(mpz_t rop, dgs_rround_mp_t *self, const mpfr_t sigma, const mpfr_t c,
                                       gmp_randstate_t state);

/**
  Sample from ``dgs_rround_mp_t`` using Karney's algorithm.

  :param rop: return value
  :param self: discrete Gaussian rounder
  :param sigma: noise parameter
  :param c: center
  :param state: state

 */

void dgs_rround_mp_call_karney(mpz_t rop, dgs_rround_mp_t *self, const mpfr_t sigma, const mpfr_t c,
                               gmp_randstate_t state);

/**
  Sample from ``dgs_rround_mp_t`` using the convolution sampler.

  :param rop: return value
  :param self: discrete Gaussian rounder
  :param sigma: noise parameter
  :param c: center
  :param state: state

 */

void dgs_rround_mp_call_convolution(mpz_t rop, dgs_rround_mp_t *self, const mpfr_t sigma, const mpfr_t c,
                                    gmp_randstate_t state);

/**
   Free memory.

   :param self: discrete Gaussian rounder

 */

void dgs_rround_mp_clear(dgs_rround_mp_t *self);

#endif  // DGS_RROUND__H
