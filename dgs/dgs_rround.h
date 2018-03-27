/**
   Generic Discrete Gaussian Rounding over the Integers.

   A discrete Gaussian distribution on the Integers is a distribution where the
   integer `x` is sampled with probability proportional to `exp(-(x-c)²/(2σ²))`.
   It is denoted by `D_{σ,c}` where `σ` is the width parameter (close to the
   standard deviation) and `c` is the center.

   AVAILABLE ALGORITHMS:

 - ``DGS_DISC_GAUSS_UNIFORM_ONLINE`` - samples are drawn from a uniform
   distribution and accepted with probability proportional to
   `\exp(-(x-c)²/(2σ²))` where `\exp(-(x-c)²/(2σ²))` is computed in each
   invocation. Typically this is very slow. Any real-valued `c` is accepted.

  AVAILABLE PRECISIONS:

  - ``mp`` - multi-precision using MPFR, cf. ``dgs_gauss_mp.c``

  - ``dp`` - double precision using machine doubles, cf. ``dgs_gauss_dp.c``.

  For readers unfamiliar with the implemented algorithms it makes sense to start
  with ``dgs_gauss_dp.c`` which implements the same algorithms as
  ``dgs_gauss_mp.c`` should be easier to read.

  TYPICAL USAGE::

      dgs_disc_gauss_dp_t *D = dgs_disc_gauss_dp_init(<sigma>, <c>, <tau>, <algorithm>);
      D->call(D); // as often as needed
      dgs_disc_gauss_dp_clear(D);

   .. author:: Michael Walter

*/

#ifndef DGS_RROUND__H
#define DGS_RROUND__H

#include "dgs_gauss.h"

/**
   Available Algorithms
*/

typedef enum {
  DGS_RROUND_DEFAULT           = 0x0, //<pick algorithm
  DGS_RROUND_UNIFORM_ONLINE    = 0x1, //<call dgs_disc_gauss_mp_call_uniform_online
} dgs_rround_alg_t;

struct _dgs_rround_mp_t;

typedef struct _dgs_rround_dp_t {

  /**
     Cutoff `τ`, samples outside the range `(⌊c⌉-⌈στ⌉,...,⌊c⌉+⌈στ⌉)` are
     considered to have probability zero. This bound applies to algorithms
     which sample from the uniform distribution.
  */

  size_t tau;

  dgs_rround_alg_t algorithm;  //<  which algorithm to use

  /**
   Return an ``long`` sampled from this sampler

   :param self: discrete Gaussian sampler.
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

 */

long dgs_rround_dp_call_uniform_online(dgs_rround_dp_t *self, double sigma, double c);

/**
   Free memory.

   :param self: discrete Gaussian rounder

 */

void dgs_rround_dp_clear(dgs_rround_dp_t *self);


//~ /**
   //~ Multi-precision Discrete Gaussians `D_{σ,c}`

   //~ Return integer `x` with probability

   //~ `ρ_{σ,c}(x) = exp(-(x-c)²/(2σ²))/exp(-(\ZZ-c)²/(2σ²))`

   //~ where `exp(-(\ZZ-c)²/(2σ²)) ≈ \sum_{i=-τσ}^{τσ} exp(-(i-c)²/(2σ^²))` is the
   //~ probability for all of the integers.

//~ */

//~ typedef struct _dgs_disc_gauss_mp_t {

  //~ /**
      //~ The width paramter `σ`, i.e. samples are accepted with probability
      //~ proportional to `\exp(-(x-c)²/(2σ²))`
   //~ */

  //~ mpfr_t sigma;

  //~ /**
     //~ The mean of the distribution `c`. The value of `c` does not have to be an
     //~ integer. However, some algorithms only support integer-valued `c`.
  //~ */

  //~ mpfr_t c;

  //~ mpfr_t c_r; //< `c_r := c % 1`
  //~ mpz_t c_z;  //< c_z := c - (c_r)

  //~ /**
     //~ Cutoff `τ`, samples outside the range `(⌊c⌉-⌈στ⌉,...,⌊c⌉+⌈στ⌉)` are
     //~ considered to have probability zero. This bound applies to algorithms
     //~ which sample from the uniform distribution.
  //~ */

  //~ size_t tau;

  //~ dgs_disc_gauss_alg_t algorithm; //<  which algorithm to use

  //~ /**
     //~ We use a uniform Bernoulli to decide signs.
   //~ */

  //~ dgs_bern_uniform_t *B;

  //~ /**
     //~ To realise rejection sampling, we call `B_{exp(-(x·x)/(2σ²))}` and accept
     //~ if it returns 1.

     //~ Used when ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE`` or
     //~ ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE`` is set.
   //~ */

  //~ dgs_bern_exp_mp_t *Bexp;


  //~ /**
     //~ `D_{σ₂,0}` which is easily sampable`

     //~ Used when ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE`` is set.
  //~ */

  //~ dgs_disc_gauss_sigma2p_t *D2;

  //~ /**
   //~ Return an ``mpz_t`` sampled from this sampler

   //~ :param rop: target value.
   //~ :param self: discrete Gaussian sampler.
   //~ :param state: entropy pool.

  //~ */

  //~ void (*call)(mpz_t rop, struct _dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

  //~ /**
   //~ We sample ``x`` with ``abs(x) < upper_bound`` in
   //~ ``DGS_DISC_GAUSS_UNIFORM_ONLINE``, ``DGS_DISC_GAUSS_UNIFORM_TABLE`` and
   //~ ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE``.
   //~ */

  //~ mpz_t upper_bound;

  //~ /**
   //~ We sample ``x`` with ``abs(x) <= upper_bound - 1`` in
   //~ ``DGS_DISC_GAUSS_UNIFORM_ONLINE``, ``DGS_DISC_GAUSS_UNIFORM_TABLE`` and
   //~ ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE``.
   //~ */

  //~ mpz_t upper_bound_minus_one;

  //~ /**
     //~ There are ``2*upper_bound -1`` elements in the range
     //~ ``-upper_bound+1,...,upper_bound-1``.
   //~ */

  //~ mpz_t two_upper_bound_minus_one;

  //~ /**
     //~ The multiplier `k` when we sample from `D_{k·σ₂,c}` in
     //~ ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE``.
  //~ */

  //~ mpz_t k;

  //~ /**
   //~ Precomputed `-1/(2σ²)`.
  //~ */

  //~ mpfr_t f;

  //~ mpz_t x; //< space for temporary integer
  //~ mpz_t y_z; //< space for temporary integer
  //~ mpz_t x2; // space for temporary integer
  //~ mpfr_t y; // space for temporary rational number
  //~ mpfr_t z; // space for temporary rational number

  //~ /**
     //~ Precomputed values for `exp(-(x-c)²/(2σ²))` in
     //~ ``DGS_DISC_GAUSS_UNIFORM_TABLE``
  //~ */

  //~ mpfr_t *rho;
  
    //~ /**
   //~ * Tables required for alias sampling.
   //~ */
   
  //~ mpz_t* alias;
  //~ dgs_bern_mp_t** bias;

//~ } dgs_disc_gauss_mp_t;

//~ dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(const mpfr_t sigma, const mpfr_t c, size_t tau, dgs_disc_gauss_alg_t algorithm);

//~ /**
   //~ Sample from ``dgs_disc_gauss_mp_t`` by rejection sampling using the uniform
   //~ distribution and tabulated ``exp()`` evaluations.

   //~ :param self: discrete Gaussian sampler

   //~ .. note::

      //~ `c` must be an integer in this algorithm

 //~ */

//~ void dgs_disc_gauss_mp_call_uniform_table(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

//~ /**
   //~ Sample from ``dgs_disc_gauss_mp_t`` by rejection sampling using the uniform
   //~ distribution and tabulated ``exp()`` evaluations.

   //~ :param self: discrete Gaussian sampler

   //~ .. note::

      //~ This function makes no assumptions about `c` but requires more resources
      //~ than ``dgs_disc_gauss_dp_call_uniform_table()``.

 //~ */

//~ void dgs_disc_gauss_mp_call_uniform_table_offset(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

//~ /**
   //~ Sample from ``dgs_disc_gauss_mp_t`` by alias sampling. This is extremely fast, 
   //~ but requires more resources and setup cost is around (2τσ)².

   //~ :param self: discrete Gaussian sampler
 //~ */
//~ void dgs_disc_gauss_mp_call_alias(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

//~ /**
  //~ Sample from ``dgs_disc_gauss_mp_t`` by rejection sampling using the uniform
  //~ distribution replacing all ``exp()`` calls with call to Bernoulli distributions.

  //~ :param self: discrete Gaussian sampler

  //~ .. note::

     //~ `c` must be an integer in this algorithm
 //~ */

//~ void dgs_disc_gauss_mp_call_uniform_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

//~ /**
  //~ Sample from ``dgs_disc_gauss_mp_t`` by rejection sampling using the uniform distribution.

  //~ :param self: discrete Gaussian sampler

 //~ */

//~ void dgs_disc_gauss_mp_call_uniform_online(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

//~ /**
  //~ Sample from ``dgs_disc_gauss_mp_t`` by rejection sampling using the `D_{k·σ₂,0}`
  //~ distribution replacing all ``exp()`` calls with call to Bernoulli distributions.

  //~ :param self: Discrete Gaussian sampler

  //~ .. note::

     //~ `c` must be an integer in this algorithm.
 //~ */

//~ void dgs_disc_gauss_mp_call_sigma2_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self, gmp_randstate_t state);

//~ /**
   //~ Clear cache of random bits.

   //~ :param self: discrete Gaussian sampler

 //~ */

//~ static inline void dgs_disc_gauss_mp_flush_cache(dgs_disc_gauss_mp_t *self) {
  //~ self->B->count = self->B->length;
//~ }

//~ /**
   //~ Free memory.

   //~ :param self: discrete Gaussian sadpler

 //~ */

//~ void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self);

#endif //DGS_RROUND__H
