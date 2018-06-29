/******************************************************************************
*
*                        DGR - Discrete Gaussian Rounders
*
* Copyright (c) 2018, Michael Walter
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* The views and conclusions contained in the software and documentation are
* those of the authors and should not be interpreted as representing official
* policies, either expressed or implied, of the FreeBSD Project.
******************************************************************************/

#include "dgs.h"
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

static inline void _dgs_rround_mp_init_f(mpfr_t f, const mpfr_t sigma) {
  mpfr_set(f, sigma, MPFR_RNDN);
  mpfr_sqr(f, f, MPFR_RNDN);        // f = σ²
  mpfr_mul_ui(f, f, 2, MPFR_RNDN);  // f = 2 σ²
  mpfr_ui_div(f, 1, f, MPFR_RNDN);  // f = 1/(2 σ²)
  mpfr_neg(f, f, MPFR_RNDN);        // f = -1/(2 σ²)
}

static inline void _dgs_rround_mp_init_upper_bound(mpz_t upper_bound, mpz_t upper_bound_minus_one,
                                                   mpz_t two_upper_bound_minus_one, const mpfr_t sigma, size_t tailcut,
                                                   mpfr_t tmp) {
  mpfr_mul_ui(tmp, sigma, tailcut, MPFR_RNDN);        // tmp = σ·τ
  mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);                // tmp = σ·τ + 1
  mpfr_get_z(upper_bound, tmp, MPFR_RNDU);            // upper_bound = ⌈σ·τ + 1⌉
  mpz_sub_ui(upper_bound_minus_one, upper_bound, 1);  // upper_bound - 1 = ⌈σ·τ⌉
  mpz_mul_ui(two_upper_bound_minus_one, upper_bound, 2);
  mpz_sub_ui(two_upper_bound_minus_one, two_upper_bound_minus_one, 1);  // 2·upper_bound - 1
}

static inline long _dgs_rround_mp_unit_gauss(dgs_rround_mp_t *self, gmp_randstate_t state) {
  long x;
  int reject;
  do {
    reject = 0;
    x      = 0;
    while (dgs_bern_mp_call(self->B_half_exp, state)) ++x;
    if (x < 2) return x;

    for (int i = 0; i < x * (x - 1); ++i) {
      if (dgs_bern_mp_call(self->B_half_exp, state) == 0) {
        reject = 1;
        break;
      }
    }
  } while (reject);

  return x;
}

dgs_rround_mp_t *dgs_rround_mp_init(size_t tau, dgs_rround_alg_t algorithm, mpfr_prec_t prec) {
  if (tau == 0) dgs_die("tau must be > 0");

  dgs_rround_mp_t *self = (dgs_rround_mp_t *)calloc(sizeof(dgs_rround_mp_t), 1);
  if (!self) dgs_die("out of memory");

  mpz_init(self->x);
  mpfr_init2(self->y, prec);
  mpfr_init2(self->z, prec);

  mpz_init(self->c_z);
  mpfr_init2(self->c_r, prec);

  mpfr_init2(self->tmp, prec);

  self->tau = tau;

  mpz_init(self->sigma_z);

  mpfr_init2(self->f, prec);
  mpz_init(self->upper_bound);
  mpz_init(self->upper_bound_minus_one);
  mpz_init(self->two_upper_bound_minus_one);

  if (algorithm == DGS_RROUND_DEFAULT) { algorithm = DGS_RROUND_KARNEY; }
  self->algorithm = algorithm;

  switch (algorithm) {

  case DGS_RROUND_UNIFORM_ONLINE: {
    self->call = dgs_rround_mp_call_uniform_online;
    break;
  }
  case DGS_RROUND_KARNEY: {
    self->call = dgs_rround_mp_call_karney;

    self->B = dgs_bern_uniform_init(0);

    char half_exp_str[] = "0.60653065971263342360379953499118045344191813548718695568289";
    mpfr_t half_exp;
    mpfr_init2(half_exp, prec);
    mpfr_set_str(half_exp, half_exp_str, 10, MPFR_RNDN);

    self->B_half_exp = dgs_bern_mp_init(half_exp);

    break;
  }
  default: free(self); dgs_die("unknown algorithm %d", algorithm);
  }
  return self;
}

void dgs_rround_mp_call_uniform_online(mpz_t rop, dgs_rround_mp_t *self, const mpfr_t sigma, const mpfr_t c,
                                       gmp_randstate_t state) {
  if (mpfr_cmp_ui(sigma, 0) <= 0) dgs_die("sigma must be > 0");

  mpfr_get_z(self->c_z, c, MPFR_RNDN);
  mpfr_sub_z(self->c_r, c, self->c_z, MPFR_RNDN);

  _dgs_rround_mp_init_upper_bound(self->upper_bound, self->upper_bound_minus_one, self->two_upper_bound_minus_one,
                                  sigma, self->tau, self->tmp);
  _dgs_rround_mp_init_f(self->f, sigma);

  do {
    mpz_urandomm(self->x, state, self->two_upper_bound_minus_one);
    mpz_sub(self->x, self->x, self->upper_bound_minus_one);
    mpfr_set_z(self->z, self->x, MPFR_RNDN);
    mpfr_sub(self->z, self->z, self->c_r, MPFR_RNDN);
    mpfr_mul(self->z, self->z, self->z, MPFR_RNDN);
    mpfr_mul(self->z, self->z, self->f, MPFR_RNDN);
    mpfr_exp(self->z, self->z, MPFR_RNDN);
    mpfr_urandomb(self->y, state);
  } while (mpfr_cmp(self->y, self->z) >= 0);

  mpz_set(rop, self->x);
  mpz_add(rop, rop, self->c_z);
}

void dgs_rround_mp_call_karney(mpz_t rop, dgs_rround_mp_t *self, const mpfr_t sigma, const mpfr_t c,
                               gmp_randstate_t state) {
  mpfr_get_z(self->sigma_z, sigma, MPFR_RNDU);
  do {
    long k = _dgs_rround_mp_unit_gauss(self, state);
    long s = 1;
    if (dgs_bern_uniform_call_libc(self->B)) s *= -1;

    // double tmp = k*self->sigma + s*self->c; tmp = self->y
    mpfr_mul_si(self->y, sigma, k, MPFR_RNDN);
    mpfr_mul_si(self->z, c, s, MPFR_RNDN);
    mpfr_add(self->y, self->y, self->z, MPFR_RNDN);

    // long i0 = (long)ceil(tmp); i0 = rop
    mpfr_get_z(rop, self->y, MPFR_RNDU);

    // double x0 = (i0 - tmp)/self->sigma; x0 = self->y
    mpfr_z_sub(self->y, rop, self->y, MPFR_RNDN);
    mpfr_div(self->y, self->y, sigma, MPFR_RNDN);

    // long j = _dgs_randomm_libc((unsigned long)ceil(self->sigma)); j = self->x
    mpz_urandomm(self->x, state, self->sigma_z);

    // double x = x0 + ((double)j)/self->sigma; x = self->z
    mpfr_si_div(self->z, 1, sigma, MPFR_RNDN);
    mpfr_mul_z(self->z, self->z, self->x, MPFR_RNDN);
    mpfr_add(self->z, self->z, self->y, MPFR_RNDN);

    //~ if (x >= 1) {
    if (mpfr_cmp_si(self->z, 1) >= 0) { continue; }

    //~ if (x == 0) {
    if (mpfr_zero_p(self->z)) {
      if (k == 0 && s < 0) {
        continue;
      } else {
        mpz_add(rop, rop, self->x);
        mpz_mul_si(rop, rop, s);
        break;
      }
    }

    //~ double bias = exp(-.5*x*(2*k+x)); bias = self->y
    mpfr_add_si(self->y, self->z, 2 * k, MPFR_RNDN);
    mpfr_mul(self->y, self->z, self->y, MPFR_RNDN);
    mpfr_mul_d(self->y, self->y, -.5, MPFR_RNDN);
    mpfr_exp(self->y, self->y, MPFR_RNDN);

    mpfr_urandomb(self->z, state);

    //~ if (drand48() <= bias) {
    if (mpfr_cmp(self->z, self->y) <= 0) {
      mpz_add(rop, rop, self->x);
      mpz_mul_si(rop, rop, s);
      break;
    }

  } while (1);
}

void dgs_rround_mp_clear(dgs_rround_mp_t *self) {
  mpz_clear(self->x);
  mpfr_clear(self->y);
  mpfr_clear(self->z);
  mpfr_clear(self->c_r);
  mpz_clear(self->c_z);
  mpfr_clear(self->tmp);

  if (self->sigma_z) mpz_clear(self->sigma_z);
  if (self->f) mpfr_clear(self->f);
  if (self->upper_bound) mpz_clear(self->upper_bound);
  if (self->upper_bound_minus_one) mpz_clear(self->upper_bound_minus_one);
  if (self->two_upper_bound_minus_one) mpz_clear(self->two_upper_bound_minus_one);

  if (self->B) dgs_bern_uniform_clear(self->B);
  if (self->B_half_exp) dgs_bern_mp_clear(self->B_half_exp);

  free(self);
}
