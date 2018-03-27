/******************************************************************************
*
*                        DGS - Discrete Gaussian Samplers
*
* Copyright (c) 2014, Martin Albrecht  <martinralbrecht+dgs@googlemail.com>
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
#include <stdlib.h>
#include <limits.h>
#include <math.h>

/** GENERAL SIGMA :: INIT **/

static inline void _dgs_rround_mp_init_f(mpfr_t f, const mpfr_t sigma) {
  mpfr_set(f, sigma, MPFR_RNDN);
  mpfr_sqr(f, f, MPFR_RNDN); // f = σ²
  mpfr_mul_ui(f, f, 2, MPFR_RNDN); // f = 2 σ²
  mpfr_ui_div(f, 1, f, MPFR_RNDN);  // f = 1/(2 σ²)
  mpfr_neg(f, f, MPFR_RNDN); // f = -1/(2 σ²)
}

static inline void _dgs_rround_mp_init_upper_bound(mpz_t upper_bound,
                                                       mpz_t upper_bound_minus_one,
                                                       mpz_t two_upper_bound_minus_one,
                                                       const mpfr_t sigma, size_t tailcut) {
  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(sigma));
  mpfr_mul_ui(tmp, sigma, tailcut, MPFR_RNDN); // tmp = σ·τ
  mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN); // tmp = σ·τ + 1
  mpfr_get_z(upper_bound, tmp, MPFR_RNDU); // upper_bound = ⌈σ·τ + 1⌉
  mpz_sub_ui(upper_bound_minus_one, upper_bound, 1); // upper_bound - 1 = ⌈σ·τ⌉
  mpz_mul_ui(two_upper_bound_minus_one, upper_bound, 2);
  mpz_sub_ui(two_upper_bound_minus_one, two_upper_bound_minus_one, 1); // 2·upper_bound - 1
  mpfr_clear(tmp);
}

dgs_rround_mp_t *dgs_rround_mp_init(size_t tau, dgs_rround_alg_t algorithm, mpfr_prec_t prec) {
  if (tau == 0)
    dgs_die("tau must be > 0");

  dgs_rround_mp_t *self = (dgs_rround_mp_t*)calloc(sizeof(dgs_rround_mp_t),1);
  if (!self) dgs_die("out of memory");

  mpz_init(self->x);
  mpfr_init2(self->y, prec);
  mpfr_init2(self->z, prec);

  mpz_init(self->c_z);
  mpfr_init2(self->c_r, prec);

  self->tau = tau;

  if (algorithm == DGS_RROUND_DEFAULT) {
    algorithm = DGS_RROUND_UNIFORM_ONLINE;
  }
  self->algorithm = algorithm;

  switch(algorithm) {

  case DGS_RROUND_UNIFORM_ONLINE: {
    self->call = dgs_rround_mp_call_uniform_online;
    mpfr_init2(self->f, prec);
    mpz_init(self->upper_bound);
    mpz_init(self->upper_bound_minus_one);
    mpz_init(self->two_upper_bound_minus_one);
   break;
  }

  default:
    free(self);
    dgs_die("unknown algorithm %d", algorithm);
  }
  return self;
}

/** GENERAL SIGMA :: CALL **/

void dgs_rround_mp_call_uniform_online(mpz_t rop, dgs_rround_mp_t *self, const mpfr_t sigma, const mpfr_t c, gmp_randstate_t state) {
  if (mpfr_cmp_ui(sigma,0)<= 0)
    dgs_die("sigma must be > 0");
  
  mpfr_get_z(self->c_z, c, MPFR_RNDN);
  mpfr_sub_z(self->c_r, c, self->c_z, MPFR_RNDN);
  
  _dgs_rround_mp_init_upper_bound(self->upper_bound,
                                        self->upper_bound_minus_one,
                                        self->two_upper_bound_minus_one,
                                        sigma, self->tau);
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


/** GENERAL SIGMA :: CLEAR **/

void dgs_rround_mp_clear(dgs_rround_mp_t *self) {
  mpz_clear(self->x);
  mpfr_clear(self->y);
  mpfr_clear(self->f);
  mpfr_clear(self->z);
  mpfr_clear(self->c_r);
  mpz_clear(self->c_z);
  
  if (self->upper_bound)
    mpz_clear(self->upper_bound);
  if (self->upper_bound_minus_one)
    mpz_clear(self->upper_bound_minus_one);
  if (self->two_upper_bound_minus_one)
    mpz_clear(self->two_upper_bound_minus_one);

  free(self);
}
