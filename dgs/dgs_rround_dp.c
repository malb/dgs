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
#include <math.h>
#include <stdlib.h>

static inline long _dgs_rround_dp_unit_gauss(dgs_rround_dp_t *self) {
  long x;
  int reject;
  do {
    reject = 0;
    x      = 0;
    while (dgs_bern_dp_call(self->B_half_exp)) ++x;
    if (x < 2) return x;

    for (int i = 0; i < x * (x - 1); ++i) {
      if (dgs_bern_dp_call(self->B_half_exp) == 0) {
        reject = 1;
        break;
      }
    }
  } while (reject);

  return x;
}

dgs_rround_dp_t *dgs_rround_dp_init(size_t tau, dgs_rround_alg_t algorithm) {
  if (tau == 0) dgs_die("tau must be > 0");

  dgs_rround_dp_t *self = (dgs_rround_dp_t *)calloc(sizeof(dgs_rround_dp_t), 1);
  if (!self) dgs_die("out of memory");

  self->tau = tau;

  if (algorithm == DGS_RROUND_DEFAULT) { algorithm = DGS_RROUND_KARNEY; }
  self->algorithm = algorithm;

  switch (algorithm) {

  case DGS_RROUND_UNIFORM_ONLINE: self->call = dgs_rround_dp_call_uniform_online; break;
  case DGS_RROUND_KARNEY: {
    self->call = dgs_rround_dp_call_karney;

    self->B          = dgs_bern_uniform_init(0);
    self->B_half_exp = dgs_bern_dp_init(exp(-.5));

    break;
  }
  default: dgs_rround_dp_clear(self); dgs_die("unknown algorithm %d", algorithm);
  }
  return self;
}

long dgs_rround_dp_call_uniform_online(dgs_rround_dp_t *self, double sigma, double c) {
  if (sigma <= 0.0) dgs_die("sigma must be > 0");

  size_t upper_bound               = ceil(sigma * self->tau) + 1;
  size_t upper_bound_minus_one     = upper_bound - 1;
  size_t two_upper_bound_minus_one = 2 * upper_bound - 1;
  double f                         = -1.0 / (2.0 * (sigma * sigma));

  long x;
  double y, z;
  do {
    x = ((long)c) + _dgs_randomm_libc(two_upper_bound_minus_one) - upper_bound_minus_one;
    z = exp(((double)x - c) * ((double)x - c) * f);
    y = drand48();
  } while (y >= z);

  return x;
}

long dgs_rround_dp_call_karney(dgs_rround_dp_t *self, double sigma, double c) {
  do {
    long k = _dgs_rround_dp_unit_gauss(self);

    long s = 1;
    if (dgs_bern_uniform_call_libc(self->B)) s *= -1;

    double tmp = k * sigma + s * c;
    long i0    = (long)ceil(tmp);
    double x0  = (i0 - tmp) / sigma;
    long j     = _dgs_randomm_libc((unsigned long)ceil(sigma));
    double x   = x0 + ((double)j) / sigma;

    if (x >= 1) { continue; }

    if (x == 0) {
      if (k == 0 && s < 0) {
        continue;
      } else {
        return s * (i0 + j);
      }
    }

    double bias = exp(-.5 * x * (2 * k + x));
    if (drand48() <= bias) { return s * (i0 + j); }
  } while (1);
}

void dgs_rround_dp_clear(dgs_rround_dp_t *self) {
  assert(self != NULL);
  if (self->B) dgs_bern_uniform_clear(self->B);
  if (self->B_half_exp) dgs_bern_dp_clear(self->B_half_exp);

  free(self);
}
