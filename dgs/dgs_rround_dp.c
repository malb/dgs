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

double _box_muller(dgs_rround_dp_t *self) {
  if (self->pool) {
    self->pool--;
    return self->bm_sample;
  }

  double r1 = drand48();
  double r2 = drand48();

  double R        = sqrt(-2 * log(r1));
  self->bm_sample = R * cos(2 * M_PI * r2);
  self->pool      = 1;
  return R * sin(2 * M_PI * r2);
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
  case DGS_RROUND_CONVOLUTION: {
    self->call = dgs_rround_dp_call_convolution;

    self->pool = 0;

    // we set parameters so that the memory does not exceed
    // DGS_DISC_GAUSS_MAX_TABLE_SIZE_BYTES
    double eta        = 2;    // smoothing paramter
    double base_sigma = 2.5;  // sigma for 2^b base samplers

    long base_sampler_size = 2 * ceil(base_sigma * tau) * (sizeof(dgs_bern_dp_t) + sizeof(long));
    int base               = (DGS_DISC_GAUSS_MAX_TABLE_SIZE_BYTES) / base_sampler_size;
    self->log_base         = 0;
    while (base >>= 1) { ++self->log_base; }  // we want a power of 2
    base       = 1 << self->log_base;
    self->mask = __DGS_LSB_BITMASK(self->log_base);
    // we can now actually reduce base_sigma a little to save some
    // more memory and increase the range of this function w.r.t. sigma
    base_sigma          = sqrt(((double)(base + 1)) / base) * eta;
    self->base_samplers = (dgs_disc_gauss_dp_t **)malloc(sizeof(dgs_disc_gauss_dp_t *) * base);
    if (!self->base_samplers) {
      dgs_rround_dp_clear(self);
      dgs_die("out of memory");
    }

    for (int i = 0; i < base; ++i) {
      self->base_samplers[i] = dgs_disc_gauss_dp_init(base_sigma, ((double)i) / base, tau, DGS_DISC_GAUSS_ALIAS);
    }

    // we assume for the precision we're targeting here using
    // gaussian rounding on the first 25 bits is sufficient
    // do the rest using bernoulli (linear interpolation of gaussian)
    self->digits = (int)ceil(25.0 / self->log_base);

    // compute rr_sigma2
    self->s_bar2  = 1;
    long double t = 1.0 / (base * base);
    long double s = 1;
    for (int i = 1; i < self->digits; ++i) {
      s *= t;
      self->s_bar2 += s;
    }
    self->s_bar2 *= (base_sigma * base_sigma);

    // we use karney as fallback for small sigma, so we initialize it
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

long dgs_rround_dp_call_convolution(dgs_rround_dp_t *self, double sigma, double c) {
  double sigma2 = sigma * sigma;

  // we need sigma to be larger than s_bar
  // otherwise fall back to karney
  if (self->s_bar2 > sigma2) { return dgs_rround_dp_call_karney(self, sigma, c); }
  double K = sqrt(sigma2 - self->s_bar2);
  // we use continuous gaussians instead of wide samplers
  // this is faster and (provably) works the same way
  double xr = _box_muller(self);
  double c1 = c + K * xr;

  long c1_z = (long)floor(c1);
  c1 -= floor(c1);  // 0 <= c1 < 1

  c1 *= (1UL << self->digits * self->log_base);
  int64_t center = (int64_t)c1;
  c1 -= center;  // 0 <= c1 < 1

  if (drand48() < c1) ++center;

  for (int i = 0; i < self->digits; ++i) {
    long x = self->base_samplers[center & self->mask]->call(self->base_samplers[center & self->mask]);
    if ((self->mask & center) > 0 && center < 0) x -= 1;

    for (int j = 0; j < self->log_base; ++j) { center /= 2; }
    center += x;
  }

  return center + c1_z;
}

void dgs_rround_dp_clear(dgs_rround_dp_t *self) {
  assert(self != NULL);
  if (self->B) dgs_bern_uniform_clear(self->B);
  if (self->B_half_exp) dgs_bern_dp_clear(self->B_half_exp);
  if (self->base_samplers) {
    for (int i = 0; i < (1 << self->log_base); ++i) {
      if (self->base_samplers[i]) { dgs_disc_gauss_dp_clear(self->base_samplers[i]); }
    }
    free(self->base_samplers);
  }

  free(self);
}
