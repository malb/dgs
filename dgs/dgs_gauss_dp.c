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
#include <math.h>
#include <stdlib.h>

static inline void _dgs_disc_gauss_dp_init_bexp(dgs_disc_gauss_dp_t *self, double sigma, long upper_bound) {
  self->f    = (2 * sigma * sigma);
  size_t l   = 2 * ceil(log2(upper_bound));
  self->Bexp = dgs_bern_exp_dp_init(self->f, l);
}

static inline void _dgs_disc_gauss_dp_init_rho(dgs_disc_gauss_dp_t *self) {
  self->rho = (double *)malloc(sizeof(double) * self->two_upper_bound_minus_one);
  if (!self->rho) {
    dgs_disc_gauss_dp_clear(self);
    dgs_die("out of memory");
  }
  long absmax = self->upper_bound_minus_one;
  for (long x = -absmax; x <= absmax; x++) {
    self->rho[x + self->upper_bound_minus_one] = exp((((double)x) - self->c_r) * (((double)x) - self->c_r) * self->f);
  }
}

static inline long _dgs_disc_gauss_dp_min_in_rho(dgs_disc_gauss_dp_t *self) {
  long mi  = 0;
  double m = self->rho[mi];
  for (long x = 1; x < self->two_upper_bound_minus_one; ++x) {
    if (self->rho[x] < m) {
      mi = x;
      m  = self->rho[mi];
    }
  }
  return mi;
}

static inline long _dgs_disc_gauss_dp_max_in_rho(dgs_disc_gauss_dp_t *self) {
  long mi  = 0;
  double m = self->rho[mi];
  for (long x = 1; x < self->two_upper_bound_minus_one; ++x) {
    if (self->rho[x] > m) {
      mi = x;
      m  = self->rho[mi];
    }
  }
  return mi;
}

dgs_disc_gauss_dp_t *dgs_disc_gauss_dp_init(double sigma, double c, size_t tau, dgs_disc_gauss_alg_t algorithm) {
  if (sigma <= 0.0) dgs_die("sigma must be > 0");
  if (tau == 0) dgs_die("tau must be > 0");

  size_t upper_bound;

  dgs_disc_gauss_dp_t *self = (dgs_disc_gauss_dp_t *)calloc(sizeof(dgs_disc_gauss_dp_t), 1);
  if (!self) dgs_die("out of memory");

  self->sigma = sigma;
  self->c     = c;
  self->c_z   = (long)c;
  self->c_r   = self->c - ((double)self->c_z);
  self->tau   = tau;

  double sigma2 = sqrt(1.0 / (2 * log(2.0)));
  double k      = sigma / sigma2;

  if (algorithm == DGS_DISC_GAUSS_DEFAULT) {
    /* 1. try the uniform algorithm */
    if (2 * ceil(self->sigma * tau) * sizeof(double) <= DGS_DISC_GAUSS_MAX_TABLE_SIZE_BYTES) {
      algorithm = DGS_DISC_GAUSS_UNIFORM_TABLE;
      /* 2. see if sigma2 is close enough */
    } else if (fabs(round(k) - k) < DGS_DISC_GAUSS_EQUAL_DIFF) {
      algorithm = DGS_DISC_GAUSS_SIGMA2_LOGTABLE;
      /* 3. do logtables */
    } else {
      algorithm = DGS_DISC_GAUSS_UNIFORM_LOGTABLE;
    }
  }
  self->algorithm = algorithm;

  switch (algorithm) {

  case DGS_DISC_GAUSS_UNIFORM_ONLINE:
    self->call = dgs_disc_gauss_dp_call_uniform_online;

    upper_bound                     = ceil(self->sigma * tau) + 1;
    self->upper_bound               = upper_bound;
    self->upper_bound_minus_one     = upper_bound - 1;
    self->two_upper_bound_minus_one = 2 * upper_bound - 1;
    self->f                         = -1.0 / (2.0 * (self->sigma * self->sigma));
    break;

  case DGS_DISC_GAUSS_UNIFORM_TABLE:
    self->call = dgs_disc_gauss_dp_call_uniform_table;

    upper_bound                     = ceil(self->sigma * tau) + 1;
    self->upper_bound               = upper_bound;
    self->upper_bound_minus_one     = upper_bound - 1;
    self->two_upper_bound_minus_one = 2 * upper_bound - 1;
    self->B                         = dgs_bern_uniform_init(0);
    self->f                         = -1.0 / (2.0 * (sigma * sigma));

    if (self->c_r == 0) {
      self->call = dgs_disc_gauss_dp_call_uniform_table;
      self->rho  = (double *)malloc(sizeof(double) * self->upper_bound);
      if (!self->rho) {
        dgs_disc_gauss_dp_clear(self);
        dgs_die("out of memory");
      }
      for (unsigned long x = 0; x < self->upper_bound; x++) {
        self->rho[x] = exp((((double)x) - self->c_r) * (((double)x) - self->c_r) * self->f);
      }
      self->rho[0] /= 2.0;
    } else {
      self->call = dgs_disc_gauss_dp_call_uniform_table_offset;
      _dgs_disc_gauss_dp_init_rho(self);
    }
    break;

  case DGS_DISC_GAUSS_UNIFORM_LOGTABLE:
    self->call = dgs_disc_gauss_dp_call_uniform_logtable;

    if (fabs(self->c_r) > DGS_DISC_GAUSS_INTEGER_CUTOFF) {
      dgs_disc_gauss_dp_clear(self);
      dgs_die("algorithm DGS_DISC_GAUSS_UNIFORM_LOGTABLE requires c%1 == 0");
    }
    upper_bound                     = ceil(self->sigma * tau) + 1;
    self->upper_bound               = upper_bound;
    self->upper_bound_minus_one     = upper_bound - 1;
    self->two_upper_bound_minus_one = 2 * upper_bound - 1;

    _dgs_disc_gauss_dp_init_bexp(self, self->sigma, self->upper_bound);
    break;

  case DGS_DISC_GAUSS_SIGMA2_LOGTABLE: {
    self->call = dgs_disc_gauss_dp_call_sigma2_logtable;

    if (fabs(self->c_r) > DGS_DISC_GAUSS_INTEGER_CUTOFF) {
      dgs_disc_gauss_dp_clear(self);
      dgs_die("algorithm DGS_DISC_GAUSS_SIGMA2_LOGTABLE requires c%1 == 0");
    }

    self->k     = round(k);
    self->sigma = self->k * sigma2;

    upper_bound                     = ceil(self->sigma * tau) + 1;
    self->upper_bound               = upper_bound;
    self->upper_bound_minus_one     = upper_bound - 1;
    self->two_upper_bound_minus_one = 2 * upper_bound - 1;

    _dgs_disc_gauss_dp_init_bexp(self, self->sigma, self->upper_bound);
    self->B  = dgs_bern_uniform_init(0);
    self->D2 = dgs_disc_gauss_sigma2p_init();
    break;
  }

  case DGS_DISC_GAUSS_ALIAS: {
    self->call = dgs_disc_gauss_dp_call_alias;

    upper_bound                     = ceil(self->sigma * tau) + 1;
    self->upper_bound               = upper_bound;
    self->upper_bound_minus_one     = upper_bound - 1;
    self->two_upper_bound_minus_one = 2 * upper_bound - 1;
    self->B                         = dgs_bern_uniform_init(0);
    self->f                         = -1.0 / (2.0 * (sigma * sigma));

    _dgs_disc_gauss_dp_init_rho(self);

    // convert rho to probabilities
    double sum = 0;
    for (long x = 0; x < self->two_upper_bound_minus_one; x++) { sum += self->rho[x]; }
    sum         = 1 / sum;
    for (long x = 0; x < self->two_upper_bound_minus_one; x++) { self->rho[x] *= sum; }

    // compute bias and alias
    self->alias = (long*)malloc(sizeof(long)*self->two_upper_bound_minus_one);
    self->bias = (dgs_bern_dp_t**)malloc(sizeof(dgs_bern_dp_t*)*self->two_upper_bound_minus_one);
    
    if (!self->alias || !self->bias){
      dgs_disc_gauss_dp_clear(self);
      dgs_die("out of memory");
    }
    
    // simple robin hood strategy approximates good alias
    // this precomputation takes ~n^2, but could be reduced by
    // using better data structures to compute min and max
    // (instead of just linear search each time)
    double avg = 1.0 / ((double)self->two_upper_bound_minus_one);
    long low   = _dgs_disc_gauss_dp_min_in_rho(self);
    long high;
    while (avg - self->rho[low] > DGS_DISC_GAUSS_STRONG_EQUAL_DIFF) {
      high = _dgs_disc_gauss_dp_max_in_rho(self);

      self->bias[low]  = dgs_bern_dp_init(self->two_upper_bound_minus_one * self->rho[low]);
      self->alias[low] = high;
      self->rho[high] -= (avg - self->rho[low]);
      self->rho[low] = avg;

      low = _dgs_disc_gauss_dp_min_in_rho(self);
    }

    break;
  }

  case DGS_DISC_GAUSS_CONVOLUTION: {
    self->call = dgs_disc_gauss_dp_call_convolution;
    
    double eta = 2;
    double sigma1 = sigma;
    double coset_sigma = sqrt(2)*eta;
    if (fabs(self->c_r) > 0) {
      // we might need to adjust the center
      sigma1 = sqrt(sigma*sigma - coset_sigma*coset_sigma);
    }
    
    long table_size = 2*ceil(sigma1*tau) * (sizeof(dgs_bern_dp_t) + sizeof(long));
    int recursion_level = 0;
    double current_sigma = sigma1;
    long z1 = 1;
    long z2 = 1;
    
    // compute recursion level for convolution
    while (table_size > DGS_DISC_GAUSS_MAX_TABLE_SIZE_BYTES) {
      recursion_level++;
      z1 = floor(sqrt(current_sigma/(eta*2)));
      if (z1 == 0) {
        dgs_disc_gauss_dp_clear(self);
        dgs_die("MAX_TABLE_SIZE too small to store alias sampler!");
      }
      z2 = (z1 > 1) ? z1 - 1 : 1;
      current_sigma /= (sqrt(z1*z1 + z2*z2));
      table_size = 2*ceil(current_sigma*tau) * (sizeof(dgs_bern_dp_t) + sizeof(long));
    }
    
    self->n_coefficients = 1 << recursion_level;
    self->coefficients = (long*)malloc(sizeof(long)*self->n_coefficients);
    for (int i = 0; i < self->n_coefficients; ++i) {
      self->coefficients[i] = 1;
    }
    
    // if there is no convolution, we simply forward to alias and
    // so we won't need adjustment of sigma, even if sampling off center
    current_sigma = (recursion_level == 0)? sigma : sigma1;
    
    // redo above computation to store coefficients
    for (int i = 0; i < recursion_level; ++i) {
      z1 = floor(sqrt(current_sigma/(eta*2)));
      z2 = (z1 > 1) ? z1 - 1 : 1;
      
      // we unroll the recursion on the coefficients on the fly
      // so we don't have to use recursion during the call
      int off = (1 << recursion_level - i - 1);
      for (int j = 0; j < (1 << i); ++j) {
        for (int k = 0; k < off;++k) {
          self->coefficients[2*j*off + k] *= z1;
        }
      }
      
      for (int j = 0; j < (1 << i); ++j) {
        for (int k = 0; k < off;++k) {
          self->coefficients[(2*j + 1)*off + k] *= z2;
        }
      }
      
      current_sigma /= (sqrt(z1*z1 + z2*z2));
    }
    
    double base_c = self->c_r;
    self->coset_sampler = NULL;
    if (recursion_level > 0 && fabs(self->c_r) > 0) {
      // we'll need to adjust the center
      base_c = 0.0;
      self->coset_sampler = dgs_disc_gauss_dp_init(coset_sigma, self->c_r, tau, DGS_DISC_GAUSS_ALIAS);
    }
    self->base_sampler = dgs_disc_gauss_dp_init(current_sigma, base_c, tau, DGS_DISC_GAUSS_ALIAS);
    
    break;
  }
  
  default:
    dgs_disc_gauss_dp_clear(self);
    dgs_die("unknown algorithm %d", algorithm);
  }
  return self;
}

long dgs_disc_gauss_dp_call_uniform_online(dgs_disc_gauss_dp_t *self) {
  long x;
  double y, z;
  double c = self->c;
  do {
    x = self->c_z + _dgs_randomm_libc(self->two_upper_bound_minus_one) - self->upper_bound_minus_one;
    z = exp(((double)x - c) * ((double)x - c) * self->f);
    y = drand48();
  } while (y >= z);

  return x;
}

long dgs_disc_gauss_dp_call_uniform_table(dgs_disc_gauss_dp_t *self) {
  long x;
  double y;
  do {
    x = _dgs_randomm_libc(self->upper_bound);
    y = drand48();
  } while (y >= self->rho[x]);

  if (dgs_bern_uniform_call_libc(self->B)) x = -x;
  return x + self->c_z;
}

long dgs_disc_gauss_dp_call_uniform_table_offset(dgs_disc_gauss_dp_t *self) {
  long x;
  double y;
  do {
    x = _dgs_randomm_libc(self->two_upper_bound_minus_one);
    y = drand48();
  } while (y >= self->rho[x]);

  return x + self->c_z - self->upper_bound_minus_one;
}

long dgs_disc_gauss_dp_call_alias(dgs_disc_gauss_dp_t *self) {
  long x = _dgs_randomm_libc(self->two_upper_bound_minus_one);
  if (self->bias[x]) {
    if (!dgs_bern_dp_call(self->bias[x])) { x = self->alias[x]; }
  }

  return x + self->c_z - self->upper_bound_minus_one;
}

long dgs_disc_gauss_dp_call_convolution(dgs_disc_gauss_dp_t *self) {
  long x = 0;
  for (int i = 0; i < self->n_coefficients; ++i) {
    x += self->coefficients[i]*self->base_sampler->call(self->base_sampler);
  }
  
  // adjust center if necessary
  if (self->coset_sampler) {
    x += self->coset_sampler->call(self->coset_sampler);
  }
  return x + self->c_z;
}

long dgs_disc_gauss_dp_call_uniform_logtable(dgs_disc_gauss_dp_t *self) {
  long x;
  do {
    x = _dgs_randomm_libc(self->two_upper_bound_minus_one) - self->upper_bound_minus_one;
  } while (dgs_bern_exp_dp_call(self->Bexp, x * x) == 0);
  return x + self->c_z;
}

long dgs_disc_gauss_dp_call_sigma2_logtable(dgs_disc_gauss_dp_t *self) {
  long x, y, z;
  long k = self->k;

  do {
    do {
      x = dgs_disc_gauss_sigma2p_dp_call(self->D2);
      y = _dgs_randomm_libc(self->k);
    } while (dgs_bern_exp_dp_call(self->Bexp, y * (y + 2 * k * x)) == 0);
    z = k * x + y;
    if (!z) {
      if (dgs_bern_uniform_call_libc(self->B)) break;
    } else {
      break;
    }
  } while (1);
  if (dgs_bern_uniform_call_libc(self->B)) z = -z;
  return z + self->c_z;
}

void dgs_disc_gauss_dp_clear(dgs_disc_gauss_dp_t *self) {
  assert(self != NULL);
  if (self->B) dgs_bern_uniform_clear(self->B);
  if (self->Bexp) dgs_bern_exp_dp_clear(self->Bexp);
  if (self->rho) free(self->rho);
  if (self->alias) free(self->alias);
  if (self->bias) {
    for (long x = 0; x < self->two_upper_bound_minus_one; x++) {
      if (self->bias[x]) { free(self->bias[x]); }
    }
    free(self->bias);
  }
  if (self->base_sampler) {
    dgs_disc_gauss_dp_clear(self->base_sampler);
  }
  if (self->coefficients) {
    free(self->coefficients);
  }
  if (self->coset_sampler) {
    dgs_disc_gauss_dp_clear(self->coset_sampler);
  }
  
  free(self);
}
