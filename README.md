# Discrete Gaussians over the Integers #

A discrete Gaussian distribution on the Integers is a distribution where the
integer $x$ is sampled with probability proportional to $exp(-(x-c)²/(2σ²))$.
It is denoted by $D_{σ,c}$ where `σ` is the width parameter (close to the
standard deviation) and $c$ is the center.

This library samples from this distributions.

## Algorithms ##

- `DGS_DISC_GAUSS_UNIFORM_TABLE` - classical rejection sampling, sampling from
  the uniform distribution and accepted with probability proportional to
  $\exp(-(x-c)²/(2σ²))$ where $\exp(-(x-c)²/(2σ²))$ is precomputed and stored in
  a table. Any real-valued `c` is supported.

- `DGS_DISC_GAUSS_UNIFORM_LOGTABLE` - samples are drawn from a uniform
  distribution and accepted with probability proportional to
  $\exp(-(x-c)²/(2σ²))$ where $\exp(-(x-c)²/(2σ²))$ is computed using
  logarithmically many calls to Bernoulli distributions. Only integer-valued $c$
  are supported.

- `DGS_DISC_GAUSS_UNIFORM_ONLINE` - samples are drawn from a uniform
  distribution and accepted with probability proportional to
  $\exp(-(x-c)²/(2σ²))$ where $\exp(-(x-c)²/(2σ²))$ is computed in each
  invocation. Typically this is very slow. Any real-valued $c$ is accepted.

- `DGS_DISC_SIGMA2_LOGTABLE` - samples are drawn from an easily samplable
  distribution with $σ = k·σ₂$ where $σ₂ := \sqrt{1/(2\log 2)}$ and accepted
  with probability proportional to $\exp(-(x-c)²/(2σ²))$ where
  $\exp(-(x-c)²/(2σ²))$ is computed using logarithmically many calls to
  Bernoulli distributions (but no calls to $\exp$). Note that this sampler
  adjusts sigma to match $σ₂·k$ for some integer $k$.  Only integer-valued $c$
  are supported.

Algorithm 2-4 are described in:

  Léo Ducas, Alain Durmus, Tancrède Lepoint and Vadim Lyubashevsky. *Lattice
  Signatures and Bimodal Gaussians*; in Advances in Cryptology – CRYPTO 2013;
  Lecture Notes in Computer Science Volume 8042, 2013, pp 40-56
  [(PDF)](http://www.di.ens.fr/~lyubash/papers/bimodal.pdf)

## Precisions ##

- `mp` - multi-precision using [MPFR](http://www.mpfr.org/),
  cf. `dgs_gauss_mp.c`

- `dp` - double precision using machine doubles, cf. `dgs_gauss_dp.c`.

For readers unfamiliar with the implemented algorithms it makes sense to start
with ``dgs_gauss_dp.c`` which implements the same algorithms as
``dgs_gauss_mp.c`` but should be easier to read.

## Typical Usage ##

    dgs_disc_gauss_dp_t *D = dgs_disc_gauss_dp_init(<sigma>, <c>, <tau>, <algorithm>);
    D->call(D); // as often as needed
    dgs_disc_gauss_dp_clear(D);
