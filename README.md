# Discrete Gaussians over the Integers #

[![Build Status](https://travis-ci.com/malb/dgs.svg?branch=master)](https://travis-ci.com/malb/dgs) [![Build Status](https://img.shields.io/bitbucket/pipelines/malb/dgs.svg)](https://bitbucket.org/malb/dgs/addon/pipelines/home) 



A discrete Gaussian distribution on the Integers is a distribution where the integer $x$ is sampled with probability proportional to $exp(-(x-c)²/(2σ²))$. It is denoted by $D_{σ,c}$ where `σ` is the width parameter (close to the standard deviation) and $c$ is the center.

This library samples from this distributions.

**WARNING: This library does not provide constant time implementations. Pull requests welcome.**

## Installation ##

Clone this repository then do

    autoreconf -i
    ./configure
    make
    make check

## Algorithms ##

The library provides two types of algorithms to sample from discrete Gaussians.

### Samplers

This type of algorithm is useful to produce a large number of samples from the same distribution, i.e. with fixed parameters. The parameters need to be provided during initialization. Available algorithms are:

  - `DGS_DISC_GAUSS_UNIFORM_TABLE` - classical rejection sampling, sampling from the uniform distribution and accepted with probability proportional to $\exp(-(x-c)²/(2σ²))$ where $\exp(-(x-c)²/(2σ²))$ is precomputed and stored in a table. Any real-valued `c` is supported.

  - `DGS_DISC_GAUSS_UNIFORM_LOGTABLE` - samples are drawn from a uniform distribution and accepted with probability proportional to $\exp(-(x-c)²/(2σ²))$ where $\exp(-(x-c)²/(2σ²))$ is computed using logarithmically many calls to Bernoulli distributions. Only integer-valued $c$ are supported.

  - `DGS_DISC_GAUSS_UNIFORM_ONLINE` - samples are drawn from a uniform distribution and accepted with probability proportional to $\exp(-(x-c)²/(2σ²))$ where $\exp(-(x-c)²/(2σ²))$ is computed in each invocation. Typically this is very slow. Any real-valued $c$ is accepted.

  - `DGS_DISC_SIGMA2_LOGTABLE` - samples are drawn from an easily samplable distribution with $σ = k·σ₂$ where $σ₂ := \sqrt{1/(2\log 2)}$ and accepted with probability proportional to $\exp(-(x-c)²/(2σ²))$ where $\exp(-(x-c)²/(2σ²))$ is computed using logarithmically many calls to Bernoulli distributions (but no calls to $\exp$). Note that this sampler adjusts sigma to match $σ₂·k$ for some integer $k$. Only integer-valued $c$ are supported.

  - `DGS_DISC_GAUSS_ALIAS` - uses the [alias method](https://en.wikipedia.org/wiki/Alias_method). Setup costs are roughly $σ²$ (as currently implemented) and table sizes linear in $σ$, but sampling is then just a randomized lookup. Any real-valued $c$ is accepted.
  
  - `DGS_DISC_GAUSS_CONVOLUTION` - Applies the convolution technique to alias sampling in order to reduce memory overhead and setup cost at the cost of running time. This is suitable for large $σ$. Any real-valued $c$ is accepted.
    
Algorithm 2-4 are described in:

  Léo Ducas, Alain Durmus, Tancrède Lepoint and Vadim Lyubashevsky. *Lattice Signatures and Bimodal Gaussians*; in Advances in Cryptology – CRYPTO 2013; Lecture Notes in Computer Science Volume 8042, 2013, pp 40-56 [(PDF)](http://www.di.ens.fr/~lyubash/papers/bimodal.pdf)
  
Algorithm 6 is described in:
  
  Thomas Pöppelmann, Léo Ducas, Tim Güneysu. *Enhanced Lattice-Based Signatures on Reconfigurable Hardware*; in Cryptographic Hardware and Embedded Systems – CHES 2014; Lecture Notes in Computer Science Volume 8731, 2014, pp 353-370 [(PDF)](https://eprint.iacr.org/2014/254.pdf)
  
  Daniele Micciancio, Michael Walter. *Gaussian Sampling over the Integers: Efficient, Generic, Constant-Time*; in Advances in Cryptology – CRYPTO 2017; Lecture Notes in Computer Science Volume 10402, 2017, pp 455-485 [(PDF)](https://eprint.iacr.org/2017/259.pdf)

### Randomized rounders

This type of algorithm is useful to produce samples where the parameters of the discrete Gaussian can change for every query, i.e. the center is randomly rounded to a nearby integer. Parameters are provided for every query. Available algorithms are:

  - `DGS_RROUND_UNIFORM_ONLINE` - essentially the same as the sampler `DGS_DISC_GAUSS_UNIFORM_ONLINE`, where the very little parameter dependent precomputation (upper bounds on samples, etc) is now done online.
  
  - `DGS_RROUND_KARNEY` - Karney's algorithm, similar in spirit to the sampler `DGS_DISC_SIGMA2_LOGTABLE`, but without the need to precompute log tables and without restriction on the center.
  
  - `DGS_RROUND_CONVOLUTION` - Reduces the rounding problem to the sampling problem and invokes alias sampling.
  
  Karney's algorithm is described in:
  
  Charles F. F. Karney. *Sampling Exactly from the Normal Distribution*; in ACM Trans. Mathematical Software 42(1), 3:1-14 (Jan. 2016) [(PDF)](https://arxiv.org/pdf/1303.6257)
  
The convolution approach to randomized rounding is described in
  
  Daniele Micciancio, Michael Walter. *Gaussian Sampling over the Integers: Efficient, Generic, Constant-Time*; in Advances in Cryptology – CRYPTO 2017; Lecture Notes in Computer Science Volume 10402, 2017, pp 455-485 [(PDF)](https://eprint.iacr.org/2017/259.pdf)

## Precisions ##

- `mp` - multi-precision using [MPFR](http://www.mpfr.org/), cf. `dgs_gauss_mp.c`

- `dp` - double precision using machine doubles, cf. `dgs_gauss_dp.c`.

For readers unfamiliar with the implemented algorithms it makes sense to start with ``dgs_gauss_dp.c`` which implements the same algorithms as ``dgs_gauss_mp.c`` but should be easier to read.

## Typical usage ##

    dgs_disc_gauss_dp_t *D = dgs_disc_gauss_dp_init(<sigma>, <c>, <tau>, <algorithm>);
    D->call(D); // as often as needed
    dgs_disc_gauss_dp_clear(D);

    dgs_rround_dp_t *D = dgs_rround_dp_init(<tau>, <algorithm>);
    D->call(D, <sigma>, <c>); // as often as needed
    dgs_rround_dp_clear(D);

## Contributors

The following people have contributed to dgs:

- Martin Albrecht
- Shai Halevi
- Michael Walter

## How to cite

	@unpublished{dgs,
	    author = {Martin R. Albrecht and Michael Walter},
	    title = {{dgs}, {D}iscrete {G}aussians over the {I}ntegers},
	    year = 2018,
	    note = {Available at \url{https://bitbucket.org/malb/dgs}},
	    url = {https://bitbucket.org/malb/dgs}
	}
