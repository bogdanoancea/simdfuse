# simdfuse

Fused SIMD evaluation of array expressions and parallel skeletons for R.

`a * b + sqrt(c)` in base R allocates three temporary vectors of length `n` and
walks memory four times. `sf(a * b + sqrt(c))` compiles the expression once and
walks the data a single time, keeping every intermediate inside a
cache-resident block. The saving grows with the number of operations and with
`n`: once the operands no longer fit in cache, ordinary evaluation is limited by
memory bandwidth, and removing traversals is the only thing that helps.

The expression is compiled to a compact instruction stream and executed by a
kernel built once per instruction set -- AVX-512, AVX2 with FMA, SSE2, NEON and
a scalar fallback -- and selected at run time from the host's reported features.
One binary vectorises on whatever hardware it finds.

## Installation

```r
# install.packages("remotes")
remotes::install_github("bogdanoancea/simdfuse")
```

A C compiler is required. OpenMP is used when available and the package works
without it, single-threaded.

## Fused elementwise evaluation

```r
library(simdfuse)

a <- runif(1e7); b <- runif(1e7); c <- runif(1e7)

z <- sf(a * b + sqrt(c) - 2.5)
identical(z, a * b + sqrt(c) - 2.5)

# standardising without temporaries
mu <- mean(a); s <- sd(a)
sf((a - mu) / s)
```

Variables are looked up in `data` first, then in the calling frame. Only the
functions listed by `sf_supported()` are allowed; anything else is an error
rather than a silent fallback to base R. `dim`, `dimnames` and `names` carry
over from the first non-scalar input, so a matrix in gives a matrix out.

## Fused reductions

Reductions consume an expression rather than a materialised vector, so they run
in constant auxiliary memory -- there is no length-`n` temporary anywhere.

```r
x <- rnorm(1e7); mu <- mean(x); s <- sd(x)

sf_sum(~ (x - mu)^2)

# a Gaussian log-likelihood with no intermediate vectors at all
sf_sum(~ -0.5 * log(2 * pi) - log(s) - 0.5 * ((x - mu) / s)^2)
```

`sf_mean()`, `sf_max()`, `sf_min()`, `sf_prod()`, `sf_any()`, `sf_all()` and
`sf_nobs()` follow the same shape.

The parallel form is **deterministic in the thread count**. Task boundaries
depend only on `n`, never on how many threads happen to be available, and task
partials combine serially in index order in extended precision. The answer on
one thread and on sixty is the same bits -- a guarantee a plain
`omp parallel for reduction(+:)` cannot make.

## Parallel skeletons

Constructs base R has no vectorised equivalent for, each consuming an
expression rather than a materialised vector.

```r
# associative prefix scan, including a numerically stable log-domain scan
x <- c(3, 1, 4, 1, 5)
sf_scan(~ x, "sum")
identical(sf_scan(~ x, "sum", parallel = FALSE), cumsum(x))

lp <- rnorm(10) * 50
sf_scan(~ lp, "logsumexp")        # a stable running log-normaliser

# weighted stencil, with a choice of boundary rule
y <- cumsum(rnorm(100))
sf_stencil(~ y, rep(1/5, 5))
sf_stencil(~ log(abs(y) + 1), c(-0.5, 0, 0.5), boundary = "replicate")

# reduction within contiguous segments
val  <- c(1, 2, 3, 10, 20, 100)
firm <- c(1, 1, 1,  2,  2,   3)
sf_segreduce(~ val^2, firm, "sum")
```

## Compiling once, calling many times

Inside an optimiser the same expression is evaluated thousands of times with
different parameters. `sf_kernel()` compiles it once and returns a function
whose formals are the expression's free variables.

```r
nll <- sf_kernel(~ 0.5 * ((x - mu) / s)^2 + log(s), reduce = "sum")
names(formals(nll))        # the expression's free variables

x <- rnorm(1e6)
optim(c(0, 1), function(p) nll(x = x, mu = p[1], s = exp(p[2])))
```

`sf_compile()` and `sf_disassemble()` expose the instruction stream if you want
to see what the compiler produced.

## What is guaranteed

In the default `math = "accurate"` mode every elementwise expression returns
the same bits as base R. Elementary functions go through the same platform
`libm` R calls, one lane at a time; `pnorm`, `dnorm`, `plogis` and `qlogis` call
R's own routines rather than being rebuilt from `erfc`.

There is one documented exception. `sqrt` is taken from the hardware on every
backend, which IEEE 754 requires to be correctly rounded, whereas base R
reaches `sqrt` through a function pointer and so always gets the C library
routine -- and not every C library is correctly rounded. On R 4.5.2 for Windows
the two disagree by one ulp on roughly one argument in two thousand, and the
hardware value is the nearer of the two. simdfuse returns the correctly rounded
value there, which means it is the more accurate of the pair and not
bit-identical to base R on that platform.

`math = "fast"` opts into vectorised polynomial `exp` and `log` with measured
error bounds (2 ulp and 4 ulp respectively), per function, and only where
measurement showed the polynomial actually wins.

Reductions promise invariance in the thread count, not bit-identity with
`sum()`: base R accumulates sequentially in long double, this is blocked
summation, and the two differ in the last bits. `vignette("accuracy")` sets out
every guarantee and every deliberate disagreement, with the figures computed on
your machine when the vignette is built.

## Backends and tuning

```r
str(sf_info())          # which kernel resolved here, lanes, OpenMP, threads
sf_supported()          # which functions may appear in an expression
sf_options()            # block size, thread count, math mode
sf_tune()               # measure this machine and pick the block size for it
```

`sf_options(isa = "sse2")` pins a narrower backend, which is useful for testing
and for reproducing results across heterogeneous machines.

## Measured performance

Speedup over base R depends on how many operations the expression contains, on
the lane width, and on whether the problem is bandwidth-bound. These are
measured figures from two machines, not projections:

| | Apple Silicon, NEON, 2 lanes, 10 cores | AMD EPYC 7713, AVX2, 4 lanes, 64 cores |
|---|---|---|
| elementwise, growing expressions | 1.40x at 2 operations, 2.58x at 8, 3.06x at 16 | 0.94x to 3.21x, rising with operation count |
| reductions, one thread | 2.78x to 3.47x | -- |
| reductions, all threads | 6.5x to 12.8x | 36x to 114x at n = 5e8, at 146-159 GB/s |
| scalar to AVX2, compute-bound `exp` | -- | 5.89x |
| scalar to AVX2, bandwidth-bound | -- | 1.03x |

Two things are worth reading off that table rather than glossed over. A
single-operator expression gains nothing from fusion and can lose slightly:
there is no intermediate to eliminate, and the block machinery is not free.
And vector width only helps when the work is compute-bound -- the same
instruction set gives 5.89x on an `exp`-heavy expression and 1.03x when the
expression is limited by memory bandwidth. On a dual-socket machine, NUMA
placement matters more than either: 93 GB/s unbound against 160 GB/s bound to
one socket, which the package cannot fix because R allocates the vectors
single-threaded before simdfuse is ever called.

`vignette("performance")` measures all of this on the machine that builds it,
so the numbers describe your hardware rather than ours.

## Limitations

- Expressions are elementwise only. No matrix algebra, no sorting, no joins.
- Inputs must all be the same length, or length 1. Partial recycling is an
  error, not a warning.
- Only `sf_supported()` functions are allowed, by design: a silent fallback to
  base R would hide the fact that nothing was fused.
- Gains need `n` large enough to matter. At small `n` the compile and dispatch
  overhead dominates; the expression cache amortises the compile across
  repeated calls, and `sf_kernel()` removes it from the loop entirely.
- Reductions are not bit-identical to base R's, by construction.

## License

GPL-3. Bug reports and feature requests at
<https://github.com/bogdanoancea/simdfuse/issues>.
