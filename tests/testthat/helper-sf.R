# Run a block of expectations once per instruction set the host supports, so a
# single CI machine exercises the scalar, SSE2, AVX2 and AVX-512 kernels.
for_each_isa <- function(code) {
  old <- sf_options()$isa
  on.exit(sf_options(isa = old), add = TRUE)
  for (isa in sf_info()$available) {
    sf_options(isa = isa)
    withCallingHandlers(code(isa), error = function(e) {
      stop(sprintf("failure under isa = '%s': %s", isa, conditionMessage(e)),
           call. = FALSE)
    })
  }
}

sf_testdata <- function(n = 5000L, seed = 20260101L) {
  set.seed(seed)
  list(
    n = n,
    a = runif(n, 0.1, 10),
    b = runif(n, 0.1, 10),
    c = runif(n, -5, 5),
    p = runif(n, 0.01, 0.99),
    g = sort(sample.int(max(1L, n %/% 50L), n, replace = TRUE))
  )
}

edge_values <- c(0, -0, 1, -1, 0.5, -0.5, NA_real_, NaN, Inf, -Inf,
                 .Machine$double.xmin, .Machine$double.xmax,
                 .Machine$double.eps, 1e-320, -1e-320, 2^53, -2^53)

# |r^2 - x| computed exactly, as the unevaluated sum of a Dekker two-product:
# Veltkamp splitting gives two half-width parts whose products are each exact,
# and hi - xx is exact by Sterbenz because the two are within a factor of two.
# Normalising by 2^k first keeps r^2 in range for operands as extreme as
# DBL_MAX or the smallest subnormal, and is itself exact. The result is in
# units of 2^(-2k), so a caller comparing two candidates must pass one common
# k for both -- scaling each back separately would underflow to zero for a
# root near 1e-160 and make two unequal residuals compare equal.
# r^2 - x is monotone in |r - sqrt(x)|, so of two candidate doubles the one
# with the smaller residual is the nearer square root.
sq_resid <- function(r, x, k = sq_scale(r)) {
  p  <- 2^k
  rr <- r * p
  xx <- x * p * p
  hi <- rr * rr
  sp <- rr * 134217729                 # 2^27 + 1, Veltkamp's splitting constant
  rh <- sp - (sp - rr)
  rl <- rr - rh
  lo <- ((rh * rh - hi) + 2 * rh * rl) + rl * rl
  abs((hi - xx) + lo)
}

sq_scale <- function(r) ifelse(is.finite(r) & r != 0, -floor(log2(abs(r))), 0)

# Bit-identity with base R's sqrt, with one documented exception. R reaches the
# C library sqrt through a function pointer in math1(), so it always gets the
# library routine; every simdfuse backend uses the hardware square root, which
# IEEE 754 requires to be correctly rounded. On R 4.5.2 for Windows the two
# disagree by one ulp on roughly one argument in two thousand -- in every
# observed case the true root lies within 2e-4 of the midpoint between two
# doubles, the signature of an extended-precision intermediate being rounded
# twice, and the hardware is the value that is right. Such an element is
# accepted only when our result is provably the nearer root; everything else
# must still match base R exactly.
expect_sqrt_identical <- function(got, want, x, info = NULL) {
  d <- which(got != want)              # NA and NaN positions compare as NA
  if (length(d)) {
    k <- sq_scale(want[d])             # one common scale for both candidates
    testthat::expect_true(
      all(sq_resid(got[d], x[d], k) < sq_resid(want[d], x[d], k)),
      info = paste("platform sqrt differs and ours is not the nearer root:", info))
    want[d] <- got[d]
  }
  testthat::expect_identical(got, want, info = info)
}

# Thread counts for the determinism sweeps. CRAN allows a check at most two
# cores, so under R CMD check the ladder stops at two; a local run -- where
# devtools and testthat both set NOT_CRAN -- sweeps up to the real maximum,
# which is the configuration the guarantee actually needs testing in.
# `from` is the lowest thread count the caller's reference is comparable with:
# the parallel scan reassociates `sum`, so its serial path is a different
# answer by design and that sweep starts at two.
sf_test_threads <- function(from = 1L) {
  ks <- if (!identical(Sys.getenv("NOT_CRAN"), "true") ||
            nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")) ||
            nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")))
    c(1L, 2L)
  else
    unique(c(1L, 2L, 3L, 4L, as.integer(sf_info()$max_threads)))
  ks[ks >= from]
}
