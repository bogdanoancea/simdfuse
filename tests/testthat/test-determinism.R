test_that("parallel reductions are deterministic in the thread count", {
  # This is the guarantee a plain `omp parallel for reduction(+:)` cannot make:
  # task boundaries depend only on n, and partials combine in index order.
  n <- 1e6L
  set.seed(7)
  x <- rnorm(n) * 10^runif(n, -8, 8)   # wide exponent range, worst case for summation
  y <- runif(n)
  ref_sum  <- sf_sum(quote(x * y + sqrt(y)), threads = 1L)
  ref_mean <- sf_mean(quote(x * y), threads = 1L)
  ref_prod <- sf_prod(quote(1 + y / 1e9), threads = 1L)

  for (th in sf_test_threads()) {
    expect_identical(sf_sum(quote(x * y + sqrt(y)), threads = th), ref_sum,
                     info = paste("threads =", th))
    expect_identical(sf_mean(quote(x * y), threads = th), ref_mean,
                     info = paste("threads =", th))
    expect_identical(sf_prod(quote(1 + y / 1e9), threads = th), ref_prod,
                     info = paste("threads =", th))
  }
})

test_that("elementwise results are deterministic in threads and block size", {
  # This expression calls log() and exp() per element, so it parallelises from
  # 2 tasks (131072 elements) up; n is comfortably above that.
  n <- 3e5L
  set.seed(8)
  x <- rnorm(n); y <- runif(n, 0.1, 5)
  ref <- sf_eval(quote(x * y + log(y) - exp(-x)), threads = 1L, block = 512L)
  for (th in sf_test_threads())
    for (blk in c(64L, 512L, 4096L))
      expect_identical(sf_eval(quote(x * y + log(y) - exp(-x)),
                               threads = th, block = blk), ref,
                       info = paste("threads", th, "block", blk))
})

test_that("scan is deterministic in the thread count", {
  # `x` alone is a streaming expression, so it needs 4 tasks per worker --
  # 8 tasks, 524288 elements -- before the engine will parallelise at all.
  # Below that the scan runs serially whatever the thread count and this test
  # would pass vacuously, comparing a value with itself.
  n <- 2e6L
  set.seed(9)
  x <- rnorm(n)
  ref <- sf_scan(quote(x), "sum", threads = 2L)
  for (th in sf_test_threads(2L))
    expect_identical(sf_scan(quote(x), "sum", threads = th), ref)
  # max and min are order-independent, so the parallel scan reproduces base R
  expect_identical(sf_scan(quote(x), "max"), cummax(x))
  expect_identical(sf_scan(quote(x), "min"), cummin(x))
})

test_that("reductions stay close to base R despite a different summation order", {
  # Blocked summation is not bit-identical to R's sequential long-double
  # accumulation; it should be at least as accurate, not merely close.
  set.seed(10)
  x <- rnorm(1e6) * 10^runif(1e6, -10, 10)
  exact <- sum(as.numeric(sort(x)))          # not exact, but a neutral referee
  r_err  <- abs(sum(x) - exact)
  sf_err <- abs(sf_sum(quote(x)) - exact)
  expect_equal(sf_sum(quote(x)), sum(x), tolerance = 1e-10)
  expect_lt(sf_err, max(r_err * 100, 1e-6) )
})
