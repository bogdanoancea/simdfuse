test_that("fused reductions agree with base R", {
  d <- sf_testdata()
  a <- d$a; b <- d$b; cc <- d$c
  mu <- 2.5
  for_each_isa(function(isa) {
    expect_equal(sf_sum(quote(a)), sum(a))
    expect_equal(sf_sum(quote((a - mu)^2)), sum((a - mu)^2))
    expect_equal(sf_mean(quote(a * b)), mean(a * b))
    expect_identical(sf_max(quote(a * b)), max(a * b))
    expect_identical(sf_min(quote(a * b)), min(a * b))
    expect_equal(sf_prod(quote(1 + a / 1e6)), prod(1 + a / 1e6))
    expect_identical(sf_any(quote(a > 9.9)), any(a > 9.9))
    expect_identical(sf_all(quote(a > 0)), all(a > 0))
    expect_identical(sf_nobs(quote(a)), as.double(length(a)))
  })
})

test_that("reductions handle NA, NaN and Inf as base R does", {
  z <- c(1, 2, NA, 4, NaN, Inf, -Inf, 0)
  for_each_isa(function(isa) {
    expect_identical(sf_sum(quote(z)), sum(z))
    expect_identical(sf_sum(quote(z), na.rm = TRUE), sum(z, na.rm = TRUE))
    expect_identical(sf_max(quote(z)), max(z))
    expect_identical(sf_max(quote(z), na.rm = TRUE), max(z, na.rm = TRUE))
    expect_identical(sf_min(quote(z)), min(z))
    expect_identical(sf_min(quote(z), na.rm = TRUE), min(z, na.rm = TRUE))
    expect_identical(sf_mean(quote(z), na.rm = TRUE), mean(z, na.rm = TRUE))
    expect_identical(sf_nobs(quote(z)), sum(!is.na(z)) + 0)
  })

  w <- c(1, 0, NA)
  expect_identical(sf_any(quote(w > 0.5)), any(w > 0.5))
  expect_identical(sf_all(quote(w > 0.5)), all(w > 0.5))
  expect_identical(sf_any(quote(w > 0.5), na.rm = TRUE), any(w > 0.5, na.rm = TRUE))
  expect_identical(sf_all(quote(w > 0.5), na.rm = TRUE), all(w > 0.5, na.rm = TRUE))

  allna <- c(NA_real_, NA_real_)
  expect_identical(sf_sum(quote(allna), na.rm = TRUE), 0)
  expect_warning(v <- sf_max(quote(allna), na.rm = TRUE))
  expect_identical(v, -Inf)
})

test_that("a fused reduction allocates no vector of results", {
  # A resource threshold, not a correctness property, so it is not worth a
  # failure on one of CRAN's flavours over an allocator difference.
  skip_on_cran()

  n <- 2e6L
  x <- runif(n)
  mu <- mean(x)

  # gc()'s "max used" is a high-water mark for the whole session, so comparing
  # two of them measures the session's baseline as much as the allocation under
  # test: with baseline B and temporaries T, p_base = B + T against p_sf = B
  # only satisfies p_sf < 0.75 * p_base when T > B/3. That held in a fresh
  # session and failed on R 4.0 with more loaded. Measure the increment.
  peak_delta <- function(f) {
    gc(reset = TRUE, full = TRUE)
    before <- unname(gc(full = TRUE)["Vcells", "used"])
    invisible(f())
    unname(gc(full = TRUE)["Vcells", "max used"]) - before
  }

  d_base <- peak_delta(function() sum((x - mu)^2 + sqrt(x)))
  d_sf   <- peak_delta(function() sf_sum(quote((x - mu)^2 + sqrt(x))))

  # One n-length double vector is n Vcells. Base R materialises three of them;
  # the fused reduction holds only its register file, a few kilobytes, so it
  # must stay well under a single vector's worth.
  expect_true(d_base > n,
    info = sprintf("base R allocated %.3g Vcells for n = %d; expected > n",
                   d_base, n))
  expect_true(d_sf < n / 4,
    info = sprintf("simdfuse allocated %.3g Vcells (base R %.3g) for n = %d; expected < n/4",
                   d_sf, d_base, n))
})

test_that("the compiled-kernel form matches the one-shot form", {
  x <- rnorm(1000); m <- 0.1; s <- 1.3
  k <- sf_kernel(~ 0.5 * ((x - m) / s)^2 + log(s), reduce = "sum")
  expect_identical(names(formals(k)), c("x", "m", "s"))
  expect_identical(k(x = x, m = m, s = s),
                   sf_sum(~ 0.5 * ((x - m) / s)^2 + log(s)))
  ke <- sf_kernel(~ x * 2)
  expect_identical(ke(x = x), x * 2)
})
