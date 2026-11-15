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
  n <- 2e6L
  x <- runif(n)
  mu <- mean(x)
  peak <- function(f) {
    gc(reset = TRUE, full = TRUE)
    invisible(f())
    unname(gc(full = TRUE)["Vcells", "max used"])
  }
  p_base <- peak(function() sum((x - mu)^2 + sqrt(x)))
  p_sf   <- peak(function() sf_sum(quote((x - mu)^2 + sqrt(x))))
  # Base R materialises three n-length temporaries; the fused reduction holds
  # only its register file, which is a few kilobytes.
  expect_lt(p_sf, p_base * 0.75)
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
