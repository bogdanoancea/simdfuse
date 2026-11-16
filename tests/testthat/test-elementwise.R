test_that("elementwise results are bit-identical to base R on every ISA", {
  d <- sf_testdata()
  a <- d$a; b <- d$b; cc <- d$c; p <- d$p
  mu <- 2.5; s <- 1.25

  cases <- list(
    quote(a + b),                 quote(a - b),
    quote(a * b),                 quote(a / b),
    quote(-a),                    quote(a * b + cc),
    quote(a * b + sqrt(a) - 2.5), quote((a - mu) / s),
    quote(sqrt(a)),               quote(exp(cc)),
    quote(log(a)),                quote(log(a, 2)),
    quote(log(a, 10)),            quote(log(a, b)),
    quote(log1p(a)),              quote(expm1(cc)),
    quote(log2(a)),               quote(log10(a)),
    quote(abs(cc)),               quote(sign(cc)),
    quote(a^2),                   quote(a^3),
    quote(a^-2),                  quote(a^0),
    quote(a^1),                   quote(a^2.5),
    quote(a^b),                   quote(sin(cc)),
    quote(cos(cc)),               quote(tan(cc)),
    quote(atan(cc)),              quote(atan2(a, b)),
    quote(asin(p)),               quote(acos(p)),
    quote(sinh(cc)),              quote(cosh(cc)),
    quote(tanh(cc)),              quote(lgamma(a)),
    quote(gamma(p)),              quote(digamma(a)),
    quote(floor(cc)),             quote(ceiling(cc)),
    quote(trunc(cc)),             quote(round(cc)),
    quote(round(cc, 2)),          quote(a %% 3),
    quote(a %/% 3),               quote(a %% b),
    quote(a %/% b),
    quote(pmin(a, b)),            quote(pmax(a, b)),
    quote(pmin(a, b, 2)),         quote(ifelse(a > b, a, b)),
    quote(dnorm(cc)),             quote(dnorm(cc, mu, s)),
    quote(pnorm(cc)),             quote(pnorm(cc, mu, s)),
    quote(plogis(cc)),            quote(qlogis(p)),
    quote((a - mu)^2 + (a - mu)),
    quote(((a - mu) / s)^2 / 2 + log(s)),
    quote(exp(-0.5 * ((a - mu) / s)^2) / (s * sqrt(2 * pi)))
  )

  for_each_isa(function(isa) {
    # sqrt is asserted on its own terms, because the C library sqrt that base R
    # calls is not correctly rounded on every platform (see helper-sf.R). The
    # reference for the composite cases then uses the same square root simdfuse
    # produces, so those cases test the arithmetic around it rather than the
    # library. Constants are left to base R: sqrt(2 * pi) sits 0.023 ulp from
    # the midpoint between two doubles, far enough out that no plausible
    # extended-precision intermediate rounds it the other way.
    expect_sqrt_identical(sf_eval(quote(sqrt(a))), sqrt(a), a, info = "sqrt(a)")
    sq <- function(v) if (length(v) > 1L) sf_eval(quote(sqrt(v))) else base::sqrt(v)
    env <- new.env(parent = environment())
    env$sqrt <- sq
    for (ex in cases) {
      want <- eval(ex, env)
      got <- sf_eval(ex)
      expect_identical(got, want, info = paste(deparse(ex), "under", isa))
    }
  })
})

test_that("comparison and logical expressions return logicals", {
  a <- c(1, 2, 3, NA); b <- c(3, 2, 1, 1)
  expect_identical(sf_eval(quote(a > b)), a > b)
  expect_identical(sf_eval(quote(a >= b)), a >= b)
  expect_identical(sf_eval(quote(a == b)), a == b)
  expect_identical(sf_eval(quote(a != b)), a != b)
  expect_identical(sf_eval(quote(!(a > b))), !(a > b))
  expect_identical(sf_eval(quote(is.na(a))), is.na(a))
  expect_identical(sf_eval(quote(is.finite(a))), is.finite(a))
  expect_identical(sf_eval(quote(is.infinite(a))), is.infinite(a))
})

test_that("three-valued logic matches R", {
  grid <- expand.grid(x = c(1, 0, NA), y = c(1, 0, NA))
  x <- grid$x; y <- grid$y
  lx <- as.logical(x); ly <- as.logical(y)
  for_each_isa(function(isa) {
    expect_identical(sf_eval(quote(x & y)), lx & ly)
    expect_identical(sf_eval(quote(x | y)), lx | ly)
    expect_identical(sf_eval(quote(!x)), !lx)
  })
})

test_that("edge and special values propagate as in base R", {
  z <- edge_values
  w <- edge_values * 0.5
  for_each_isa(function(isa) {
    expect_identical(sf_eval(quote(z + w)), z + w)
    expect_identical(sf_eval(quote(z * w)), z * w)
    expect_identical(sf_eval(quote(z / w)), z / w)
    expect_sqrt_identical(sf_eval(quote(sqrt(z))), suppressWarnings(sqrt(z)), z,
                          info = "sqrt(edge values)")
    expect_identical(sf_eval(quote(abs(z))), abs(z))
    expect_identical(sf_eval(quote(z^2)), z^2)
    expect_identical(sf_eval(quote(pmin(z, w))), pmin(z, w))
    expect_identical(sf_eval(quote(pmax(z, w))), pmax(z, w))
    expect_identical(sf_eval(quote(ifelse(z > 0, z, w))), ifelse(z > 0, z, w))
  })
  expect_identical(sf_eval(quote(log(z))), suppressWarnings(log(z)))
})

test_that("scalars, attributes and input types behave", {
  a <- runif(100); k <- 3
  expect_identical(sf_eval(quote(a * k)), a * k)

  m <- matrix(runif(100), 10)
  expect_identical(dim(sf_eval(quote(m * 2))), c(10L, 10L))
  expect_identical(sf_eval(quote(m * 2)), m * 2)

  v <- c(first = 1, second = 2, third = 3)
  expect_identical(names(sf_eval(quote(v * 2))), names(v))

  i <- 1:100
  expect_identical(sf_eval(quote(i * 2)), as.double(i) * 2)
  l <- rep(c(TRUE, FALSE), 50)
  expect_identical(sf_eval(quote(l + 1)), as.double(l) + 1)
})

test_that("block size does not change results", {
  d <- sf_testdata(3000L)
  a <- d$a; b <- d$b
  ref <- sf_eval(quote(a * b + sqrt(a)), block = 512L)
  for (blk in c(8L, 16L, 64L, 127L, 1024L, 4096L, 100000L))
    expect_identical(sf_eval(quote(a * b + sqrt(a)), block = blk), ref)
})

test_that("vector lengths that are not block multiples are handled", {
  for (n in c(1L, 2L, 3L, 7L, 8L, 9L, 15L, 16L, 17L, 511L, 512L, 513L, 1023L)) {
    a <- seq_len(n) + 0.5
    expect_identical(sf_eval(quote(a * 2 + 1)), a * 2 + 1, info = paste("n =", n))
  }
})
