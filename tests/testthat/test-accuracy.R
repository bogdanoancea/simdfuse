# Accuracy of the optional "fast" math mode, measured rather than assumed.
# The bounds below are what the vectorised polynomials actually achieve; they
# are asserted so that a regression in the kernels shows up as a test failure.

ulp_err <- function(got, want) {
  ok <- is.finite(want) & want != 0 & is.finite(got)
  if (!any(ok)) return(0)
  max(abs(got[ok] - want[ok]) / (2^(floor(log2(abs(want[ok]))) - 52)))
}

test_that("accurate mode reproduces base R bit for bit", {
  set.seed(101)
  x <- runif(50000, -700, 700)
  y <- 10^runif(50000, -300, 300)
  for_each_isa(function(isa) {
    expect_identical(sf_eval(quote(exp(x)), math = "accurate"), exp(x))
    expect_identical(sf_eval(quote(log(y)), math = "accurate"), log(y))
    expect_identical(sf_eval(quote(x * y + sqrt(y)), math = "accurate"),
                     x * y + sqrt(y))
  })
})

test_that("fast exp is within 2 ulp over the whole representable range", {
  set.seed(102)
  for_each_isa(function(isa) {
    for (rng in list(c(-700, 700), c(-1, 1), c(-0.35, 0.35), c(-40, 40))) {
      x <- runif(20000, rng[1], rng[2])
      expect_lt(ulp_err(sf_eval(quote(exp(x)), math = "fast"), exp(x)), 2,
                label = sprintf("exp ulp on [%g,%g] under %s", rng[1], rng[2], isa))
    }
  })
})

test_that("fast log is within 4 ulp, including subnormal inputs", {
  set.seed(103)
  for_each_isa(function(isa) {
    cases <- list(10^runif(20000, -300, 300),
                  runif(20000, 0.5, 2),
                  runif(20000, 1, 1 + 1e-6),
                  10^runif(20000, -320, -308))
    for (i in seq_along(cases)) {
      y <- cases[[i]]
      expect_lt(ulp_err(sf_eval(quote(log(y)), math = "fast"), log(y)), 4,
                label = sprintf("log ulp, case %d under %s", i, isa))
    }
  })
})

test_that("fast math preserves the special values, except NA's payload", {
  # Ordinary finite values are covered by the ulp bounds above; what must be
  # exact is the behaviour at the boundaries of the representable range.
  z <- c(0, NaN, Inf, -Inf, 709.8, 710, -745, -746, -1000)
  for_each_isa(function(isa) {
    expect_identical(sf_eval(quote(exp(z)), math = "fast"), exp(z))
  })
  zf <- c(1, -1, 0.5, 100, -100, 700)
  for_each_isa(function(isa) {
    expect_lt(ulp_err(sf_eval(quote(exp(zf)), math = "fast"), exp(zf)), 2)
  })
  y <- c(0, 1, -1, NaN, Inf, -Inf, 1e-320, .Machine$double.xmin,
         .Machine$double.xmax)
  for_each_isa(function(isa) {
    g <- sf_eval(quote(log(y)), math = "fast")
    w <- suppressWarnings(log(y))
    fin <- is.finite(w)
    expect_lt(ulp_err(g, w), 4)
    expect_identical(g[!fin], w[!fin])
  })

  # Documented deviation: R distinguishes NA from NaN by the payload bits of
  # the quiet NaN. The polynomial kernels arithmetically destroy that payload,
  # so a fast-mode exp() or log() of NA returns NaN. is.na() is TRUE either
  # way, so only code that separates NA from NaN is affected. Accurate mode
  # routes through libm, which preserves the payload.
  na <- NA_real_
  expect_true(is.na(sf_eval(quote(exp(na)), math = "fast")))
  expect_true(is.nan(sf_eval(quote(exp(na)), math = "fast")))
  expect_identical(sf_eval(quote(exp(na)), math = "accurate"), exp(na))
  expect_false(is.nan(sf_eval(quote(exp(na)), math = "accurate")))
})

test_that("fused multiply-add is used only in fast mode", {
  pa <- sf_compile(~ a * b + c, math = "accurate")
  pf <- sf_compile(~ a * b + c, math = "fast")
  expect_false("FMA3" %in% sf_disassemble(pa)$op)
  expect_true("FMA3" %in% sf_disassemble(pf)$op)
  # accurate mode must keep the two separate roundings that base R performs
  set.seed(104)
  a <- runif(10000); b <- runif(10000); c <- runif(10000)
  expect_identical(sf_eval(quote(a * b + c), math = "accurate"), a * b + c)
})

test_that("integer exponents use pow() in accurate mode", {
  pa <- sf_compile(~ x^3, math = "accurate")
  pf <- sf_compile(~ x^3, math = "fast")
  expect_true("POW" %in% sf_disassemble(pa)$op)
  expect_true("POWI" %in% sf_disassemble(pf)$op)
  # x^2 is the one case base R itself turns into a multiplication
  expect_true("POWI" %in% sf_disassemble(sf_compile(~ x^2, math = "accurate"))$op)
})

test_that("%% and %/% match base R wherever base R is reliable", {
  set.seed(105)
  a <- runif(30000, -50, 50); b <- runif(30000, -50, 50)
  for_each_isa(function(isa) {
    expect_identical(sf_eval(quote(a %% b)), a %% b)
    expect_identical(sf_eval(quote(a %/% b)), a %/% b)
  })
  v <- c(5, -5, 0, Inf, -Inf, NA, NaN, 7, -7, 2, -2, 0.1, 1e16, -1e16)
  g <- expand.grid(x = v, y = v)
  x <- g$x; y <- g$y
  want <- suppressWarnings(mapply(function(p, q) p %/% q, x, y))
  expect_identical(sf_eval(quote(x %/% y)), want)
})

test_that("%% is exact where base R warns of accuracy loss", {
  # R computes x %% y as x - floor(x/y) * y, which carries no information once
  # |x/y| exceeds 1/eps; R emits "probable complete loss of accuracy in
  # modulus" and returns a meaningless value. simdfuse routes through fmod(),
  # which is exact by definition, and so returns the true remainder.
  expect_warning(r_val <- 1e300 %% 7, "loss of accuracy")
  x <- c(1e300, 1e300); y <- c(7, 7)
  sf_val <- sf_eval(quote(x %% y))[1]

  # independent reference: 1e300 is an exact integer, so long-divide its digits
  digits <- as.numeric(strsplit(format(1e300, scientific = FALSE), "")[[1]])
  ref <- 0
  for (d in digits) ref <- (ref * 10 + d) %% 7

  expect_identical(sf_val, ref)
  expect_false(isTRUE(all.equal(r_val, ref)))
})

test_that("every backend returns the same bits for the same expression", {
  # Without this the dispatcher would make answers depend on the host CPU.
  set.seed(106)
  a <- runif(20000, -30, 30); b <- runif(20000, 0.1, 30)
  exprs <- list(quote(a %% b), quote(a %/% b), quote(a * b + sqrt(b)),
                quote((a - 1) / b), quote(exp(a / 30)), quote(log(b)),
                quote(pmin(a, b)), quote(ifelse(a > 0, a, b)))
  old <- sf_options()$isa
  on.exit(sf_options(isa = old), add = TRUE)
  for (ex in exprs) {
    sf_options(isa = "scalar")
    ref <- sf_eval(ex)
    for (isa in sf_info()$available) {
      sf_options(isa = isa)
      expect_identical(sf_eval(ex), ref,
                       info = sprintf("%s under %s", deparse(ex), isa))
    }
  }
})

test_that("fast math is gated per function and falls back exactly", {
  old <- sf_options()
  on.exit(sf_options(fast_exp = old$fast_exp, fast_log = old$fast_log))
  set.seed(11)
  x <- runif(5000, 0.5, 20)

  # Switched off, fast mode must reach the same libm accurate mode does, and
  # therefore agree with base R bit for bit -- including the NA payload that
  # the polynomial path destroys.
  sf_options(fast_exp = FALSE, fast_log = FALSE)
  expect_identical(sf_eval(exp(x), math = "fast"), exp(x))
  expect_identical(sf_eval(log(x), math = "fast"), log(x))
  na <- c(NA_real_, 1)
  expect_false(is.nan(sf_eval(exp(na), math = "fast"))[1])

  # Switched on, within the documented bound but no longer bit-identical.
  sf_options(fast_exp = TRUE, fast_log = TRUE)
  expect_lt(max(abs(sf_eval(exp(x), math = "fast") - exp(x)) /
                (.Machine$double.eps * abs(exp(x)))), 2)
  expect_lt(max(abs(sf_eval(log(x), math = "fast") - log(x)) /
                (.Machine$double.eps * pmax(abs(log(x)), 1))), 4)

  # The two flags are independent.
  sf_options(fast_exp = TRUE, fast_log = FALSE)
  expect_identical(sf_eval(log(x), math = "fast"), log(x))
  expect_false(identical(sf_eval(exp(x), math = "fast"), exp(x)))

  # Accurate mode ignores them entirely.
  for (fe in c(TRUE, FALSE)) {
    sf_options(fast_exp = fe, fast_log = fe)
    expect_identical(sf_eval(exp(x) + log(x), math = "accurate"), exp(x) + log(x))
  }
})
