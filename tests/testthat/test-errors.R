test_that("unsupported constructs fail at compile time with a useful message", {
  x <- runif(10); y <- runif(10)
  expect_error(sf_eval(quote(myfun(x))), "not supported")
  expect_error(sf_eval(quote(sum(x))), "is a reduction")
  expect_error(sf_eval(quote(cumsum(x))), "is a reduction")
  expect_error(sf_eval(quote(x[1])), "not supported")
  expect_error(sf_eval(quote(sqrt(x, y))), "takes one argument")
  expect_error(sf_eval(quote(log(x, base = 2))), "named arguments")
  expect_error(sf_eval(quote(is.nan(x))), "not supported")
  expect_error(sf_compile(~ 1 + 2), "no variables")
})

test_that("length mismatches are errors, not silent recycling", {
  a <- runif(10); b <- runif(3)
  expect_error(sf_eval(quote(a + b)), "does not recycle")
  z <- numeric(0)
  expect_error(sf_eval(quote(z + 1)), "length 0")
})

test_that("unsupported input types are rejected clearly", {
  s <- letters[1:5]
  expect_error(sf_eval(quote(s + 1)), "only double, integer and logical")
  f <- factor(c("a", "b"))
  expect_error(sf_eval(quote(f + 1)), "factor")
  cx <- complex(real = 1:4, imaginary = 1)
  expect_error(sf_eval(quote(cx + 1)), "only double, integer and logical")
  expect_error(sf_eval(quote(nope + 1)), "not found")
})

test_that("options validate their arguments", {
  expect_error(sf_options(block = 0), "positive integer")
  expect_error(sf_options(threads = -1), ">= 0")
  expect_error(sf_options(math = "exact"))
  expect_error(sf_options(isa = "avx1024"))
})

# The original form of this test only exercised the fallback path on a CPU that
# happens to lack AVX-512, and asserted nothing about the honoured path. The
# version below checks the contract on every platform: `forced` is TRUE exactly
# when the requested backend is the one in use. That is what caught the aarch64
# branch reporting forced = TRUE after it had warned and fallen back.
test_that("an unavailable instruction set warns and falls back to auto", {
  old <- sf_options()$isa
  on.exit(sf_options(isa = old), add = TRUE)
  sf_options(isa = "auto")
  auto_isa <- sf_info()$isa
  absent <- setdiff(c("avx512", "avx2", "sse2", "neon", "scalar"),
                    sf_info()$available)
  skip_if(length(absent) == 0L, "every backend is available on this CPU")
  expect_warning(sf_options(isa = absent[1L]),
                 "not available|only the scalar")
  expect_false(sf_info()$forced)
  expect_identical(sf_info()$isa, auto_isa)
})

test_that("an available instruction set is pinned and reported as forced", {
  old <- sf_options()$isa
  on.exit(sf_options(isa = old), add = TRUE)
  for (isa in sf_info()$available) {
    expect_silent(sf_options(isa = isa))
    expect_true(sf_info()$forced, info = isa)
    # the invariant: forced means the active backend is the requested one
    expect_true(startsWith(sf_info()$isa, sub("avx512", "avx512f", isa)),
                info = isa)
  }
  sf_options(isa = "auto")
  expect_false(sf_info()$forced)
})

test_that("compiled programs are reused from the cache", {
  x <- runif(100)
  p1 <- .sf_program_for(quote(x * 2 + 1), FALSE)
  p2 <- .sf_program_for(quote(x * 2 + 1), FALSE)
  expect_identical(p1, p2)
})

test_that("a failed instruction-set request lands on the automatic choice", {
  # Regression guard. The engine warns when it cannot honour a request, and
  # Rf_warning unwinds the C frame whenever the caller has a handler attached
  # (expect_warning, tryCatch, options(warn = 2)). Warning before the fallback
  # pick therefore skipped it and left the previously pinned backend active,
  # so the outcome depended on how the caller handled warnings. It must not.
  old <- sf_options()$isa
  on.exit(sf_options(isa = old), add = TRUE)
  sf_options(isa = "auto")
  auto_isa <- sf_info()$isa
  absent <- setdiff(c("avx512", "avx2", "sse2", "neon", "scalar"),
                    sf_info()$available)
  skip_if(length(absent) == 0L, "every backend is available on this CPU")
  skip_if(identical(auto_isa, "scalar"), "no SIMD backend to contrast with")

  for (mode in c("captured", "suppressed", "escalated")) {
    sf_options(isa = "scalar")                 # deliberately not the auto pick
    expect_identical(sf_info()$isa, "scalar")
    switch(mode,
      captured   = expect_warning(sf_options(isa = absent[1L])),
      suppressed = suppressWarnings(sf_options(isa = absent[1L])),
      escalated  = {
        op <- options(warn = 2)
        try(sf_options(isa = absent[1L]), silent = TRUE)
        options(op)
      })
    expect_identical(sf_info()$isa, auto_isa, info = mode)
    expect_false(sf_info()$forced, info = mode)
    expect_identical(sf_options()$isa, "auto", info = mode)
  }
})
