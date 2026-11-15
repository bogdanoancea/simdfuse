
test_that("integer and logical inputs are read in place, with base R's semantics", {
  set.seed(31)
  n <- 3000
  iv <- sample(c(-5L:5L, NA_integer_), n, TRUE)
  lv <- sample(c(TRUE, FALSE, NA), n, TRUE)
  dv <- runif(n, 0.5, 2)

  for_each_isa(function(isa) {
    # identical, not merely equal: NA must stay NA rather than becoming NaN or
    # the double value of INT_MIN
    expect_identical(sf_eval(iv * 2 + 1), as.double(iv) * 2 + 1)
    expect_identical(sf_eval(lv + 0), as.double(lv) + 0)
    expect_identical(sf_eval(iv * dv + lv), as.double(iv) * dv + as.double(lv))
    expect_identical(is.na(sf_eval(iv + 0)), is.na(as.double(iv)))
    expect_false(any(is.nan(sf_eval(iv + 0))))

    # integer scalars, which take the broadcast path rather than the block path
    expect_identical(sf_eval(dv * 3L), dv * 3)
    expect_identical(sf_eval(dv + TRUE), dv + 1)

    # reductions over integer storage
    expect_identical(sf_sum(quote(iv), na.rm = TRUE), sum(as.double(iv), na.rm = TRUE))
    expect_identical(sf_max(quote(iv), na.rm = TRUE), max(as.double(iv), na.rm = TRUE))
    expect_identical(sf_sum(quote(lv), na.rm = TRUE), sum(as.double(lv), na.rm = TRUE))

    # INT_MIN is NA_integer_, so it must never appear as a finite double
    edge <- c(NA_integer_, .Machine$integer.max, -.Machine$integer.max, 0L, 1L)
    expect_identical(sf_eval(edge + 0), as.double(edge) + 0)
  })

  # the input vector is not modified or coerced in place
  before <- iv
  invisible(sf_eval(iv * 2))
  expect_identical(iv, before)
  expect_true(is.integer(iv))

  # a block boundary landing mid-vector must convert the tail correctly
  old <- sf_options()$block
  on.exit(sf_options(block = old))
  for (b in c(64L, 128L)) {
    sf_options(block = b)
    k <- as.integer(b) + 7L          # deliberately not a multiple of the block
    jv <- sample(c(1L:9L, NA_integer_), k, TRUE)
    expect_identical(sf_eval(jv * 2), as.double(jv) * 2)
  }
})
