test_that("serial scan reproduces base R exactly", {
  for (n in c(1L, 2L, 17L, 513L, 5000L)) {
    set.seed(n)
    x <- rnorm(n)
    expect_identical(sf_scan(quote(x), "sum", parallel = FALSE), cumsum(x))
    expect_identical(sf_scan(quote(x), "max", parallel = FALSE), cummax(x))
    expect_identical(sf_scan(quote(x), "min", parallel = FALSE), cummin(x))
    expect_equal(sf_scan(quote(x), "prod", parallel = FALSE), cumprod(x))
  }
})

test_that("scan propagates NA like base R", {
  x <- c(1, 2, NA, 4, 5)
  expect_identical(sf_scan(quote(x), "sum", parallel = FALSE), cumsum(x))
  expect_identical(sf_scan(quote(x), "max", parallel = FALSE), cummax(x))
  expect_identical(sf_scan(quote(x), "min", parallel = FALSE), cummin(x))
})

test_that("log-domain scan is stable where the naive form overflows", {
  lp <- c(800, 801, 802)          # exp(800) overflows a double
  got <- sf_scan(quote(lp), "logsumexp", parallel = FALSE)
  want <- c(800,
            800 + log1p(exp(1)),
            802 + log(exp(-2) + exp(-1) + 1))
  expect_equal(got, want)
  expect_true(all(is.finite(got)))
  expect_true(all(is.infinite(log(cumsum(exp(lp))))))
})

test_that("scan fuses with the expression", {
  x <- runif(1000, 0.5, 2)
  expect_equal(sf_scan(~ log(x), "sum", parallel = FALSE), cumsum(log(x)))
  expect_equal(sf_scan(~ x^2, "max", parallel = FALSE), cummax(x^2))
})

test_that("stencil matches stats::filter once the weight order is matched", {
  # sf_stencil uses the correlation convention, w[1] multiplying the earliest
  # element, so c(-1/2, 0, 1/2) is a centred first difference. stats::filter
  # uses the convolution convention and reverses the weights. For symmetric
  # weights the two agree; for anything else they differ by that reversal.
  set.seed(11)
  x <- cumsum(rnorm(1000))
  for (w in list(rep(1/3, 3), rep(1/5, 5), c(-0.5, 0, 0.5), c(1, -2, 1),
                 c(0.1, 0.2, 0.7))) {
    got <- sf_stencil(quote(x), w)
    want <- as.numeric(stats::filter(x, rev(w), sides = 2))
    expect_equal(got, want, tolerance = 1e-12,
                 info = paste(deparse(w), collapse = ""))
  }
  # the stencil convention is the one that makes a derivative read correctly
  d <- sf_stencil(quote(x), c(-0.5, 0, 0.5))
  expect_equal(d[5], (x[6] - x[4]) / 2)
})

test_that("stencil boundary rules behave", {
  x <- c(1, 2, 3, 4, 5)
  w <- c(1, 1, 1)
  expect_identical(is.na(sf_stencil(quote(x), w, "NA")), c(TRUE, rep(FALSE, 3), TRUE))
  expect_equal(sf_stencil(quote(x), w, "zero")[1], 0 + 1 + 2)
  expect_equal(sf_stencil(quote(x), w, "replicate")[1], 1 + 1 + 2)
  expect_equal(sf_stencil(quote(x), w, "reflect")[1], 1 + 1 + 2)
  expect_equal(sf_stencil(quote(x), w, "circular")[1], 5 + 1 + 2)
  expect_equal(sf_stencil(quote(x), w, "circular")[5], 4 + 5 + 1)
  expect_equal(sf_rollmean(quote(x), 3, "replicate")[3], mean(c(2, 3, 4)))
  expect_error(sf_stencil(quote(x), c(1, 1)), "odd length")
})

test_that("stencil works across task boundaries", {
  set.seed(12)
  n <- 3e5L                              # several 65536-element tasks
  x <- rnorm(n)
  got <- sf_stencil(quote(x), rep(1/5, 5), "replicate")
  want <- sf_stencil(quote(x), rep(1/5, 5), "replicate", threads = 1L)
  expect_identical(got, want)
  mid <- 150000L
  expect_equal(got[mid], mean(x[(mid - 2):(mid + 2)]))
})

test_that("segmented reduce matches a grouped apply", {
  set.seed(13)
  n <- 20000L
  g <- sort(sample.int(400L, n, replace = TRUE))
  y <- rnorm(n)
  res <- sf_segreduce(quote(y^2), g, "sum")
  want <- tapply(y^2, g, sum)
  expect_identical(res$id, as.integer(names(want)))
  expect_equal(res$value, as.numeric(want))
  expect_identical(res$length, as.integer(table(g)[as.character(res$id)]))

  resm <- sf_segreduce(quote(y), g, "mean")
  expect_equal(resm$value, as.numeric(tapply(y, g, mean)))
  resx <- sf_segreduce(quote(y), g, "max")
  expect_identical(resx$value, as.numeric(tapply(y, g, max)))
})

test_that("segments are runs, not groups, and NA starts a run", {
  y <- c(1, 2, 3, 4, 5, 6)
  g <- c(1, 1, 2, 2, 1, 1)               # id 1 appears twice, non-contiguously
  res <- sf_segreduce(quote(y), g, "sum")
  expect_identical(res$id, c(1L, 2L, 1L))
  expect_equal(res$value, c(3, 7, 11))

  y4 <- y[1:4]
  gna <- c(1, NA, NA, 2)
  resna <- sf_segreduce(quote(y4), gna, "sum")
  expect_identical(nrow(resna), 4L)
})

test_that("segmented reduce fuses with the expression", {
  y <- c(1, 2, 3, 10, 20); yh <- c(1, 1, 1, 2, 2); id <- c(1, 1, 1, 2, 2)
  res <- sf_segreduce(~ (y - yh)^2, id, "sum")
  expect_equal(res$value, c(sum((y[1:3] - 1)^2), sum((y[4:5] - 2)^2)))
})
