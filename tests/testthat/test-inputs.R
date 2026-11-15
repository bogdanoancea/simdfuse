
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

test_that("block aliasing is exact at every boundary", {
  # Full blocks point straight at the input and write straight to the output,
  # skipping two copies; only a block needing lane padding still copies. The
  # boundaries are therefore the interesting cases: lengths at, either side of,
  # and well past a block and a lane edge.
  old <- sf_options(); on.exit(sf_options(block = old$block, isa = old$isa))
  set.seed(7)
  for_each_isa(function(isa) {
    w <- sf_info()$lanes
    for (blk in c(16L, 128L)) {
      sf_options(block = blk)
      lens <- unique(c(1L, w - 1L, w, w + 1L, blk - 1L, blk, blk + 1L,
                       2L * blk, 2L * blk + w - 1L, 1001L))
      for (n in lens[lens >= 1L]) {
        x <- runif(n, .5, 3); y <- runif(n, .5, 3)
        if (n > 3) { x[2] <- NA; y[3] <- NaN; x[4] <- Inf }
        # a bare variable: the final register is a load, which must not alias
        expect_identical(sf_eval(x), x)
        expect_identical(sf_eval(x + y), x + y)
        expect_identical(sf_eval((x + y * 0.5 - x * 0.25) * (x - y) + x / y - 1),
                         (x + y * 0.5 - x * 0.25) * (x - y) + x / y - 1)
        # reductions and skeletons read the result from the register arena
        # rather than an output buffer, so they exercise the other path
        expect_identical(sf_sum(quote(x), na.rm = TRUE), sum(x, na.rm = TRUE))
        expect_equal(sf_scan(quote(x), op = "sum"), cumsum(x))
      }
    }
  })
  # inputs must come back untouched, since the kernel now reads them in place
  z <- runif(500, .5, 3); before <- z
  invisible(sf_eval(z * 2 + 1)); invisible(sf_sum(quote(z)))
  expect_identical(z, before)
})

test_that("the block cap keeps a program's register arena inside L1", {
  old <- sf_options()
  on.exit(sf_options(block = old$block, block_red = old$block_red,
                     l1_bytes = old$l1_bytes))

  # the cap holds the arena at a fixed share of L1 whatever the register count
  sf_options(l1_bytes = 49152L)                        # 48 KiB
  caps <- vapply(c(1L, 2L, 4L, 8L, 16L), simdfuse:::.sf_block_cap, 1L)
  arena <- caps * c(1L, 2L, 4L, 8L, 16L) * 8
  expect_true(all(arena <= 49152 / 3 + 64 * 8))        # within rounding
  expect_true(all(diff(caps) <= 0))                    # monotone in nreg
  expect_true(all(caps >= 64L))                        # never degenerate
  expect_true(all(caps %% 8L == 0L))                   # whole lanes

  # a bigger cache permits a bigger block
  sf_options(l1_bytes = 196608L)
  expect_gt(simdfuse:::.sf_block_cap(8L), caps[4])
  sf_options(l1_bytes = 49152L)

  set.seed(5)
  x <- runif(4000, .5, 3); mu <- mean(x); s <- sd(x)
  ll <- quote(-0.5 * log(2 * pi) - log(s) - 0.5 * ((x - mu) / s)^2)

  # the cap applies to the session default
  sf_options(block = 16384L, block_red = 16384L)
  expect_lt(.sf_block_for("sum", NULL, 8L), 16384L)
  expect_identical(.sf_block_for("sum", NULL, 8L), simdfuse:::.sf_block_cap(8L))

  # but an explicit block= is an instruction and is honoured as given
  expect_identical(.sf_block_for("sum", 16384L, 8L), 16384L)
  expect_identical(.sf_block_for("none", 4096L, 20L), 4096L)

  # a one-register program is not capped below the session default
  expect_identical(.sf_block_for("none", NULL, 1L),
                   min(16384L, simdfuse:::.sf_block_cap(1L)))

  # and the results do not change because of any of it
  for_each_isa(function(isa) {
    for (b in c(128L, 1024L, 16384L)) {
      sf_options(block = b, block_red = b)
      expect_equal(sf_sum(ll), sum(-0.5*log(2*pi) - log(s) - 0.5*((x-mu)/s)^2),
                   tolerance = 1e-12)
      expect_identical(sf_eval(ll), -0.5*log(2*pi) - log(s) - 0.5*((x-mu)/s)^2)
    }
  })
})
