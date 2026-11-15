# Block-size calibration.
#
# The best block is strongly machine-specific and the spread is large: on an
# Intel Xeon with 48 KiB of L1d the optimum sat near 128 doubles, while on an
# Apple M-series with 128 KiB it sat near 2048, and using the other machine's
# value cost 12 to 18 percent. There is no formula here that holds across
# microarchitectures, so the honest answer is to measure once per machine.

.sf_time <- function(f, min_time = 0.05, batches = 3L) {
  now <- function() as.numeric(Sys.time())
  f()
  reps <- 1L
  repeat {
    t0 <- now(); for (i in seq_len(reps)) f(); el <- now() - t0
    if (el >= min_time || reps >= 1e6) break
    growth <- if (el <= 0) 20 else max(2, min(50, min_time / el))
    reps <- as.integer(min(1e6, ceiling(reps * growth)))
  }
  ts <- numeric(batches)
  for (k in seq_len(batches)) {
    t0 <- now(); for (i in seq_len(reps)) f(); ts[k] <- (now() - t0) / reps
  }
  min(ts)
}

#' Calibrate the block size for this machine
#'
#' @description
#' Measures a representative mix of expressions across candidate block sizes
#' and reports, by default also sets, the one that is fastest here.
#'
#' @details
#' The block is how many elements are evaluated per pass, and therefore how
#' much of an expression's working set stays in L1 cache. The best value
#' depends on cache size and on microarchitecture in ways that no simple
#' formula captures: measured optima have ranged from 128 doubles on an Intel
#' Xeon to 2048 on an Apple M-series, with 12 to 18 percent lost by using the
#' wrong one. The package default is a compromise, so running this once per
#' machine is worth a few seconds.
#'
#' Four expressions are timed -- shallow arithmetic, deep arithmetic, a
#' transcendental pair, and a fused reduction -- because they do not agree.
#' Elementwise expressions often prefer large blocks, while a streaming
#' reduction degrades sharply once the block outgrows L1. The chosen block
#' minimises the geometric mean of the four times relative to the best seen for
#' each, so no single shape dominates. The per-case columns are returned so you
#' can override the compromise if your own work is all of one kind.
#'
#' `sf_tune()` changes an option, so it lasts for the session. To make it
#' permanent put the result in your `.Rprofile`:
#' `options(simdfuse.block = 2048L)`.
#'
#' @param blocks Candidate block sizes.
#' @param n Vector length to measure with. Should exceed the last-level cache
#'   for the result to be representative of large-data use, and ideally match
#'   the size you actually work at. This matters more than it looks: below that
#'   the array stays in cache, the block size barely registers, and the chosen
#'   value does not transfer upward. Tuning at 2e6 and running at 1e7 picked a
#'   reduction block that measurement at 1e7 rejected on every repeat.
#' @param set Set the tuned options for the session: both block sizes and,
#'   per [sf_options()], whether each elementary function uses its vectorised
#'   form.
#' @param quiet Suppress the printed summary.
#'
#' @return A data frame with one row per candidate block: its time for each
#'   expression in milliseconds, `relative`, the geometric mean over the
#'   elementwise cases relative to the best seen for each, and `relative_red`,
#'   the same for the reduction. Returned invisibly.
#' @export
#' @examples
#' \donttest{
#' sf_tune(n = 2e5, blocks = c(128L, 512L))   # small and quick, for illustration
#' }
sf_tune <- function(blocks = c(64L, 128L, 256L, 512L, 1024L, 2048L, 4096L),
                    n = 1e7, set = TRUE, quiet = FALSE) {
  blocks <- sort(unique(as.integer(blocks)))
  if (anyNA(blocks) || any(blocks < 1L))
    stop("simdfuse: `blocks` must be positive integers", call. = FALSE)
  n <- as.integer(n)
  if (is.na(n) || n < 1000L)
    stop("simdfuse: `n` must be at least 1000", call. = FALSE)

  a <- stats::runif(n, 0.5, 2)
  b <- stats::runif(n, 0.5, 2)
  # Four shapes, because they do not agree. Elementwise expressions tolerate
  # and often prefer large blocks; a streaming reduction, whose only work per
  # block is a copy and an accumulate, degrades sharply once the block stops
  # fitting comfortably in L1. Tuning on elementwise cases alone picks a block
  # that makes sum(x) measurably slower, which is exactly what happened on the
  # first machine this was tried on.
  cases <- list(
    shallow = list(quote(a + b * 0.5), "none"),
    deep    = list(quote((a + b * 0.5 - a * 0.25) * (a - b) + a / b - 1), "none"),
    math    = list(quote(exp(a) + log(b)), "none"),
    reduce  = list(quote((a - 1.25)^2), "sum")
  )
  env <- environment()

  tm <- matrix(NA_real_, nrow = length(blocks), ncol = length(cases),
               dimnames = list(NULL, names(cases)))
  for (j in seq_along(cases)) {
    ex <- cases[[j]][[1]]; red <- cases[[j]][[2]]
    for (i in seq_along(blocks)) {
      blk <- blocks[i]
      tm[i, j] <- .sf_time(function()
        .sf_run(ex, env, NULL, red, block = blk, threads = 1L))
    }
  }

  rel <- sweep(tm, 2L, apply(tm, 2L, min), "/")
  is_red <- vapply(cases, function(z) !identical(z[[2]], "none"), TRUE)
  # Two answers rather than one compromise: the elementwise cases and the
  # reduction disagree, sometimes by a factor of sixteen in the chosen block.
  gm <- exp(rowMeans(log(rel[, !is_red, drop = FALSE])))
  gm_red <- exp(rowMeans(log(rel[, is_red, drop = FALSE])))
  # These curves are flat near their minimum and the timings are noisy, so
  # which.min() alone picks an arbitrary point among ties and the choice jumps
  # between runs. Instead: take every block within TOL of the best for one
  # criterion, and among those choose the one that is best for the other. The
  # answer is stable across runs and the cost of splitting the two is smaller.
  TOL <- 1.02
  tied <- function(mine) which(mine <= min(mine) * TOL)
  pick <- function(mine, other) {
    ok <- tied(mine)
    ok2 <- ok[other[ok] <= min(other[ok]) * TOL]   # still tied on both
    # Largest of the sizes that measure the same, not smallest. Fewer blocks
    # means the per-block dispatch is amortised over more elements, and the
    # cost of getting this wrong is not symmetric: measured on a 2-lane
    # machine, a parallel reduction at block 64 ran about 10% slower than the
    # same reduction at 128, while the two were within 2% single-threaded. The
    # single-thread sweep that chooses the block cannot see that, so where
    # measurement cannot separate two sizes, take the one with less overhead.
    blocks[max(ok2)]
  }
  best     <- pick(gm, gm_red)
  best_red <- pick(gm_red, gm)
  n_tied     <- length(tied(gm))
  n_tied_red <- length(tied(gm_red))

  out <- data.frame(block = blocks, round(tm * 1000, 3),
                    relative = round(gm, 3), relative_red = round(gm_red, 3),
                    row.names = NULL)
  if (set) sf_options(block = best, block_red = best_red)
  if (!quiet) {
    cat("simdfuse block calibration, n = ", format(n, big.mark = ","),
        ", isa = ", sf_info()$isa, "\n", sep = "")
    print(out, row.names = FALSE)
    cat("\nbest block, elementwise: ", best,
        "\nbest block, reductions : ", best_red,
        if (set) "\n(both set for this session)" else "\n(not set)", "\n", sep = "")
    if (n_tied > 1L || n_tied_red > 1L)
      cat(sprintf(paste0("(sizes within %.0f%% of the best: %d elementwise, %d reductions",
                         " -- a flat range rather than a sharp optimum, so the exact",
                         " value matters less than staying inside it)\n"),
                  100 * (TOL - 1), n_tied, n_tied_red))
    cat("to keep them, add to your .Rprofile:\n  options(simdfuse.block = ",
        best, "L, simdfuse.block_red = ", best_red, "L)\n", sep = "")
    if (best != best_red) {
      cost <- 100 * (max(gm[blocks == best_red], gm_red[blocks == best]) - 1)
      if (cost > 100 * (TOL - 1))
        cat(sprintf("using one value for both would cost about %.0f%% on whichever side lost\n", cost))
    }
  }
  ## ---- fast-math: is the vectorised polynomial actually worth it? --------
  ## Deliberately measured rather than assumed. The polynomial has a fixed
  ## operation count per vector, so it only amortises once there are enough
  ## lanes to pay for it -- but how many is enough depends on the local libm.
  ## Measured: 0.54x at 1 lane and 1.70x at 8 lanes against glibc, yet 1.71x
  ## at 2 lanes against Apple's libm. Neither lane count nor platform predicts
  ## it alone, so each function is timed both ways on this machine.
  if (!quiet) cat("\nfast math (vectorised polynomial vs system libm)\n")
  old_f <- list(exp = getOption("simdfuse.fast_exp"),
                log = getOption("simdfuse.fast_log"))
  on.exit(options(simdfuse.fast_exp = old_f$exp,
                  simdfuse.fast_log = old_f$log), add = TRUE)
  # Measured over several input distributions, not one. A system libm branches
  # on its argument -- fast paths near 1, full range reduction elsewhere --
  # while the polynomial is branchless and costs the same everywhere, so the
  # ratio between them depends on the data. Measured spread on one machine was
  # a factor of 1.8 for log. Where the polynomial wins comfortably that changes
  # nothing; close to the crossover, which is where two-lane machines sit, it
  # flips the answer. So the decision is taken on the worst of the
  # distributions rather than on a single sample: fast mode trades accuracy for
  # speed, and a gain that depends on the data is a bad trade.
  grids <- list(near1 = c(0.9, 1.1), moderate = c(0.5, 2), wide = c(1e-6, 1e6))
  fm <- c(exp = NA_real_, log = NA_real_); fm_lo <- fm; fm_hi <- fm
  for (fn in c("exp", "log")) {
    opt <- paste0("simdfuse.fast_", fn)
    rs <- vapply(grids, function(g) {
      v <- stats::runif(n, g[1], g[2])
      if (fn == "exp") v <- log(v)          # keep exp() off the overflow path
      e2 <- new.env(parent = emptyenv()); assign("v", v, envir = e2)
      ex <- if (fn == "exp") quote(exp(v)) else quote(log(v))
      tt <- function(val) {
        do.call(options, stats::setNames(list(val), opt))
        .sf_time(function()
          .sf_run(ex, e2, NULL, "none", math = "fast", block = best, threads = 1L))
      }
      tt(FALSE) / tt(TRUE)
    }, 0)
    fm[[fn]] <- min(rs); fm_lo[[fn]] <- min(rs); fm_hi[[fn]] <- max(rs)
  }
  use <- fm > 1.05
  use[is.na(use)] <- TRUE
  if (!quiet) for (fn in names(fm))
    cat(sprintf("  %-4s %5.2fx to %5.2fx over %d input ranges -> %s\n", fn,
                fm_lo[[fn]], fm_hi[[fn]], length(grids),
                if (use[[fn]]) "vectorised" else "system libm"))

  if (set) {
    sf_options(fast_exp = use[["exp"]], fast_log = use[["log"]])
    on.exit()  # keep what was just set
  }
  if (!quiet && !all(use))
    cat("note: ", paste(names(use)[!use], collapse = ", "),
        " stay on the system libm here, so math = \"fast\" leaves ",
        if (sum(!use) == 1L) "that function" else "those functions",
        " bit-identical to base R\n", sep = "")

  invisible(list(table = out, block = best, block_red = best_red,
                 fastmath = use, fastmath_gain = fm,
                 fastmath_range = rbind(lo = fm_lo, hi = fm_hi)))
}
