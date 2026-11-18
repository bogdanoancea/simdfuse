# Shared benchmark engine for simdfuse.
#
# Sourced by two consumers:
#   * vignette("performance"), which runs it live on whatever machine builds
#     the vignette, so the numbers a reader sees are their own;
#   * bench-run.R and bench-block.R, which run the "full" profile and write
#     CSVs for the author's reference figures.
#
# Three things are deliberate here.
#
# Profiles. The full profile uses 10-20 million element vectors and is far too
# slow and too memory-hungry to run inside R CMD check. The vignette therefore
# runs "quick" (or "tiny" on CRAN), which keeps every measurement real but
# bounds the cost. The sizes still straddle the cache hierarchy, which is what
# the fusion result depends on; see the note in sf_bench_fusion.
#
# No timing dependency, and deliberately none: a vignette that dies because a
# benchmarking package is missing is a nuisance, and on a cluster or a build
# farm it is the likeliest thing to be missing. sf_bench_time calibrates its
# own repetition
# count against the clock and reports the minimum of several batches, which is
# the estimator least contaminated by interference from other processes.
#
# Every function returns a data frame with the same five columns, so the
# figures and tables do not care which profile produced them.

sf_bench_profile <- function(scale = c("quick", "tiny", "full")) {
  scale <- match.arg(scale)
  switch(scale,
    tiny = list(
      scale = "tiny", ns = c(1e4, 1e5), n = 2e5, n_thread = 4e5,
      min_time = 0.02, thread_cap = NA,
      blocks = c(32L, 64L, 128L, 256L, 512L, 2048L)),
    quick = list(
      scale = "quick", ns = c(1e4, 1e5, 1e6), n = 2e6, n_thread = 4e6,
      min_time = 0.05, thread_cap = NA,
      blocks = c(16L, 32L, 64L, 128L, 256L, 512L, 1024L, 4096L)),
    full = list(
      scale = "full", ns = c(1e4, 1e5, 1e6, 1e7), n = 1e7, n_thread = 2e7,
      min_time = 0.20, thread_cap = NA,
      blocks = c(16L, 32L, 64L, 128L, 192L, 256L, 384L, 512L, 768L,
                 1024L, 2048L, 4096L, 16384L))
  )
}

# Time a zero-argument function. Calibrates the repetition count until a batch
# takes at least min_time, then reports the fastest of `batches` such batches.
#
# Sys.time() rather than proc.time(): the latter has millisecond resolution, so
# a single call to a fast expression measures as exactly zero, and dividing
# min_time by that sends the repetition count to tens of millions. The growth
# factor is clamped for the same reason -- an unlucky first measurement must
# not be able to turn one measurement into a ten-minute one.
sf_bench_time <- function(f, min_time = 0.05, batches = 3L, max_reps = 1e6) {
  now <- function() as.numeric(Sys.time())
  f()                                   # warm up: first-call compilation, caches
  reps <- 1L
  repeat {
    t0 <- now()
    for (i in seq_len(reps)) f()
    el <- now() - t0
    if (el >= min_time || reps >= max_reps) break
    growth <- if (el <= 0) 20 else max(2, min(50, min_time / el))
    reps <- as.integer(min(max_reps, ceiling(reps * growth)))
  }
  ts <- numeric(batches)
  for (b in seq_len(batches)) {
    t0 <- now()
    for (i in seq_len(reps)) f()
    ts[b] <- (now() - t0) / reps
  }
  # The minimum is the estimate, but the spread across batches is carried along
  # so that a reader can tell a real difference from scheduler noise. On a
  # laptop with asymmetric cores, repeated threaded measurements of the same
  # expression have been observed to differ by 25 percent, which is larger
  # than several of the effects being compared.
  out <- min(ts)
  attr(out, "spread") <- if (min(ts) > 0) (max(ts) - min(ts)) / min(ts) else NA_real_
  out
}

.sf_row <- function(bench, n, case, engine, seconds) {
  sp <- attr(seconds, "spread")
  data.frame(bench = bench, n = n, case = case, engine = engine,
             seconds = as.numeric(seconds),
             spread = if (is.null(sp)) NA_real_ else as.numeric(sp),
             stringsAsFactors = FALSE)
}

# ---- the machine this ran on ----------------------------------------------

sf_bench_sysinfo <- function() {
  info <- simdfuse::sf_info()
  cache <- function() {
    out <- try(switch(Sys.info()[["sysname"]],
      Darwin = {
        v <- suppressWarnings(system2("sysctl",
               c("-n", "hw.l1dcachesize", "hw.l2cachesize"),
               stdout = TRUE, stderr = FALSE))
        if (length(v) >= 2) sprintf("L1d %s KiB, L2 %s KiB",
                                    round(as.numeric(v[1]) / 1024),
                                    round(as.numeric(v[2]) / 1024)) else NA
      },
      Linux = {
        p <- "/sys/devices/system/cpu/cpu0/cache"
        if (!dir.exists(p)) NA else {
          lv <- list.dirs(p, recursive = FALSE)
          paste(vapply(lv, function(d) {
            ty <- tryCatch(readLines(file.path(d, "type"), warn = FALSE), error = function(e) "")
            sz <- tryCatch(readLines(file.path(d, "size"), warn = FALSE), error = function(e) "")
            lvl <- tryCatch(readLines(file.path(d, "level"), warn = FALSE), error = function(e) "")
            if (length(ty) && ty %in% c("Data", "Unified")) sprintf("L%s %s", lvl, sz) else ""
          }, ""), collapse = " ")
        }
      }, NA), silent = TRUE)
    if (inherits(out, "try-error") || !length(out) || is.na(out[1])) NA_character_
    else trimws(gsub("\\s+", " ", out[1]))
  }
  list(
    platform    = R.version$platform,
    r_version   = paste(R.version$major, R.version$minor, sep = "."),
    isa         = info$isa,
    lanes       = info$lanes,
    available   = paste(info$available, collapse = ", "),
    openmp      = info$openmp,
    max_threads = info$max_threads,
    cores       = tryCatch(parallel::detectCores(logical = FALSE),
                           error = function(e) NA_integer_),
    logical     = tryCatch(parallel::detectCores(logical = TRUE),
                           error = function(e) NA_integer_),
    cache       = cache(),
    # Anything in the environment that can hold OpenMP below the core count.
    # R CMD check --as-cran and some BLAS setups set these, and without
    # reporting them a capped thread sweep looks like a hardware limit.
    omp_env     = {
      v <- Sys.getenv(c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT",
                        "OMP_DYNAMIC", "R_NUM_THREADS",
                        "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))
      v <- v[nzchar(v)]
      if (length(v)) paste(sprintf("%s=%s", names(v), v), collapse = ", ")
      else NA_character_
    },
    under_check = nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")) ||
                  nzchar(Sys.getenv("_R_CHECK_TIMINGS_"))
  )
}

# ---- 1. fusion -------------------------------------------------------------
# The quantity being removed is traversals of memory, so the benefit should
# grow with the number of operators and with n. The largest n in a profile
# needs to exceed the last-level cache for the effect to be fully visible; at
# the "quick" and "tiny" sizes part of the working set still fits, so the
# measured speedup is a lower bound on what the full profile shows.
sf_bench_fusion <- function(p) {
  exprs <- list(
    "1 op"   = quote(a + b),
    "2 ops"  = quote(a + b * 0.5),
    "4 ops"  = quote(a + b * 0.5 - a * 0.25),
    "8 ops"  = quote((a + b * 0.5 - a * 0.25) * (a - b) + a / b - 1),
    "16 ops" = quote(((a + b * 0.5 - a * 0.25) * (a - b) + a / b - 1) *
                     ((a - b * 0.5 + a * 0.25) * (a + b) - a / b + 1))
  )
  rows <- list()
  for (n in p$ns) {
    a <- runif(n, 0.5, 2); b <- runif(n, 0.5, 2)
    env <- environment()
    for (lab in names(exprs)) {
      ex <- exprs[[lab]]
      tb <- sf_bench_time(function() eval(ex, env), p$min_time)
      ts <- sf_bench_time(function() simdfuse::sf_eval(ex, threads = 1L), p$min_time)
      rows[[length(rows) + 1L]] <- rbind(
        .sf_row("fusion", n, lab, "base R", tb),
        .sf_row("fusion", n, lab, "simdfuse", ts))
    }
    rm(a, b)
  }
  do.call(rbind, rows)
}

# ---- 2. elementary functions ----------------------------------------------
sf_bench_mathfun <- function(p) {
  n <- p$n
  x <- runif(n, 0.1, 10); xs <- runif(n, -5, 5); m <- 3; sg <- 1.5
  env <- environment()
  cases <- list(
    "exp(x)"          = quote(exp(xs)),
    "log(x)"          = quote(log(x)),
    "sqrt(x)"         = quote(sqrt(x)),
    "gaussian loglik" = quote(-0.5 * log(2 * pi) - log(sg) - 0.5 * ((x - m) / sg)^2)
  )
  rows <- list()
  for (lab in names(cases)) {
    ex <- cases[[lab]]
    rows[[length(rows) + 1L]] <- rbind(
      .sf_row("mathfun", n, lab, "base R",
              sf_bench_time(function() eval(ex, env), p$min_time)),
      .sf_row("mathfun", n, lab, "simdfuse accurate",
              sf_bench_time(function() simdfuse::sf_eval(ex, math = "accurate",
                                                         threads = 1L), p$min_time)),
      .sf_row("mathfun", n, lab, "simdfuse fast",
              sf_bench_time(function() simdfuse::sf_eval(ex, math = "fast",
                                                         threads = 1L), p$min_time)))
  }
  do.call(rbind, rows)
}

# ---- 3. reductions ---------------------------------------------------------
# Floor for the two benchmarks that measure threading.
#
# Threading does not engage until the engine has tasks to hand out. A
# streaming expression needs SF_MIN_TASKS_PER_THREAD (4) tasks of 65536
# elements per worker, so a twelve-thread machine needs about 3.1M elements
# before a sweep measures anything at all; below that the engine correctly
# runs serially and every "all threads" bar is the serial path wearing a
# different label. The `tiny` and `quick` profiles sit far below that, which
# is how a vignette built for CRAN came to report threading as a slowdown.
#
# These two benchmarks therefore set their own n, independent of the size
# profile. They are single-pass and cheap, so the floor costs about a second
# even at `tiny`. Expressions that call libm per element parallelise from a
# much smaller n, but the floor is set for the streaming case because that is
# what these two measure.
SF_BENCH_PAR_FLOOR <- 4e6

sf_bench_reduce <- function(p) {
  n <- max(p$n, SF_BENCH_PAR_FLOOR)
  x <- runif(n, 0.1, 10); mu <- mean(x); m <- 3; sg <- 1.5
  env <- environment()
  cases <- list(
    "sum(x)"          = list(quote(sum(x)), quote(x), simdfuse::sf_sum),
    "sum((x-mu)^2)"   = list(quote(sum((x - mu)^2)), quote((x - mu)^2), simdfuse::sf_sum),
    "mean(abs(x-mu))" = list(quote(mean(abs(x - mu))), quote(abs(x - mu)), simdfuse::sf_mean),
    "max(x*x)"        = list(quote(max(x * x)), quote(x * x), simdfuse::sf_max),
    "gaussian loglik" = list(
      quote(sum(-0.5 * log(2 * pi) - log(sg) - 0.5 * ((x - m) / sg)^2)),
      quote(-0.5 * log(2 * pi) - log(sg) - 0.5 * ((x - m) / sg)^2), simdfuse::sf_sum)
  )
  rows <- list()
  for (lab in names(cases)) {
    cs <- cases[[lab]]; bex <- cs[[1]]; sex <- cs[[2]]; f <- cs[[3]]
    rr <- rbind(
      .sf_row("reduce", n, lab, "base R",
              sf_bench_time(function() eval(bex, env), p$min_time)),
      .sf_row("reduce", n, lab, "simdfuse 1 thread",
              sf_bench_time(function() f(sex, threads = 1L), p$min_time)),
      .sf_row("reduce", n, lab, "simdfuse all threads",
              sf_bench_time(function() f(sex), p$min_time, batches = 5L)))
    if (lab == "sum(x)" && requireNamespace("collapse", quietly = TRUE))
      rr <- rbind(rr, .sf_row("reduce", n, lab, "collapse::fsum",
                   sf_bench_time(function() collapse::fsum(x), p$min_time)))
    rows[[length(rows) + 1L]] <- rr
  }
  do.call(rbind, rows)
}

# ---- 4. instruction set sweep ---------------------------------------------
# The same work through each backend the CPU supports. This is the measurement
# that separates the SIMD contribution from the fusion contribution.
sf_bench_isa <- function(p) {
  old <- simdfuse::sf_options()
  on.exit(simdfuse::sf_options(isa = old$isa, fast_exp = old$fast_exp), add = TRUE)
  # The point of the isa_exp row is to show the vectorised polynomial running
  # through each backend, so it has to force the polynomial on. Left to the
  # session default it silently measures the system libm on every backend --
  # after sf_tune() has turned fast_exp off, every row comes out identical and
  # the panel looks like the dispatch layer does nothing.
  simdfuse::sf_options(fast_exp = TRUE)
  n <- p$n
  a <- runif(n, 0.5, 2); b <- runif(n, 0.5, 2); xs <- runif(n, -5, 5)
  env <- environment()
  ex8 <- quote((a + b * 0.5 - a * 0.25) * (a - b) + a / b - 1)
  exe <- quote(exp(xs))
  rows <- list()
  for (isa in simdfuse::sf_info()$available) {
    simdfuse::sf_options(isa = isa)
    rows[[length(rows) + 1L]] <- rbind(
      .sf_row("isa", n, isa, "simdfuse",
              sf_bench_time(function() simdfuse::sf_eval(ex8, threads = 1L), p$min_time)),
      .sf_row("isa_exp", n, isa, "simdfuse fast exp",
              sf_bench_time(function() simdfuse::sf_eval(exe, math = "fast",
                                                         threads = 1L), p$min_time)))
  }
  do.call(rbind, rows)
}

# ---- 5. block size --------------------------------------------------------
sf_bench_block <- function(p) {
  n <- p$n
  a <- runif(n, 0.5, 2); b <- runif(n, 0.5, 2)
  env <- environment()
  # The reduction is here because it disagrees with the rest: it prefers a
  # small block where elementwise expressions tolerate large ones.
  cases <- list(
    "2 ops"  = list(quote(a + b * 0.5), "accurate", FALSE),
    "8 ops"  = list(quote((a + b * 0.5 - a * 0.25) * (a - b) + a / b - 1), "accurate", FALSE),
    "16 ops" = list(quote(((a + b * 0.5 - a * 0.25) * (a - b) + a / b - 1) *
                          ((a - b * 0.5 + a * 0.25) * (a + b) - a / b + 1)), "accurate", FALSE),
    "fast exp" = list(quote(exp(a)), "fast", FALSE),
    "sum reduce" = list(quote((a - 1.25)^2), "accurate", TRUE)
  )
  rows <- list()
  for (lab in names(cases)) {
    ex <- cases[[lab]][[1]]; mm <- cases[[lab]][[2]]; isred <- cases[[lab]][[3]]
    for (blk in p$blocks)
      rows[[length(rows) + 1L]] <- .sf_row("block", n, lab, as.character(blk),
        sf_bench_time(function() if (isred)
          simdfuse::sf_sum(ex, block = blk, threads = 1L, math = mm)
        else
          simdfuse::sf_eval(ex, block = blk, threads = 1L, math = mm), p$min_time))
  }
  # engine holds the block size here, which keeps every panel one shape
  do.call(rbind, rows)
}

# ---- 6. thread scaling ----------------------------------------------------
sf_bench_threads <- function(p) {
  mx <- simdfuse::sf_info()$max_threads
  # Only an actual check run restricts the sweep, and then for CRAN's sake
  # rather than ours. Tying the cap to the size profile instead -- as an
  # earlier version did -- meant that R CMD build, which is how the vignette
  # is normally produced, reported two threads on every machine regardless of
  # its core count.
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")) ||
      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) mx <- min(mx, 2L)
  if (!is.na(p$thread_cap)) mx <- min(mx, p$thread_cap)
  mx <- max(1L, as.integer(mx))
  # A geometric ladder, plus the maximum itself: 1 2 4 8 12 on twelve cores
  # rather than all twelve. The count of measurements stays bounded while the
  # coverage follows whatever hardware this is running on -- an earlier version
  # capped this at two or four and so under-reported every real machine.
  ks <- sort(unique(c(2L^(0:ceiling(log2(mx))), mx)))
  ks <- ks[ks >= 1L & ks <= mx]
  n <- max(p$n_thread, SF_BENCH_PAR_FLOOR)
  a <- runif(n, 0.5, 2); b <- runif(n, 0.5, 2)
  env <- environment()
  ex <- quote((a + b * 0.5 - a * 0.25) * (a - b) + a / b - 1)
  rows <- list(.sf_row("threads", n, "base R", "base R",
                       sf_bench_time(function() eval(ex, env), p$min_time)))
  for (th in ks)
    rows[[length(rows) + 1L]] <- .sf_row("threads", n, as.character(th), "simdfuse",
      sf_bench_time(function() simdfuse::sf_eval(ex, threads = th), p$min_time,
                    batches = 5L))
  do.call(rbind, rows)
}

# ---- 7. skeletons ---------------------------------------------------------
sf_bench_skeletons <- function(p) {
  n <- p$n
  x <- runif(n, 0.1, 10); mu <- mean(x)
  env <- environment()
  w5 <- rep(1 / 5, 5)
  g <- sort(sample.int(max(2L, as.integer(n / 100)), n, replace = TRUE))
  rows <- list(
    .sf_row("scan", n, "cumsum", "base R",
            sf_bench_time(function() cumsum(x), p$min_time)),
    .sf_row("scan", n, "cumsum", "simdfuse serial",
            sf_bench_time(function() simdfuse::sf_scan(quote(x), "sum",
                                                       parallel = FALSE), p$min_time)),
    .sf_row("scan", n, "cumsum", "simdfuse parallel",
            sf_bench_time(function() simdfuse::sf_scan(quote(x), "sum"), p$min_time)),
    .sf_row("scan", n, "cumsum(log(x))", "base R",
            sf_bench_time(function() cumsum(log(x)), p$min_time)),
    .sf_row("scan", n, "cumsum(log(x))", "simdfuse",
            sf_bench_time(function() simdfuse::sf_scan(quote(log(x)), "sum"), p$min_time)),
    .sf_row("stencil", n, "5-point mean", "stats::filter",
            sf_bench_time(function() stats::filter(x, w5, sides = 2), p$min_time)),
    .sf_row("stencil", n, "5-point mean", "simdfuse",
            sf_bench_time(function() simdfuse::sf_stencil(quote(x), w5), p$min_time)),
    .sf_row("segreduce", n, "sum by group", "base rowsum",
            sf_bench_time(function() rowsum((x - mu)^2, g, reorder = FALSE), p$min_time)),
    .sf_row("segreduce", n, "sum by group", "simdfuse",
            sf_bench_time(function() simdfuse::sf_segreduce(quote((x - mu)^2), g,
                                                            "sum"), p$min_time)))
  if (requireNamespace("data.table", quietly = TRUE))
    rows[[length(rows) + 1L]] <- .sf_row("stencil", n, "5-point mean",
      "data.table::frollmean",
      sf_bench_time(function() data.table::frollmean(x, 5, align = "center"), p$min_time))
  if (requireNamespace("collapse", quietly = TRUE))
    rows[[length(rows) + 1L]] <- .sf_row("segreduce", n, "sum by group",
      "collapse::fsum",
      sf_bench_time(function() collapse::fsum((x - mu)^2, g), p$min_time))
  do.call(rbind, rows)
}

# ---- 8. peak memory -------------------------------------------------------
# Reported in MiB rather than seconds; the column name is kept so every panel
# shares one shape.
sf_bench_memory <- function(p) {
  n <- p$n
  x <- runif(n, 0.1, 10); mu <- mean(x)
  peak <- function(f) {
    gc(reset = TRUE, full = TRUE); invisible(f())
    unname(gc(full = TRUE)["Vcells", "max used"]) * 8 / 2^20
  }
  rbind(
    .sf_row("memory", n, "sum((x-mu)^2+sqrt(x))", "base R",
            peak(function() sum((x - mu)^2 + sqrt(x)))),
    .sf_row("memory", n, "sum((x-mu)^2+sqrt(x))", "simdfuse",
            peak(function() simdfuse::sf_sum(quote((x - mu)^2 + sqrt(x))))))
}

# ---- everything -----------------------------------------------------------
sf_bench_all <- function(p = sf_bench_profile("quick"), which = NULL,
                         verbose = FALSE) {
  panels <- list(
    fusion    = sf_bench_fusion,
    mathfun   = sf_bench_mathfun,
    reduce    = sf_bench_reduce,
    isa       = sf_bench_isa,
    block     = sf_bench_block,
    threads   = sf_bench_threads,
    skeletons = sf_bench_skeletons,
    memory    = sf_bench_memory
  )
  if (!is.null(which)) panels <- panels[which]
  out <- list()
  for (nm in names(panels)) {
    if (verbose) cat("  ", nm, " ... ", sep = "")
    t0 <- proc.time()[["elapsed"]]
    out[[nm]] <- panels[[nm]](p)
    if (verbose) cat(sprintf("%.1fs\n", proc.time()[["elapsed"]] - t0))
  }
  res <- do.call(rbind, out)
  attr(res, "profile") <- p$scale
  attr(res, "sysinfo") <- sf_bench_sysinfo()
  rownames(res) <- NULL
  res
}
