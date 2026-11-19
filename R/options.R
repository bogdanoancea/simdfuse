#' Engine configuration
#'
#' @description
#' `sf_options()` reads or sets the engine's defaults. Called with no
#' arguments it returns the current settings.
#'
#' @param block Elements per evaluation block. The whole working set of an
#'   expression -- inputs, temporaries and result -- is this many doubles per
#'   register, so small values keep everything in L1 cache at the cost of more
#'   loop overhead. Rounded up to a multiple of the SIMD width. The default of
#'   128 (1 KiB per register) is a compromise: measured optima have ranged from
#'   64 doubles on an Intel Xeon to 2048 on an Apple M-series, and using the
#'   wrong machine's value costs 12 to 18 percent. Run [sf_tune()] once to find
#'   the best value here. The penalty is not monotone in the block size -- some
#'   sizes alias badly in a set-associative L1 cache -- so measuring beats
#'   reasoning about it.
#' @param block_red Block size for fused reductions, which prefer a smaller
#'   block than elementwise expressions because they stream without writing an
#'   output vector. `NA` (the default) means use `block` for both. [sf_tune()]
#'   measures and sets the two separately.
#' @param threads Worker threads. `0` means every thread OpenMP reports.
#'   Fewer than requested are used when the data are too small to keep them
#'   busy: a worker is created only once it would receive several blocks of
#'   work, because below that the cost of starting the team exceeds what the
#'   extra thread saves. Reductions return identical results at any setting,
#'   requested or actual.
#' @param isa Instruction set: `"auto"` (widest the CPU supports),
#'   or one of `"avx512"`, `"avx2"`, `"sse2"`, `"neon"`, `"scalar"`. A request
#'   the hardware cannot satisfy falls back to `"auto"` with a warning. Pinning
#'   a narrower backend is mainly useful for testing and for reproducing
#'   results across heterogeneous machines.
#' @param math `"accurate"` (default) routes every elementary function through
#'   the platform `libm` one lane at a time, reproducing base R bit for bit.
#'   The one exception is `sqrt()`, which every backend takes from the hardware:
#'   IEEE 754 requires that to be correctly rounded, while base R calls the C
#'   library, and on R 4.5.2 for Windows the library result is one ulp away from
#'   the true root for roughly one argument in two thousand. simdfuse returns
#'   the correctly rounded value there and so differs from base R.
#'   `"fast"` uses vectorised polynomial `exp()` and `log()` and allows
#'   fused multiply-add contraction, which is faster but not bit-identical;
#'   see `vignette("accuracy", package = "simdfuse")` for measured error.
#' @param fast_exp,fast_log Whether `math = "fast"` uses the vectorised
#'   polynomial for that function, rather than the platform `libm`. Whether the
#'   polynomial wins depends on the lane width and on the quality of the local
#'   `libm`, and the two do not point the same way. At four lanes or more it has
#'   won everywhere measured, by 1.5x to 13x; at two lanes it has mostly lost,
#'   and against Apple's `libm` on NEON the same `exp()` measured 1.71x on one
#'   M-series machine and 0.76x on an M1 Max. The default is therefore `TRUE`
#'   on AVX2 and AVX-512 and `FALSE` on narrower backends, and [sf_tune()]
#'   measures both on the machine in hand and sets them.
#' @param cache_max Number of compiled programs to keep.
#'
#' @return A list of the settings in effect, invisibly when setting.
#' @export
sf_options <- function(block = NULL, block_red = NULL, threads = NULL,
                       isa = NULL, math = NULL, fast_exp = NULL,
                       fast_log = NULL, cache_max = NULL) {
  changed <- FALSE
  if (!is.null(block)) {
    block <- as.integer(block)
    if (is.na(block) || block < 1L) stop("simdfuse: `block` must be a positive integer")
    options(simdfuse.block = block); changed <- TRUE
  }
  if (!is.null(block_red)) {
    block_red <- if (length(block_red) == 1L && is.na(block_red)) NA_integer_
                 else as.integer(block_red)
    if (!is.na(block_red) && block_red < 1L)
      stop("simdfuse: `block_red` must be a positive integer or NA")
    options(simdfuse.block_red = block_red); changed <- TRUE
  }
  if (!is.null(threads)) {
    threads <- as.integer(threads)
    if (is.na(threads) || threads < 0L) stop("simdfuse: `threads` must be >= 0")
    options(simdfuse.threads = threads); changed <- TRUE
  }
  if (!is.null(isa)) {
    isa <- match.arg(isa, c("auto", "avx512", "avx2", "sse2", "neon", "scalar"))
    # Bookkeeping runs on exit because the engine warns when a request cannot
    # be met, and a caller's warning handler unwinds this frame. It also
    # records what actually happened rather than what was asked for: a request
    # the hardware could not satisfy leaves the engine on its automatic choice.
    on.exit({
      .sf$info <- .Call(C_sf_info)
      options(simdfuse.isa = if (isTRUE(.sf$info$forced)) isa else "auto")
    }, add = TRUE)
    .Call(C_sf_set_isa, isa)
    changed <- TRUE
  }
  if (!is.null(math)) {
    math <- match.arg(math, c("accurate", "fast"))
    options(simdfuse.math = math); changed <- TRUE
  }
  for (nm in c("fast_exp", "fast_log")) {
    v <- get(nm)
    if (!is.null(v)) {
      v <- isTRUE(as.logical(v)[1])
      do.call(options, stats::setNames(list(v), paste0("simdfuse.", nm)))
      changed <- TRUE
    }
  }
  if (!is.null(cache_max)) {
    cache_max <- as.integer(cache_max)
    if (is.na(cache_max) || cache_max < 0L) stop("simdfuse: `cache_max` must be >= 0")
    options(simdfuse.cache_max = cache_max); changed <- TRUE
  }
  cur <- list(
    block     = getOption("simdfuse.block"),
    block_red = getOption("simdfuse.block_red"),
    threads   = getOption("simdfuse.threads"),
    isa       = getOption("simdfuse.isa"),
    math      = getOption("simdfuse.math"),
    fast_exp  = getOption("simdfuse.fast_exp"),
    fast_log  = getOption("simdfuse.fast_log"),
    cache_max = getOption("simdfuse.cache_max")
  )
  if (changed) invisible(cur) else cur
}

#' Report the active engine
#'
#' @description
#' Describes what the compiled engine resolved to on this machine: the
#' instruction set in use, its lane count, which backends the CPU supports,
#' and whether OpenMP is available.
#'
#' @return A list with elements `isa`, `lanes`, `forced`, `available`,
#'   `openmp`, `max_threads`, `task_elems` and `instr_width`.
#' @export
#' @examples
#' sf_info()
sf_info <- function() .Call(C_sf_info)

#' Functions the expression compiler understands
#'
#' @description
#' Any expression built only from these functions, numeric literals and
#' double, integer or logical vectors can be compiled. Anything else is
#' rejected at compile time with a message naming the offending call, so an
#' unsupported expression fails loudly rather than falling back to slow code.
#'
#' @return A data frame with one row per supported function.
#' @export
#' @examples
#' head(sf_supported(), 20)
sf_supported <- function() {
  tab <- .sf_fun_table()
  data.frame(
    fun   = names(tab),
    arity = vapply(tab, function(z) paste(z$arity, collapse = " or "), ""),
    exact = vapply(tab, function(z) !isTRUE(z$inexact), TRUE),
    row.names = NULL, stringsAsFactors = FALSE
  )
}
