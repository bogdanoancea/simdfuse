.sf_scan_ops <- c(sum = 0L, prod = 1L, max = 2L, min = 3L, logsumexp = 4L)
.sf_boundaries <- c(NA_ = 0L, zero = 1L, replicate = 2L, reflect = 3L,
                    circular = 4L)

#' Parallel associative scan over an expression
#'
#' @description
#' A prefix scan under any of several associative operators, computed over the
#' *result of an expression* rather than over a stored vector. R has
#' `cumsum()`, `cumprod()`, `cummax()` and `cummin()`, all strictly serial and
#' all requiring their input to exist first; `sf_scan(~ log(x / y))` scans the
#' transformed series without ever building it, and spreads the work over
#' several cores.
#'
#' @details
#' The parallel algorithm is reduce-then-scan: one sweep computes a total per
#' task, a serial exclusive prefix over those totals fixes each task's starting
#' offset, and a second sweep writes the output. That costs a second evaluation
#' of the expression, so for a cheap expression on a short vector the serial
#' path wins -- hence `parallel = FALSE`, which is also the only mode that
#' reproduces `cumsum()` bit for bit (it carries the running total in long
#' double exactly as base R does).
#'
#' `op = "logsumexp"` accumulates in the log domain with the usual
#' max-shift, giving a numerically stable running normaliser. Base R has no
#' vectorised equivalent.
#'
#' @param expr An elementwise expression or one-sided formula.
#' @param op One of `"sum"`, `"prod"`, `"max"`, `"min"`, `"logsumexp"`.
#' @param data Optional list, data frame or environment of variables.
#' @param parallel `TRUE`, `FALSE`, or `NA` (the default) to let the engine
#'   decide from the vector length.
#' @param math,block,threads See [sf()].
#'
#' @return A double vector of inclusive prefix values.
#' @export
#' @examples
#' x <- c(3, 1, 4, 1, 5)
#' sf_scan(~ x, "sum")
#' identical(sf_scan(~ x, "sum", parallel = FALSE), cumsum(x))
#'
#' # a stable running log-normaliser
#' lp <- rnorm(10) * 50
#' sf_scan(~ lp, "logsumexp")
sf_scan <- function(expr, op = c("sum", "prod", "max", "min", "logsumexp"),
                    data = NULL, parallel = NA, math = NULL, block = NULL,
                    threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_scan_impl(r$ex, r$env, match.arg(op), data, parallel, math, block, threads)
}

.sf_scan_impl <- function(ex, env, op, data = NULL, parallel = NA, math = NULL,
                          block = NULL, threads = NULL) {
  r <- list(ex = ex, env = env)
  math <- match.arg(math %||% getOption("simdfuse.math"), c("accurate", "fast"))
  prog <- .sf_program_for(r$ex, identical(math, "fast"))
  inp <- .sf_fetch(prog, data, r$env)
  prog <- .sf_specialise(prog, inp, identical(math, "fast"))
  thr <- as.integer(threads %||% getOption("simdfuse.threads"))
  if (identical(parallel, FALSE)) thr <- 1L
  if (identical(parallel, TRUE) && thr == 1L) thr <- 0L
  out <- .Call(C_sf_scan, prog$code, prog$consts, inp$vars, as.double(inp$n),
               prog$nreg, prog$nprologue,
               .sf_block_for("none", block, prog$nreg),
               .sf_fast_mask(prog$fast), .sf_scan_ops[[op]], thr)
  .sf_restore_attrs(out, inp$proto)
}

#' @rdname sf_scan
#' @param ... Further arguments passed to [sf_scan()].
#' @export
sf_cumsum <- function(expr, data = NULL, ...) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_scan_impl(r$ex, r$env, "sum", data = data, ...)
}

#' Weighted stencil over an expression
#'
#' @description
#' Applies a fixed set of weights across a sliding window:
#' `out[i] = sum(w * expr[(i - h):(i + h)])` where `h` is half the window.
#' Note the orientation: `weights[1]` multiplies the *earliest* element of the
#' window, so `c(-1/2, 0, 1/2)` is a centred first difference. This is the
#' correlation convention; `stats::filter()` uses the convolution convention
#' and so takes the weights in the reverse order. For symmetric weights the
#' two agree.
#' This is a one-dimensional convolution, evaluated over the result of an
#' expression without materialising it, in parallel, with the window
#' arithmetic vectorised across output positions.
#'
#' @details
#' `stats::filter()` and `data.table::frollmean()` cover the common cases on a
#' stored vector; what is different here is that the input is an expression, so
#' smoothing a transformed series costs one pass rather than two plus a
#' temporary. Note that the inner product is evaluated with fused
#' multiply-add, so results can differ from `filter()` in the last bit.
#'
#' @param expr An elementwise expression or one-sided formula.
#' @param weights A numeric vector of odd length; the middle element is the
#'   weight on the current position.
#' @param boundary How to treat positions whose window runs off the end:
#'   `"NA"` (the default, matching `filter()`), `"zero"`, `"replicate"`,
#'   `"reflect"` or `"circular"`.
#' @param data Optional list, data frame or environment of variables.
#' @param math,block,threads See [sf()].
#'
#' @return A double vector the same length as the input.
#' @export
#' @examples
#' x <- cumsum(rnorm(100))
#' sm <- sf_stencil(~ x, rep(1/5, 5))
#' head(sm, 6)
#'
#' # a first difference of a log transform, smoothed, in one pass
#' sf_stencil(~ log(abs(x) + 1), c(-0.5, 0, 0.5), boundary = "replicate")[1:5]
sf_stencil <- function(expr, weights,
                       boundary = c("NA", "zero", "replicate", "reflect",
                                    "circular"),
                       data = NULL, math = NULL, block = NULL,
                       threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_stencil_impl(r$ex, r$env, weights, match.arg(boundary), data, math,
                   block, threads)
}

.sf_stencil_impl <- function(ex, env, weights, boundary, data = NULL,
                             math = NULL, block = NULL, threads = NULL) {
  r <- list(ex = ex, env = env)
  weights <- as.double(weights)
  if (length(weights) %% 2L != 1L)
    stop("simdfuse: `weights` must have odd length so the window has a centre",
         call. = FALSE)
  if (anyNA(weights)) stop("simdfuse: `weights` must not contain NA", call. = FALSE)
  math <- match.arg(math %||% getOption("simdfuse.math"), c("accurate", "fast"))
  prog <- .sf_program_for(r$ex, identical(math, "fast"))
  inp <- .sf_fetch(prog, data, r$env)
  prog <- .sf_specialise(prog, inp, identical(math, "fast"))
  bnd <- .sf_boundaries[[if (boundary == "NA") "NA_" else boundary]]
  out <- .Call(C_sf_stencil, prog$code, prog$consts, inp$vars,
               as.double(inp$n), prog$nreg, prog$nprologue,
               .sf_block_for("none", block, prog$nreg),
               .sf_fast_mask(prog$fast), weights, bnd,
               as.integer(threads %||% getOption("simdfuse.threads")))
  .sf_restore_attrs(out, inp$proto)
}

#' @rdname sf_stencil
#' @param k Window width; must be odd.
#' @param ... Further arguments passed to [sf_stencil()].
#' @export
sf_rollmean <- function(expr, k, boundary = c("NA", "zero", "replicate",
                                              "reflect", "circular"),
                        data = NULL, ...) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  k <- as.integer(k)
  if (is.na(k) || k < 1L || k %% 2L != 1L)
    stop("simdfuse: `k` must be a positive odd integer", call. = FALSE)
  .sf_stencil_impl(r$ex, r$env, rep(1 / k, k), match.arg(boundary),
                   data = data, ...)
}

#' Reduce an expression within contiguous segments
#'
#' @description
#' Collapses an expression within runs of equal values in `by`, returning one
#' value per run. This is the segmented-reduction primitive from data-parallel
#' computing, and it fits the shape of a sorted panel: with `by` the unit
#' identifier, `sf_segreduce(~ (y - yhat)^2, id, "sum")` gives the
#' within-unit residual sum of squares in a single pass over the data, with no
#' intermediate vector and no grouping table.
#'
#' @details
#' Segments must be *contiguous*: a new run starts whenever `by[i]` differs
#' from `by[i - 1]`, and every `NA` starts a run of its own. Sort by the
#' grouping variable first if it is not already grouped. This is deliberate --
#' it is what keeps the operation single-pass and lets runs be processed
#' independently, which also makes the result deterministic.
#'
#' @param expr An elementwise expression or one-sided formula.
#' @param by A vector of segment identifiers, the same length as the
#'   expression. Coerced to integer; factors use their codes.
#' @param reduce One of `"sum"`, `"mean"`, `"max"`, `"min"`, `"prod"`,
#'   `"any"`, `"all"`, `"nobs"`.
#' @param data Optional list, data frame or environment of variables.
#' @param na.rm Drop `NA` and `NaN` within each segment.
#' @param math,block,threads See [sf()].
#'
#' @return A data frame with columns `id`, `value` and `length`, one row per
#'   run, in order of first appearance.
#' @export
#' @examples
#' y    <- c(1, 2, 3, 10, 20, 100)
#' firm <- c(1, 1, 1,  2,  2,   3)
#' sf_segreduce(~ y^2, firm, "sum")
sf_segreduce <- function(expr, by,
                         reduce = c("sum", "mean", "max", "min", "prod",
                                    "any", "all", "nobs"),
                         data = NULL, na.rm = FALSE, math = NULL,
                         block = NULL, threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  reduce <- match.arg(reduce)
  math <- match.arg(math %||% getOption("simdfuse.math"), c("accurate", "fast"))
  prog <- .sf_program_for(r$ex, identical(math, "fast"))
  inp <- .sf_fetch(prog, data, r$env)
  prog <- .sf_specialise(prog, inp, identical(math, "fast"))
  if (is.factor(by)) by <- as.integer(by) else by <- as.integer(by)
  if (length(by) != inp$n)
    stop(sprintf("simdfuse: `by` has length %d but the expression has length %d",
                 length(by), inp$n), call. = FALSE)
  res <- .Call(C_sf_segreduce, prog$code, prog$consts, inp$vars,
               as.double(inp$n), prog$nreg, prog$nprologue,
               .sf_block_for("none", block, prog$nreg),
               .sf_fast_mask(prog$fast), by, .sf_redcode(reduce), isTRUE(na.rm),
               as.integer(threads %||% getOption("simdfuse.threads")))
  out <- data.frame(id = res$id, value = res$value, length = res$length)
  if (reduce %in% c("any", "all")) out$value <- as.logical(out$value)
  out
}
