`%||%` <- function(a, b) if (is.null(a)) b else a

# Accept either a bare expression, a one-sided formula written inline, or a
# formula held in a variable. A formula carries its own environment, which is
# the right place to look up its variables.
.sf_rhs <- function(ex, env) {
  if (is.call(ex) && identical(ex[[1L]], quote(`~`)))
    return(list(ex = ex[[length(ex)]], env = env))
  # quote(...) is unwrapped so the functions can be driven programmatically
  if (is.call(ex) && identical(ex[[1L]], quote(quote)) && length(ex) == 2L)
    return(list(ex = ex[[2L]], env = env))
  if (is.name(ex)) {
    val <- tryCatch(get(as.character(ex), envir = env), error = function(e) NULL)
    if (inherits(val, "formula"))
      return(list(ex = val[[length(val)]], env = environment(val) %||% env))
    # A variable holding a quoted expression is used as the expression itself,
    # which is what makes the functions usable programmatically as well as
    # interactively. A numeric variable is of course still just data.
    if (is.call(val) || (is.name(val) && !identical(val, ex)))
      return(list(ex = val, env = env))
  }
  list(ex = ex, env = env)
}

.sf_program_for <- function(ex, fast, scalars = character(0)) {
  key <- paste0(if (fast) "F:" else "A:",
                paste(sort(scalars), collapse = ","), ":",
                paste(deparse(ex, width.cutoff = 500L), collapse = "\n"))
  hit <- .sf$cache[[key]]
  if (!is.null(hit)) return(hit)
  prog <- .sf_compile_ast(ex, fast = fast, scalars = scalars)
  cm <- getOption("simdfuse.cache_max", 256L)
  if (cm > 0L) {
    assign(key, prog, envir = .sf$cache)
    .sf$cache_order <- c(.sf$cache_order, key)
    if (length(.sf$cache_order) > cm) {
      drop <- .sf$cache_order[1L]
      .sf$cache_order <- .sf$cache_order[-1L]
      if (exists(drop, envir = .sf$cache, inherits = FALSE))
        rm(list = drop, envir = .sf$cache)
    }
  }
  prog
}

# Coerce, check lengths, and pick the input whose attributes the elementwise
# result should inherit. simdfuse deliberately does not recycle: silent
# recycling of a mismatched length is a bug amplifier, not a convenience.
.sf_prepare <- function(prog, vals) {
  nv <- length(prog$varnames)
  vars <- vector("list", nv)
  n <- NA_integer_
  proto <- NULL
  protonm <- NULL
  for (i in seq_len(nv)) {
    nm <- prog$varnames[i]
    v <- vals[[i]]
    if (is.null(v))
      stop(sprintf("simdfuse: `%s` not found", nm), call. = FALSE)
    if (is.factor(v))
      stop(sprintf("simdfuse: `%s` is a factor; convert it explicitly", nm),
           call. = FALSE)
    # Integer and logical vectors are passed through in their own storage.
    # The block loader converts them as it copies into a register, so there is
    # no second vector the size of the input and no extra pass over it.
    if (!is.double(v) && !is.integer(v) && !is.logical(v))
      stop(sprintf("simdfuse: `%s` is of class %s; only double, integer and logical vectors are supported",
                   nm, paste(class(v), collapse = "/")), call. = FALSE)
    len <- length(v)
    if (len == 0L)
      stop(sprintf("simdfuse: `%s` has length 0", nm), call. = FALSE)
    if (len != 1L) {
      if (is.na(n)) {
        n <- len; proto <- v; protonm <- nm
      } else if (len != n) {
        stop(sprintf("simdfuse: `%s` has length %d but `%s` has length %d; simdfuse does not recycle",
                     nm, len, protonm, n), call. = FALSE)
      }
    }
    vars[[i]] <- v
  }
  if (is.na(n)) { n <- 1L; proto <- vars[[1L]] }
  list(vars = vars, n = n, proto = proto)
}

.sf_fetch <- function(prog, data, env) {
  vals <- vector("list", length(prog$varnames))
  for (i in seq_along(prog$varnames)) {
    nm <- prog$varnames[i]
    v <- NULL
    if (!is.null(data)) v <- data[[nm]]
    if (is.null(v)) {
      if (!exists(nm, envir = env))
        stop(sprintf("simdfuse: `%s` not found", nm), call. = FALSE)
      v <- get(nm, envir = env)
    }
    vals[[i]] <- v
  }
  .sf_prepare(prog, vals)
}

.sf_restore_attrs <- function(out, proto) {
  if (is.null(proto) || length(out) != length(proto)) return(out)
  d <- attr(proto, "dim")
  if (!is.null(d)) {
    dim(out) <- d
    dn <- attr(proto, "dimnames")
    if (!is.null(dn)) dimnames(out) <- dn
  } else {
    nmv <- attr(proto, "names")
    if (!is.null(nmv)) names(out) <- nmv
  }
  out
}

# Reductions and elementwise expressions want different block sizes, and by
# enough to matter: on an Apple M-series the elementwise optimum was 2048 while
# sum(x) wanted 128, and forcing one compromise value cost about 15 percent on
# whichever side lost. So a reduction uses simdfuse.block_red when that is set.
.sf_block_for <- function(red, block) {
  if (!is.null(block)) return(as.integer(block))
  if (!identical(red, "none")) {
    br <- getOption("simdfuse.block_red")
    if (!is.null(br) && !is.na(br)) return(as.integer(br))
  }
  as.integer(getOption("simdfuse.block"))
}

# Which elementary functions actually use their vectorised form. Only consulted
# when the program was compiled in fast mode; see sf_options(fast_exp=).
.sf_fast_mask <- function(fast) {
  if (!isTRUE(fast)) return(0L)
  m <- 0L
  if (!identical(getOption("simdfuse.fast_exp"), FALSE)) m <- m + 1L
  if (!identical(getOption("simdfuse.fast_log"), FALSE)) m <- m + 2L
  m
}

.sf_dispatch <- function(prog, inp, red, na.rm, block, threads) {
  blk <- .sf_block_for(red, block)
  thr <- as.integer(threads %||% getOption("simdfuse.threads"))
  out <- .Call(C_sf_eval, prog$code, prog$consts, inp$vars,
               as.double(inp$n), prog$nreg, prog$nprologue, blk,
               .sf_fast_mask(prog$fast),
               .sf_redcode(red), isTRUE(na.rm), thr)
  if (identical(red, "none")) {
    # A comparison or logical expression should come back logical, as it would
    # from base R; the engine carries logicals as 0/1 doubles with NA as NaN.
    if (isTRUE(prog$logical)) {
      att <- attributes(out)
      out <- as.logical(out)
      attributes(out) <- att
    }
    out <- .sf_restore_attrs(out, inp$proto)
  }
  out
}

.sf_run <- function(ex, env, data, red, na.rm = FALSE,
                    math = NULL, block = NULL, threads = NULL) {
  math <- match.arg(math %||% getOption("simdfuse.math"), c("accurate", "fast"))
  prog <- .sf_program_for(ex, identical(math, "fast"))
  inp <- .sf_fetch(prog, data, env)
  .sf_dispatch(.sf_specialise(prog, inp, identical(math, "fast")), inp,
               red, na.rm, block, threads)
}

# Which inputs are length one is only known once they have been fetched, so the
# program is compiled generically first and then respecialised with that
# knowledge. The specialised form is cached, so the second compile happens only
# the first time a given call site is seen with a given shape.
.sf_specialise <- function(prog, inp, fast) {
  sc <- prog$varnames[vapply(inp$vars, length, 1L) == 1L]
  if (length(sc) == 0L) return(prog)
  out <- .sf_program_for(prog$expr, fast, sc)
  if (!identical(out$varnames, prog$varnames))
    stop("simdfuse: internal error, input order changed on respecialisation")
  out
}

#' Evaluate an array expression in one fused pass
#'
#' @description
#' `sf()` evaluates an elementwise R expression without materialising the
#' intermediate vectors that ordinary R evaluation creates. `a * b + sqrt(c)`
#' in base R allocates three temporaries of length `n` and traverses memory
#' four times; `sf(a * b + sqrt(c))` compiles the expression once and walks the
#' data a single time, keeping every temporary inside a cache-resident block.
#'
#' The saving grows with the number of operations and with `n`: once the
#' operands no longer fit in cache, ordinary evaluation is limited by memory
#' bandwidth, and removing traversals is the only thing that helps.
#'
#' @param expr An elementwise expression, or a one-sided formula. Variables are
#'   looked up in `data` first, then in the calling frame. Only the functions
#'   listed by [sf_supported()] are allowed; anything else is an error rather
#'   than a silent fallback.
#' @param data Optional list, data frame or environment holding the variables.
#' @param math `"accurate"` (bit-identical to base R, except for `sqrt()` on a
#'   platform whose C library square root is not correctly rounded -- see
#'   [sf_options()]) or `"fast"`. Defaults to [sf_options()]`$math`.
#' @param block Elements per block; defaults to [sf_options()]`$block`.
#' @param threads Worker threads; defaults to [sf_options()]`$threads`.
#'
#' @return A double vector the length of the non-scalar inputs. `dim`,
#'   `dimnames` and `names` are carried over from the first non-scalar input,
#'   so a matrix in gives a matrix out.
#'
#' @section Length handling:
#' Inputs must all have the same length, or length 1. simdfuse does not recycle
#' partial lengths: a mismatch is an error, not a warning.
#'
#' @seealso [sf_sum()] and friends for fused reductions, [sf_scan()],
#'   [sf_stencil()] and [sf_segreduce()] for the parallel skeletons.
#' @export
#' @examples
#' a <- runif(1000); b <- runif(1000); c <- runif(1000)
#' all.equal(sf(a * b + sqrt(c) - 2.5), a * b + sqrt(c) - 2.5)
#'
#' # standardising without temporaries
#' mu <- mean(a); s <- sd(a)
#' z <- sf((a - mu) / s)
sf <- function(expr, data = NULL, math = NULL, block = NULL, threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_run(r$ex, r$env, data, "none", math = math, block = block,
          threads = threads)
}

#' @rdname sf
#' @export
sf_eval <- sf

#' Fused reductions
#'
#' @description
#' These collapse an expression to a scalar without ever allocating a vector of
#' results: `sf_sum(~ (x - mu)^2)` runs in constant auxiliary memory whatever
#' the length of `x`, where `sum((x - mu)^2)` allocates two full-length
#' temporaries.
#'
#' @section Reproducibility:
#' Parallel reductions here are deterministic in the thread count. Work is cut
#' into fixed-size tasks whose boundaries depend only on `n`, and task partials
#' are combined serially in index order in extended precision. The same call on
#' the same data therefore returns the same bits on one thread and on sixty,
#' which a plain OpenMP `reduction(+:)` clause cannot guarantee. The value is
#' not necessarily bit-identical to base R's `sum()`, which accumulates
#' sequentially in long double; the blocked summation used here is usually the
#' more accurate of the two. See `vignette("accuracy", package = "simdfuse")`.
#'
#' @param expr An elementwise expression or one-sided formula.
#' @param data Optional list, data frame or environment of variables.
#' @param na.rm Drop `NA` and `NaN` before reducing.
#' @param math,block,threads See [sf()].
#'
#' @return A length-one double, or logical for `sf_any()` and `sf_all()`.
#' @export
#' @examples
#' x <- rnorm(10000); mu <- mean(x)
#' sf_sum(~ (x - mu)^2)
#' sf_mean(~ abs(x))
#' sf_max(~ x^2)
#'
#' # a Gaussian log-likelihood with no intermediate vectors at all
#' s <- sd(x)
#' sf_sum(~ -0.5 * log(2 * pi) - log(s) - 0.5 * ((x - mu) / s)^2)
sf_sum <- function(expr, data = NULL, na.rm = FALSE, math = NULL,
                   block = NULL, threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_run(r$ex, r$env, data, "sum", na.rm, math, block, threads)
}

#' @rdname sf_sum
#' @export
sf_mean <- function(expr, data = NULL, na.rm = FALSE, math = NULL,
                    block = NULL, threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_run(r$ex, r$env, data, "mean", na.rm, math, block, threads)
}

#' @rdname sf_sum
#' @export
sf_max <- function(expr, data = NULL, na.rm = FALSE, math = NULL,
                   block = NULL, threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_run(r$ex, r$env, data, "max", na.rm, math, block, threads)
}

#' @rdname sf_sum
#' @export
sf_min <- function(expr, data = NULL, na.rm = FALSE, math = NULL,
                   block = NULL, threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_run(r$ex, r$env, data, "min", na.rm, math, block, threads)
}

#' @rdname sf_sum
#' @export
sf_prod <- function(expr, data = NULL, na.rm = FALSE, math = NULL,
                    block = NULL, threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_run(r$ex, r$env, data, "prod", na.rm, math, block, threads)
}

#' @rdname sf_sum
#' @export
sf_any <- function(expr, data = NULL, na.rm = FALSE, math = NULL,
                   block = NULL, threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_run(r$ex, r$env, data, "any", na.rm, math, block, threads)
}

#' @rdname sf_sum
#' @export
sf_all <- function(expr, data = NULL, na.rm = FALSE, math = NULL,
                   block = NULL, threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_run(r$ex, r$env, data, "all", na.rm, math, block, threads)
}

#' @rdname sf_sum
#' @export
sf_nobs <- function(expr, data = NULL, math = NULL, block = NULL,
                    threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  .sf_run(r$ex, r$env, data, "nobs", FALSE, math, block, threads)
}

#' Compile an expression into a reusable function
#'
#' @description
#' Returns a closure whose formal arguments are the expression's variables, in
#' the order they first appear. The expression is compiled once, so this is the
#' form to use inside an optimiser or a bootstrap loop where the same
#' expression is evaluated many times on different data.
#'
#' @param expr An elementwise expression or one-sided formula.
#' @param reduce One of `"none"`, `"sum"`, `"mean"`, `"max"`, `"min"`,
#'   `"prod"`, `"any"`, `"all"`, `"nobs"`.
#' @param na.rm Passed to the reduction.
#' @param math,block,threads See [sf()].
#'
#' @return A function of the expression's variables.
#' @export
#' @examples
#' nll <- sf_kernel(~ 0.5 * ((x - mu) / s)^2 + log(s), reduce = "sum")
#' names(formals(nll))
#' x <- rnorm(1000)
#' nll(x = x, mu = 0, s = 1)
sf_kernel <- function(expr,
                      reduce = c("none", "sum", "mean", "max", "min", "prod",
                                 "any", "all", "nobs"),
                      na.rm = FALSE, math = NULL, block = NULL,
                      threads = NULL) {
  r <- .sf_rhs(substitute(expr), parent.frame())
  reduce <- match.arg(reduce)
  math <- match.arg(math %||% getOption("simdfuse.math"), c("accurate", "fast"))
  prog <- .sf_compile_ast(r$ex, fast = identical(math, "fast"))

  shell <- function() NULL
  fargs <- rep(list(quote(expr = )), length(prog$varnames))
  names(fargs) <- prog$varnames
  formals(shell) <- fargs
  body(shell) <- quote({
    .vals <- mget(.sf_prog$varnames, envir = environment(), inherits = FALSE)
    .inp <- .sf_prepare(.sf_prog, .vals)
    .sf_dispatch(.sf_specialise(.sf_prog, .inp, .sf_prog$fast), .inp,
                 .sf_red, .sf_narm, .sf_block, .sf_threads)
  })
  e <- new.env(parent = asNamespace("simdfuse"))
  assign(".sf_prog", prog, envir = e)
  assign(".sf_red", reduce, envir = e)
  assign(".sf_narm", isTRUE(na.rm), envir = e)
  assign(".sf_block", block, envir = e)
  assign(".sf_threads", threads, envir = e)
  environment(shell) <- e
  shell
}
