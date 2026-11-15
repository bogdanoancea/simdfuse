# The expression compiler: an R call tree becomes a straight-line program for
# the block VM.
#
# Three things happen here that matter for performance and for correctness:
#
#  * Common subexpressions are shared. The call tree is interned into a DAG
#    keyed on canonical form, so (x-m)/s appearing twice is computed once.
#  * Registers are allocated by use counting. A node's register is returned to
#    the free list only when its last consumer has been emitted, which keeps
#    the register file -- and therefore the cache footprint -- minimal.
#  * Fused multiply-add is emitted only in "fast" math mode. Contracting
#    a*b+c into one instruction removes an intermediate rounding, so it would
#    silently break bit-identity with base R if it were always on.
#  * Block-invariant subexpressions are hoisted into a prologue that runs once
#    per call instead of once per block. This matters more than it sounds:
#    `pi` is an ordinary length-one variable in R, not a literal, so without
#    hoisting an expression like -0.5*log(2*pi) - log(s) would recompute two
#    logarithms for every element of the data.

.sf_fun_table <- function() {
  if (!is.null(.sf$funtab)) return(.sf$funtab)
  t <- list()
  bin <- function(nm, op) t[[nm]] <<- list(op = op, arity = 2L, kind = "bin")
  un  <- function(nm, op, inexact = FALSE)
    t[[nm]] <<- list(op = op, arity = 1L, kind = "un", inexact = inexact)

  bin("+", "ADD"); bin("-", "SUB"); bin("*", "MUL"); bin("/", "DIV")
  bin("%%", "MOD"); bin("%/%", "IDIV")
  bin("<", "LT"); bin("<=", "LE"); bin(">", "GT"); bin(">=", "GE")
  bin("==", "EQ"); bin("!=", "NE")
  bin("&", "AND"); bin("|", "OR")
  bin("atan2", "ATAN2")

  un("sqrt", "SQRT"); un("exp", "EXP"); un("abs", "ABS"); un("sign", "SIGN")
  un("floor", "FLOOR"); un("ceiling", "CEILING"); un("trunc", "TRUNC")
  un("log1p", "LOG1P"); un("expm1", "EXPM1")
  un("log2", "LOG2"); un("log10", "LOG10")
  un("sin", "SIN"); un("cos", "COS"); un("tan", "TAN")
  un("asin", "ASIN"); un("acos", "ACOS"); un("atan", "ATAN")
  un("sinh", "SINH"); un("cosh", "COSH"); un("tanh", "TANH")
  un("lgamma", "LGAMMA"); un("gamma", "GAMMA"); un("digamma", "DIGAMMA")
  un("!", "NOT")
  un("is.na", "ISNA"); un("is.finite", "ISFINITE"); un("is.infinite", "ISINF")

  # special forms, handled before the table lookup
  t[["^"]]      <- list(arity = 2L, kind = "special")
  t[["log"]]    <- list(arity = c(1L, 2L), kind = "special")
  t[["round"]]  <- list(arity = c(1L, 2L), kind = "special")
  t[["ifelse"]] <- list(arity = 3L, kind = "special")
  t[["pmin"]]   <- list(arity = "2+", kind = "special")
  t[["pmax"]]   <- list(arity = "2+", kind = "special")
  t[["pnorm"]]  <- list(arity = c(1L, 2L, 3L), kind = "special")
  t[["dnorm"]]  <- list(arity = c(1L, 2L, 3L), kind = "special")
  t[["plogis"]] <- list(arity = c(1L, 2L, 3L), kind = "special")
  t[["qlogis"]] <- list(arity = c(1L, 2L, 3L), kind = "special")
  t[["("]]          <- list(arity = 1L, kind = "special")
  t[["as.numeric"]] <- list(arity = 1L, kind = "special")
  t[["as.double"]]  <- list(arity = 1L, kind = "special")

  .sf$funtab <- t
  t
}

# Reductions are not elementwise; naming one inside an expression is a mistake
# worth catching with a pointed message rather than a generic one.
.sf_reduction_names <- c("sum", "mean", "max", "min", "prod", "any", "all",
                         "var", "sd", "median", "cumsum", "cumprod",
                         "cummax", "cummin", "length", "which")

.sf_is_literal <- function(e) {
  if (is.numeric(e) || is.logical(e)) return(length(e) == 1L && !is.na(e))
  if (is.call(e)) {
    fn <- e[[1L]]
    if (!is.name(fn)) return(FALSE)
    if (is.null(.sf_fun_table()[[as.character(fn)]])) return(FALSE)
    if (length(e) == 1L) return(FALSE)
    return(all(vapply(as.list(e)[-1L], .sf_is_literal, TRUE)))
  }
  FALSE
}

#' Compile an expression without evaluating it
#'
#' @description
#' Exposed mainly for inspection and teaching: it shows how many instructions
#' and registers an expression costs, which is the quantity that decides
#' whether the working set stays in L1 cache.
#'
#' @param expr An expression, or a one-sided formula.
#' @param math `"accurate"` or `"fast"`; see [sf_options()]. Affects whether
#'   fused multiply-add is emitted.
#' @param scalars Names of inputs that will be length one. Supplying them lets
#'   the compiler hoist the subexpressions that depend only on those inputs out
#'   of the per-block loop. `sf()` and friends work this out for themselves at
#'   call time; the argument exists so that `sf_compile()` can show the same
#'   program they will actually run.
#' @return An object of class `"sf_program"`.
#' @export
#' @examples
#' sf_compile(~ a * b + sqrt(c) - 2.5)
sf_compile <- function(expr, math = getOption("simdfuse.math"),
                       scalars = character(0)) {
  ex <- .sf_rhs(substitute(expr), parent.frame())$ex
  math <- match.arg(math, c("accurate", "fast"))
  .sf_compile_ast(ex, fast = identical(math, "fast"), scalars = scalars)
}

.sf_compile_ast <- function(ex, fast = FALSE, scalars = character(0)) {
  tab <- .sf_fun_table()

  nodes <- list()
  keymap <- new.env(parent = emptyenv(), hash = TRUE)
  consts <- numeric(0)
  constmap <- new.env(parent = emptyenv(), hash = TRUE)
  varnames <- character(0)

  const_idx <- function(v) {
    v <- as.double(v)
    key <- sprintf("%a", v)
    hit <- constmap[[key]]
    if (!is.null(hit)) return(hit)
    consts[[length(consts) + 1L]] <<- v
    idx <- length(consts) - 1L
    assign(key, idx, envir = constmap)
    idx
  }

  var_idx <- function(nm) {
    hit <- match(nm, varnames)
    if (!is.na(hit)) return(hit - 1L)
    varnames[[length(varnames) + 1L]] <<- nm
    length(varnames) - 1L
  }

  mknode <- function(key, op, kind, args = integer(0), imm = 0L) {
    hit <- keymap[[key]]
    if (!is.null(hit)) return(hit)
    nodes[[length(nodes) + 1L]] <<- list(op = .sf_opcode(op), kind = kind,
                                         args = args, imm = imm)
    id <- length(nodes)
    assign(key, id, envir = keymap)
    id
  }

  n_loadc <- function(v) {
    v <- as.double(v)
    mknode(paste0("c:", sprintf("%a", v)), "LOADC", "loadc", imm = const_idx(v))
  }
  n_loadv <- function(nm)
    mknode(paste0("v:", nm), "LOADV", "loadv", imm = var_idx(nm))

  n_un <- function(op, a)
    mknode(paste0(op, "(", a, ")"), op, "un", args = c(a))
  n_bin <- function(op, a, b)
    mknode(paste0(op, "(", a, ",", b, ")"), op, "bin", args = c(a, b))
  n_ter <- function(op, a, b, c)
    mknode(paste0(op, "(", a, ",", b, ",", c, ")"), op, "ter", args = c(a, b, c))
  n_powi <- function(a, k)
    mknode(paste0("POWI(", a, ",", k, ")"), "POWI", "powi", args = c(a),
           imm = const_idx(k))

  bad <- function(e, why)
    stop(sprintf("simdfuse: cannot compile `%s`: %s",
                 paste(deparse(e), collapse = " "), why), call. = FALSE)

  walk <- function(e) {
    if (is.numeric(e) || is.logical(e)) {
      if (length(e) != 1L) bad(e, "literal vectors must be length 1; pass them as variables")
      if (is.na(e)) return(n_loadc(NA_real_))
      return(n_loadc(e))
    }
    if (is.name(e)) {
      nm <- as.character(e)
      if (!nzchar(nm)) bad(e, "empty symbol")
      return(n_loadv(nm))
    }
    if (!is.call(e)) bad(e, "unsupported expression type")

    fn <- e[[1L]]
    if (!is.name(fn))
      bad(e, "only calls to named functions can be compiled")
    nm <- as.character(fn)
    args <- as.list(e)[-1L]

    if (nm %in% .sf_reduction_names)
      bad(e, sprintf("`%s()` is a reduction, not an elementwise operation; use sf_%s() or sf_scan()",
                     nm, nm))

    spec <- tab[[nm]]
    if (is.null(spec))
      bad(e, sprintf("`%s()` is not supported; see sf_supported()", nm))

    if (any(vapply(args, function(z) is.name(z) && !nzchar(as.character(z)), TRUE)))
      bad(e, "missing arguments are not allowed")
    if (!is.null(names(args)) && any(nzchar(names(args))))
      bad(e, "named arguments are not supported; pass them positionally")

    if (.sf_is_literal(e)) return(n_loadc(eval(e, baseenv())))

    # ---- special forms
    if (nm %in% c("(", "as.numeric", "as.double")) return(walk(args[[1L]]))

    if (nm %in% c("+", "-") && length(args) == 1L) {
      a <- walk(args[[1L]])
      if (nm == "+") return(a)
      return(n_un("NEG", a))
    }

    if (nm == "^") {
      if (length(args) != 2L) bad(e, "`^` takes two arguments")
      k <- args[[2L]]
      lit <- is.numeric(k) && length(k) == 1L && !is.na(k) && k == trunc(k)
      # R's arithmetic special-cases only y == 2 (as x * x) and otherwise
      # calls the platform pow(), so repeated squaring would silently diverge
      # in the last bits. In "accurate" mode POWI is therefore used only for
      # the exponents where it provably agrees; "fast" mode takes the win.
      if (lit && (abs(k) <= 1024) && (fast || k %in% c(0, 1, 2)))
        return(n_powi(walk(args[[1L]]), as.double(k)))
      return(n_bin("POW", walk(args[[1L]]), walk(k)))
    }

    if (nm == "log") {
      a <- walk(args[[1L]])
      if (length(args) == 1L) return(n_un("LOG", a))
      bs <- args[[2L]]
      # R's logbase() routes base 2 and base 10 to log2() and log10() rather
      # than to a ratio of logarithms, and the two differ in the last bits.
      if (is.numeric(bs) && length(bs) == 1L && !is.na(bs)) {
        if (bs == 2)  return(n_un("LOG2", a))
        if (bs == 10) return(n_un("LOG10", a))
      }
      return(n_bin("DIV", n_un("LOG", a), n_un("LOG", walk(bs))))
    }

    if (nm == "round") {
      a <- walk(args[[1L]])
      d <- if (length(args) == 2L) walk(args[[2L]]) else n_loadc(0)
      return(n_bin("ROUND2", a, d))
    }

    if (nm == "ifelse") {
      if (length(args) != 3L) bad(e, "`ifelse()` takes three arguments")
      return(n_ter("SELECT", walk(args[[1L]]), walk(args[[2L]]), walk(args[[3L]])))
    }

    if (nm %in% c("pmin", "pmax")) {
      if (length(args) < 2L) bad(e, sprintf("`%s()` needs at least two arguments", nm))
      op <- if (nm == "pmin") "MIN2" else "MAX2"
      acc <- walk(args[[1L]])
      for (i in 2:length(args)) acc <- n_bin(op, acc, walk(args[[i]]))
      return(acc)
    }

    if (nm %in% c("pnorm", "dnorm", "plogis", "qlogis")) {
      # Location and scale are passed through to R's own routines rather than
      # folded into a standardised argument: R's dnorm() reorders the division
      # by sigma for large |z|, so standardising first would not be exact.
      q <- walk(args[[1L]])
      loc <- if (length(args) >= 2L) walk(args[[2L]]) else n_loadc(0)
      scl <- if (length(args) >= 3L) walk(args[[3L]]) else n_loadc(1)
      return(n_ter(toupper(nm), q, loc, scl))
    }

    # ---- table-driven forms
    if (spec$kind == "un") {
      if (length(args) != 1L) bad(e, sprintf("`%s()` takes one argument", nm))
      return(n_un(spec$op, walk(args[[1L]])))
    }
    if (spec$kind == "bin") {
      if (length(args) != 2L) bad(e, sprintf("`%s` takes two arguments", nm))
      a <- walk(args[[1L]]); b <- walk(args[[2L]])
      if (fast && spec$op == "ADD") {
        # contract x*y + z into a single fused instruction
        na <- nodes[[a]]; nb <- nodes[[b]]
        if (identical(na$op, .sf_opcode("MUL")))
          return(n_ter("FMA3", na$args[1L], na$args[2L], b))
        if (identical(nb$op, .sf_opcode("MUL")))
          return(n_ter("FMA3", nb$args[1L], nb$args[2L], a))
      }
      return(n_bin(spec$op, a, b))
    }
    bad(e, "unsupported form")
  }

  root <- walk(ex)

  if (length(varnames) == 0L)
    stop("simdfuse: the expression has no variables, so its length is undefined",
         call. = FALSE)

  # ---- mark block-invariant nodes
  # A node is invariant if every leaf below it is a literal or a length-one
  # input, so its value is the same in every block of every call.
  nn <- length(nodes)
  inv <- logical(nn)
  for (id in seq_len(nn)) {
    nd <- nodes[[id]]
    inv[id] <- switch(nd$kind,
      loadc = TRUE,
      loadv = varnames[nd$imm + 1L] %in% scalars,
      all(inv[nd$args]))
  }
  # Invariance is closed downwards, so emitting all invariant nodes first
  # preserves every dependency.
  order <- c(which(inv), which(!inv))
  nprologue <- sum(inv)

  # ---- linearise with use-counted register allocation
  # An invariant register must stay live for the whole call if any per-block
  # instruction reads it, because the prologue runs only once. If its consumers
  # are all inside the prologue, the register can be recycled like any other --
  # which keeps the working set small, and the working set is the whole point.
  inv_consumers <- integer(nn)
  pinned <- logical(nn)
  for (id in seq_len(nn)) {
    for (a in nodes[[id]]$args) {
      if (inv[id]) inv_consumers[a] <- inv_consumers[a] + 1L
      else pinned[a] <- TRUE
    }
  }
  rem <- integer(nn)
  for (nd in nodes) for (a in nd$args) rem[a] <- rem[a] + 1L

  reg <- integer(nn)
  free <- integer(0)
  nreg <- 0L
  code <- integer(5L * nn)
  pos <- 0L

  for (id in order) {
    nd <- nodes[[id]]
    if (length(free)) {
      d <- free[1L]; free <- free[-1L]
    } else {
      d <- nreg; nreg <- nreg + 1L
    }
    reg[id] <- d
    a <- b <- cc <- 0L
    switch(nd$kind,
      loadv = { a <- nd$imm },
      loadc = { a <- nd$imm },
      un    = { a <- reg[nd$args[1L]]; b <- a; cc <- a },
      powi  = { a <- reg[nd$args[1L]]; b <- nd$imm; cc <- a },
      bin   = { a <- reg[nd$args[1L]]; b <- reg[nd$args[2L]]; cc <- a },
      ter   = { a <- reg[nd$args[1L]]; b <- reg[nd$args[2L]]
                cc <- reg[nd$args[3L]] },
      stop("simdfuse: internal error, unknown node kind")
    )
    o <- 5L * pos
    pos <- pos + 1L
    code[o + 1L] <- nd$op
    code[o + 2L] <- as.integer(a)
    code[o + 3L] <- as.integer(b)
    code[o + 4L] <- as.integer(cc)
    code[o + 5L] <- as.integer(d)

    for (ch in nd$args) {
      rem[ch] <- rem[ch] - 1L
      if (inv[ch]) {
        inv_consumers[ch] <- inv_consumers[ch] - 1L
        if (inv_consumers[ch] == 0L && !pinned[ch]) free <- c(reg[ch], free)
      } else if (rem[ch] == 0L) {
        free <- c(reg[ch], free)
      }
    }
  }

  logical_ops <- c("LT", "LE", "GT", "GE", "EQ", "NE", "AND", "OR", "NOT",
                   "ISNA", "ISFINITE", "ISINF")
  root_op <- names(.sf$ops)[nodes[[root]]$op + 1L]

  structure(list(
    code      = code,
    logical   = root_op %in% logical_ops,
    consts    = as.double(consts),
    varnames  = varnames,
    scalars   = scalars,
    nreg      = as.integer(max(nreg, 1L)),
    ninstr    = nn,
    nprologue = as.integer(nprologue),
    fast      = fast,
    expr      = ex
  ), class = "sf_program")
}

#' @export
print.sf_program <- function(x, ...) {
  cat("<sf_program>\n")
  cat("  expression : ", paste(deparse(x$expr), collapse = " "), "\n", sep = "")
  cat("  inputs     : ", paste(x$varnames, collapse = ", "), "\n", sep = "")
  cat("  instructions: ", x$ninstr,
      if (x$nprologue > 0L) sprintf(" (%d hoisted out of the block loop)",
                                    x$nprologue) else "", "\n", sep = "")
  cat("  registers  : ", x$nreg, "\n", sep = "")
  cat("  constants  : ", length(x$consts), "\n", sep = "")
  cat("  math       : ", if (x$fast) "fast" else "accurate", "\n", sep = "")
  blk <- getOption("simdfuse.block")
  cat("  working set: ", format(x$nreg * blk * 8 / 1024, digits = 3),
      " KiB at block = ", blk, "\n", sep = "")
  invisible(x)
}

#' Disassemble a compiled program
#'
#' @param x An `"sf_program"` from [sf_compile()].
#' @return A data frame with one row per instruction, invisibly printed.
#' @export
#' @examples
#' sf_disassemble(sf_compile(~ (x - mu)^2 / s))
sf_disassemble <- function(x) {
  stopifnot(inherits(x, "sf_program"))
  nms <- names(.sf$ops)
  m <- matrix(x$code, ncol = 5L, byrow = TRUE)
  op <- nms[m[, 1L] + 1L]
  arg <- character(nrow(m))
  for (i in seq_len(nrow(m))) {
    arg[i] <- switch(op[i],
      LOADV = x$varnames[m[i, 2L] + 1L],
      LOADC = format(x$consts[m[i, 2L] + 1L]),
      POWI  = sprintf("r%d ^ %s", m[i, 2L], format(x$consts[m[i, 3L] + 1L])),
      SELECT = sprintf("r%d ? r%d : r%d", m[i, 2L], m[i, 3L], m[i, 4L]),
      FMA3   = sprintf("r%d * r%d + r%d", m[i, 2L], m[i, 3L], m[i, 4L]),
      if (m[i, 2L] == m[i, 3L]) sprintf("r%d", m[i, 2L])
      else sprintf("r%d, r%d", m[i, 2L], m[i, 3L])
    )
  }
  data.frame(i = seq_len(nrow(m)), op = op, args = arg,
             dst = paste0("r", m[, 5L]),
             once = seq_len(nrow(m)) <= x$nprologue,
             row.names = NULL, stringsAsFactors = FALSE)
}
