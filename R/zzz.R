# Package state. Opcode and reduction tables are read from the compiled engine
# at load time rather than duplicated in R, so the two can never drift apart.

.sf <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  .sf$ops  <- .Call(C_sf_op_table)
  .sf$reds <- .Call(C_sf_red_table)
  .sf$info <- .Call(C_sf_info)
  .sf$cache <- new.env(parent = emptyenv())
  .sf$cache_order <- character(0)

  op <- options()
  defaults <- list(
    simdfuse.block     = 128L,
    simdfuse.block_red = NA_integer_,
    simdfuse.threads   = 0L,       # 0 = every thread OpenMP offers
    simdfuse.isa       = "auto",
    simdfuse.math      = "accurate",
    simdfuse.fast_exp  = TRUE,
    simdfuse.fast_log  = TRUE,
    simdfuse.cache_max = 256L
  )
  toset <- !(names(defaults) %in% names(op))
  if (any(toset)) options(defaults[toset])
  invisible()
}

.sf_opcode <- function(name) {
  v <- .sf$ops[[name]]
  if (is.null(v)) stop(sprintf("simdfuse: internal error, unknown opcode '%s'", name))
  as.integer(v)
}

.sf_redcode <- function(name) {
  v <- .sf$reds[[name]]
  if (is.null(v)) stop(sprintf("simdfuse: internal error, unknown reduction '%s'", name))
  as.integer(v)
}
