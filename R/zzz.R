# Size of this CPU's L1 data cache, in bytes, best effort. Used to cap the
# block size by the number of registers a program needs: see .sf_block_cap().
# A wrong answer costs some performance and nothing else, so an undetectable
# machine gets the smallest common value rather than an error.
.sf_detect_l1 <- function() {
  v <- tryCatch({
    if (file.exists("/sys/devices/system/cpu/cpu0/cache")) {
      # pick the index whose level is 1 and type is Data (or Unified)
      idx <- list.files("/sys/devices/system/cpu/cpu0/cache",
                        pattern = "^index[0-9]+$", full.names = TRUE)
      sz <- NA_real_
      for (d in idx) {
        lv <- tryCatch(readLines(file.path(d, "level"), warn = FALSE)[1], error = function(e) NA)
        ty <- tryCatch(readLines(file.path(d, "type"),  warn = FALSE)[1], error = function(e) NA)
        if (identical(lv, "1") && ty %in% c("Data", "Unified")) {
          raw <- readLines(file.path(d, "size"), warn = FALSE)[1]
          num <- as.numeric(sub("[KkMm]?B?$", "", raw))
          mult <- if (grepl("[Mm]", raw)) 1024^2 else if (grepl("[Kk]", raw)) 1024 else 1
          sz <- num * mult
          break
        }
      }
      sz
    } else if (Sys.info()[["sysname"]] == "Darwin") {
      as.numeric(system2("sysctl", c("-n", "hw.l1dcachesize"), stdout = TRUE, stderr = FALSE))
    } else NA_real_
  }, error = function(e) NA_real_, warning = function(w) NA_real_)
  if (length(v) != 1L || is.na(v) || !is.finite(v) || v < 4096 || v > 4 * 1024^2)
    return(32768L)                 # smallest in common use; conservative
  as.integer(v)
}

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
    simdfuse.l1_bytes  = .sf_detect_l1(),
    simdfuse.l1_share  = 1/3,
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
