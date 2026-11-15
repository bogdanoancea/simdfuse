# Run a block of expectations once per instruction set the host supports, so a
# single CI machine exercises the scalar, SSE2, AVX2 and AVX-512 kernels.
for_each_isa <- function(code) {
  old <- sf_options()$isa
  on.exit(sf_options(isa = old), add = TRUE)
  for (isa in sf_info()$available) {
    sf_options(isa = isa)
    withCallingHandlers(code(isa), error = function(e) {
      stop(sprintf("failure under isa = '%s': %s", isa, conditionMessage(e)),
           call. = FALSE)
    })
  }
}

sf_testdata <- function(n = 5000L, seed = 20260101L) {
  set.seed(seed)
  list(
    n = n,
    a = runif(n, 0.1, 10),
    b = runif(n, 0.1, 10),
    c = runif(n, -5, 5),
    p = runif(n, 0.01, 0.99),
    g = sort(sample.int(max(1L, n %/% 50L), n, replace = TRUE))
  )
}

edge_values <- c(0, -0, 1, -1, 0.5, -0.5, NA_real_, NaN, Inf, -Inf,
                 .Machine$double.xmin, .Machine$double.xmax,
                 .Machine$double.eps, 1e-320, -1e-320, 2^53, -2^53)
