#' simdfuse: fused SIMD evaluation of array expressions
#'
#' @description
#' Evaluates elementwise R expressions in a single pass over the data, keeping
#' every intermediate value inside a cache-resident block instead of allocating
#' a temporary vector per operation. See [sf()] to get started, [sf_sum()] for
#' fused reductions, and [sf_scan()], [sf_stencil()] and [sf_segreduce()] for
#' the parallel skeletons.
#'
#' @section Entry points:
#' `vignette("simdfuse")` for an introduction, `vignette("performance")` for
#' measured timings, `vignette("accuracy")` for exactly what is and is not
#' guaranteed about the numbers.
#'
#' @keywords internal
#' @useDynLib simdfuse, .registration = TRUE, .fixes = ""
"_PACKAGE"
