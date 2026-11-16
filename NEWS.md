# simdfuse 0.1.0

First release.

* `sf()` and `sf_eval()` evaluate an elementwise R expression in a single fused
  pass, keeping every intermediate inside a cache-resident block instead of
  allocating a temporary vector per operation. The expression is compiled to a
  compact instruction stream and executed by a kernel built once per
  instruction set -- AVX-512, AVX2 with FMA, SSE2, NEON and a scalar
  fallback -- and selected at run time from the host's reported features.

* Fused reductions (`sf_sum()`, `sf_mean()`, `sf_max()`, `sf_min()`,
  `sf_prod()`, `sf_any()`, `sf_all()`, `sf_nobs()`) consume an expression
  rather than a materialised vector, so they run in constant auxiliary memory.
  Their parallel form is deterministic in the thread count: task boundaries
  depend only on `n`, and task partials combine serially in index order.

* Parallel skeletons for which base R has no construct: `sf_scan()` for an
  associative prefix scan, including a numerically stable log-domain scan;
  `sf_stencil()` for a weighted stencil; and `sf_segreduce()` for reduction
  within contiguous segments.

* `sf_kernel()` and `sf_compile()` expose a compiled expression as a reusable
  function; `sf_disassemble()` prints the instruction stream.

* `sf_tune()` measures the host and selects the block size and the per-function
  fast-math settings for it.

* `math = "accurate"` reproduces base R bit for bit. The one exception is
  `sqrt()`, which every backend takes from the hardware: IEEE 754 requires that
  to be correctly rounded, whereas base R calls the C library, and not every C
  library is. Where the two disagree simdfuse returns the correctly rounded
  value. See `vignette("accuracy")`.
