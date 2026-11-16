## Test environments

* macOS Tahoe 26.6.2, aarch64-apple-darwin23, R 4.6.1, Apple clang 21
  -- NEON and scalar backends compiled and exercised
* Windows 11, x86_64-w64-mingw32, R 4.5.2 (Rtools45)
  -- scalar, SSE2 and AVX2 backends compiled and exercised
* Linux, x86_64 (AMD EPYC 7713), R 4.5.x
  -- scalar, SSE2 and AVX2 backends compiled and exercised
* win-builder, R-devel and R-release  <RUN BEFORE SUBMITTING>

## R CMD check results

macOS, R 4.6.1, `--as-cran` with `_R_CHECK_LIMIT_CORES_=TRUE`:

    0 errors | 0 warnings | 1 note

    * checking CRAN incoming feasibility ... NOTE
      New submission

Other platforms: <FILL IN -- this file should not be submitted with a placeholder>

## Notes for the reviewer

This is a new submission.

The package compiles one copy of its kernel per instruction set using
`__attribute__((target(...)))` and chooses between them at run time from
`__builtin_cpu_supports()`. Every x86 path is behind
`#if defined(__x86_64__) || defined(__i386__)`, NEON behind `__aarch64__`, and
a plain scalar kernel is always built, so a platform with neither set of
intrinsics still gets a working package.

Examples, tests and vignettes use at most two cores. The default thread count
is capped at two when `_R_CHECK_PACKAGE_NAME_` or `_R_CHECK_LIMIT_CORES_` is
set; the benchmark vignette's thread sweep stops at two under the same
condition; and the determinism tests sweep beyond two threads only when
`NOT_CRAN` is set.

The performance vignette measures on the machine that builds it rather than
reporting stored numbers. It uses its smallest size profile unless `NOT_CRAN`
is set, which holds it to a few seconds.
