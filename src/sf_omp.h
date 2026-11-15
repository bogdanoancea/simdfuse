/* sf_omp.h -- OpenMP include, which must precede every R header.
 *
 * R's Rinternals.h defines an unprefixed macro `match` (to Rf_match), and
 * libomp's omp.h uses match(...) as an OpenMP clause keyword in its
 * `#pragma omp begin declare variant match(device={kind(host)})` lines.
 * Clang macro-expands inside OpenMP pragmas, so if omp.h is preprocessed
 * after Rinternals.h the clause becomes Rf_match(...) and compilation fails
 * with "expected 'match', 'adjust_args', or 'append_args' clause".
 *
 * Observed with Apple clang 21 and Homebrew libomp on R 4.6, aarch64.
 * GCC's own omp.h does not use that keyword, which is why the problem does
 * not show up on a Linux build.
 *
 * Include this file FIRST in any translation unit that uses OpenMP, ahead of
 * R.h, Rinternals.h and anything that might pull them in.
 */
#ifndef SF_OMP_H
#define SF_OMP_H

#ifdef _OPENMP
#include <omp.h>
#endif

#endif /* SF_OMP_H */
