/* sf_eval.c -- elementwise evaluation and fused reductions.
 *
 * A reduction never allocates an output vector: the expression is evaluated
 * one cache-resident block at a time and collapsed in place, so
 * sf_sum(~(x - mu)^2) runs in O(1) auxiliary memory regardless of n.
 *
 * Work is cut into tasks of SF_TASK_ELEMS elements. The boundaries depend only
 * on n, never on the thread count, and task partials are combined serially in
 * index order in long double -- so a reduction returns the same bits on one
 * thread and on sixty. `omp parallel for reduction(+:)` cannot offer that.
 */

#include "sf_omp.h"   /* first: must precede any R header */
#include "sf_prog.h"
#include <math.h>
#include <string.h>


SEXP C_sf_eval(SEXP s_code, SEXP s_consts, SEXP s_vars, SEXP s_n, SEXP s_nreg,
               SEXP s_nprologue, SEXP s_block, SEXP s_fast, SEXP s_red,
               SEXP s_narm, SEXP s_threads)
{
  sf_prog P;
  sf_prog_init(&P, s_code, s_consts, s_vars, s_n, s_nreg, s_nprologue,
               s_block, s_fast);

  const int red = INTEGER(s_red)[0];
  const int na_rm = (Rf_asLogical(s_narm) == TRUE);
  if (red < 0 || red >= SF_RED_COUNT) Rf_error("simdfuse: bad reduction code");

  const R_xlen_t n = P.n;
  const R_xlen_t ntask = (n + SF_TASK_ELEMS - 1) / SF_TASK_ELEMS;
  const int nthreads = sf_prog_threads(n, ntask, INTEGER(s_threads)[0]);
  double *rf = sf_prog_regfile(&P, nthreads);
  const size_t per = (size_t)P.nreg * (size_t)P.B;
  (void)per;   /* only read in the OpenMP branch */

  SEXP out = R_NilValue;
  double *outp = NULL;
  double *part = NULL;
  int nprot = 0;

  if (red == SF_RED_NONE) {
    out = PROTECT(Rf_allocVector(REALSXP, n)); nprot++;
    outp = REAL(out);
  } else {
    part = (double *)R_alloc((size_t)(ntask > 0 ? ntask : 1) * 3, sizeof(double));
  }

  if (n > 0) {
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
#ifdef _OPENMP
      double *regs = rf + (size_t)omp_get_thread_num() * per;
#else
      double *regs = rf;
#endif
      double *fin = sf_prog_final(&P, regs);
      double blk[2];
      sf_prog_prologue(&P, regs);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (R_xlen_t t = 0; t < ntask; t++) {
        const R_xlen_t lo = t * SF_TASK_ELEMS;
        R_xlen_t hi = lo + SF_TASK_ELEMS;
        if (hi > n) hi = n;

        sf_acc A; sf_acc_init(&A, red);

        for (R_xlen_t s = lo; s < hi; s += P.B) {
          const int blen = (int)((hi - s < P.B) ? (hi - s) : P.B);
          sf_prog_block(&P, regs, s, blen);
          if (red == SF_RED_NONE) {
            memcpy(outp + s, fin, (size_t)blen * sizeof(double));
          } else {
            sf_prog_reduce(fin, blen, red, na_rm, blk);
            sf_acc_push(&A, red, blk);
          }
        }
        if (red != SF_RED_NONE) {
          part[t * 3]     = (double)A.v;
          part[t * 3 + 1] = A.aux;
          part[t * 3 + 2] = A.aux2;
        }
      }
    }
  }

  if (red == SF_RED_NONE) {
    UNPROTECT(nprot);
    return out;
  }

  /* Serial combine in task order: the step that buys reproducibility. */
  sf_acc T; sf_acc_init(&T, red);
  for (R_xlen_t t = 0; t < ntask; t++) {
    double blk[2];
    blk[0] = part[t * 3];
    blk[1] = (red == SF_RED_ANY || red == SF_RED_ALL) ? part[t * 3 + 2]
                                                      : part[t * 3 + 1];
    if (red == SF_RED_ANY || red == SF_RED_ALL) {
      T.aux  += part[t * 3 + 1];
      T.aux2 += part[t * 3 + 2];
    } else {
      sf_acc_push(&T, red, blk);
    }
  }

  switch (red) {
  case SF_RED_SUM:  return Rf_ScalarReal((double)T.v);
  case SF_RED_PROD: return Rf_ScalarReal((double)T.v);
  case SF_RED_NOBS: return Rf_ScalarReal((double)T.v);
  case SF_RED_MEAN: {
    double cnt = 0.0;
    for (R_xlen_t t = 0; t < ntask; t++) cnt += part[t * 3 + 1];
    if (cnt <= 0.0) return Rf_ScalarReal(R_NaN);
    return Rf_ScalarReal((double)(T.v / (long double)cnt));
  }
  case SF_RED_MAX:
  case SF_RED_MIN: {
    if (T.aux != 0.0) return Rf_ScalarReal(NA_REAL);
    double v = (double)T.v;
    if (n == 0 || (v == R_NegInf && red == SF_RED_MAX) ||
        (v == R_PosInf && red == SF_RED_MIN))
      Rf_warning("no non-missing arguments to %s; returning %s",
                 red == SF_RED_MAX ? "max" : "min",
                 red == SF_RED_MAX ? "-Inf" : "Inf");
    return Rf_ScalarReal(v);
  }
  case SF_RED_ANY:
    if (T.aux  > 0.0) return Rf_ScalarLogical(TRUE);
    if (T.aux2 > 0.0 && !na_rm) return Rf_ScalarLogical(NA_LOGICAL);
    return Rf_ScalarLogical(FALSE);
  case SF_RED_ALL:
    if (T.aux  > 0.0) return Rf_ScalarLogical(FALSE);
    if (T.aux2 > 0.0 && !na_rm) return Rf_ScalarLogical(NA_LOGICAL);
    return Rf_ScalarLogical(TRUE);
  default:
    Rf_error("simdfuse: unreachable reduction %d", red);
  }
  return R_NilValue;
}
