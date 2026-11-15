/* sf_skeletons.c -- parallel algorithmic skeletons that base R has no
 * construct for: an associative scan, a weighted stencil, and a reduction
 * over contiguous segments.
 *
 * What makes these more than reimplementations of cumsum() and filter() is
 * that each consumes a *program*, not a vector. sf_scan(~ log(x / lag)) or
 * sf_segreduce(~ (y - yhat)^2, firm) never materialises the transformed
 * vector: the expression is evaluated one cache-resident block at a time and
 * consumed by the skeleton in the same pass.
 *
 * The scan is a textbook reduce-then-scan -- one sweep for per-task totals, a
 * serial exclusive prefix over those totals, a second sweep to write output.
 * Across tasks that costs a second evaluation of the expression, which is the
 * price of parallelism here and why parallel = FALSE stays available (and is
 * also the only mode that reproduces cumsum() bit for bit).
 */

#include "sf_omp.h"   /* first: must precede any R header */
#include "sf_prog.h"
#include <math.h>
#include <string.h>


enum { SF_SC_SUM = 0, SF_SC_PROD, SF_SC_MAX, SF_SC_MIN, SF_SC_LSE, SF_SC_COUNT };

static double sf_identity(int op) {
  switch (op) {
  case SF_SC_SUM:  return 0.0;
  case SF_SC_PROD: return 1.0;
  case SF_SC_MAX:  return R_NegInf;
  case SF_SC_MIN:  return R_PosInf;
  case SF_SC_LSE:  return R_NegInf;
  default:         return R_NaN;
  }
}

/* Combine propagates NA explicitly: the hardware min/max instructions do not,
 * and R's cummax/cummin do. */
static double sf_combine(int op, double a, double b) {
  if (ISNAN(a) || ISNAN(b)) return a + b;      /* carries the NA payload */
  switch (op) {
  case SF_SC_SUM:  return a + b;
  case SF_SC_PROD: return a * b;
  case SF_SC_MAX:  return a > b ? a : b;
  case SF_SC_MIN:  return a < b ? a : b;
  case SF_SC_LSE: {
    if (a == R_NegInf) return b;
    if (b == R_NegInf) return a;
    const double m = a > b ? a : b;
    const double d = a > b ? b - a : a - b;
    return m + log1p(exp(d));
  }
  default: return R_NaN;
  }
}

/* ------------------------------------------------------------------- scan */

SEXP C_sf_scan(SEXP s_code, SEXP s_consts, SEXP s_vars, SEXP s_n, SEXP s_nreg,
               SEXP s_nprologue, SEXP s_block, SEXP s_fast, SEXP s_op,
               SEXP s_threads)
{
  sf_prog P;
  sf_prog_init(&P, s_code, s_consts, s_vars, s_n, s_nreg, s_nprologue,
               s_block, s_fast);
  const int op = INTEGER(s_op)[0];
  if (op < 0 || op >= SF_SC_COUNT) Rf_error("simdfuse: unknown scan operator");

  const R_xlen_t n = P.n;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *o = REAL(out);
  if (n == 0) { UNPROTECT(1); return out; }

  const R_xlen_t ntask = (n + SF_TASK_ELEMS - 1) / SF_TASK_ELEMS;
  const int nthreads = sf_prog_threads(n, ntask, INTEGER(s_threads)[0]);
  const double id = sf_identity(op);
  const size_t per = (size_t)P.nreg * (size_t)P.B;
  (void)per;   /* only read in the OpenMP branch */

  if (nthreads == 1) {
    double *rf = sf_prog_regfile(&P, 1);
    double *fin = sf_prog_final(&P, rf);
    sf_prog_prologue(&P, rf);
    if (op == SF_SC_SUM) {
      /* long double carry, matching base R's cumsum() */
      long double acc = 0.0L;
      for (R_xlen_t s = 0; s < n; s += P.B) {
        const int blen = (int)((n - s < P.B) ? (n - s) : P.B);
        sf_prog_block(&P, rf, s, blen);
        for (int j = 0; j < blen; j++) {
          acc += (long double)fin[j];
          o[s + j] = (double)acc;
        }
      }
    } else {
      double acc = id;
      for (R_xlen_t s = 0; s < n; s += P.B) {
        const int blen = (int)((n - s < P.B) ? (n - s) : P.B);
        sf_prog_block(&P, rf, s, blen);
        for (int j = 0; j < blen; j++) {
          acc = sf_combine(op, acc, fin[j]);
          o[s + j] = acc;
        }
      }
    }
    UNPROTECT(1);
    return out;
  }

  double *tot = (double *)R_alloc((size_t)ntask, sizeof(double));
  double *off = (double *)R_alloc((size_t)ntask, sizeof(double));
  double *rf = sf_prog_regfile(&P, nthreads);

  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1) {
      double run = id;
      for (R_xlen_t t = 0; t < ntask; t++) {
        off[t] = run;
        run = sf_combine(op, run, tot[t]);
      }
    }
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
      sf_prog_prologue(&P, regs);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (R_xlen_t t = 0; t < ntask; t++) {
        const R_xlen_t lo = t * SF_TASK_ELEMS;
        R_xlen_t hi = lo + SF_TASK_ELEMS; if (hi > n) hi = n;
        double acc = (pass == 0) ? id : off[t];
        for (R_xlen_t s = lo; s < hi; s += P.B) {
          const int blen = (int)((hi - s < P.B) ? (hi - s) : P.B);
          sf_prog_block(&P, regs, s, blen);
          if (pass == 0) {
            for (int j = 0; j < blen; j++) acc = sf_combine(op, acc, fin[j]);
          } else {
            for (int j = 0; j < blen; j++) {
              acc = sf_combine(op, acc, fin[j]);
              o[s + j] = acc;
            }
          }
        }
        if (pass == 0) tot[t] = acc;
      }
    }
  }

  UNPROTECT(1);
  return out;
}

/* ---------------------------------------------------------------- stencil */

enum { SF_BND_NA = 0, SF_BND_ZERO, SF_BND_REPLICATE, SF_BND_REFLECT,
       SF_BND_CIRCULAR, SF_BND_COUNT };

/* Map an out-of-range index onto a real one per the boundary rule.
 * Returns -1 when the rule supplies a constant instead, in *cval. */
static R_xlen_t sf_wrap(R_xlen_t k, R_xlen_t n, int bnd, double *cval) {
  if (k >= 0 && k < n) return k;
  switch (bnd) {
  case SF_BND_NA:        *cval = NA_REAL; return -1;
  case SF_BND_ZERO:      *cval = 0.0;     return -1;
  case SF_BND_REPLICATE: return (k < 0) ? 0 : n - 1;
  case SF_BND_REFLECT: {
    R_xlen_t r = (k < 0) ? (-k - 1) : (2 * n - k - 1);
    if (r < 0) r = 0;
    if (r >= n) r = n - 1;
    return r;
  }
  case SF_BND_CIRCULAR: {
    R_xlen_t r = k % n;
    if (r < 0) r += n;
    return r;
  }
  default: *cval = NA_REAL; return -1;
  }
}

SEXP C_sf_stencil(SEXP s_code, SEXP s_consts, SEXP s_vars, SEXP s_n, SEXP s_nreg,
                  SEXP s_nprologue, SEXP s_block, SEXP s_fast, SEXP s_w,
                  SEXP s_bnd, SEXP s_threads)
{
  sf_prog P;
  sf_prog_init(&P, s_code, s_consts, s_vars, s_n, s_nreg, s_nprologue,
               s_block, s_fast);
  const double *w = REAL(s_w);
  const int nw = (int)XLENGTH(s_w);
  const int bnd = INTEGER(s_bnd)[0];
  if (nw < 1) Rf_error("simdfuse: stencil needs at least one weight");
  if (nw % 2 == 0) Rf_error("simdfuse: stencil weights must have odd length");
  if (bnd < 0 || bnd >= SF_BND_COUNT) Rf_error("simdfuse: unknown boundary rule");

  const int half = nw / 2;
  const R_xlen_t n = P.n;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *o = REAL(out);
  if (n == 0) { UNPROTECT(1); return out; }

  /* Edge caches: the halo of the first and last task can reference indices
   * near either end of the vector (circular wraps to the far end), so those
   * few values are computed once, serially, up front. */
  const R_xlen_t edge = (half < n) ? half : n;
  double *head = (double *)R_alloc((size_t)(edge > 0 ? edge : 1), sizeof(double));
  double *tail = (double *)R_alloc((size_t)(edge > 0 ? edge : 1), sizeof(double));
  {
    double *rf1 = sf_prog_regfile(&P, 1);
    double *fin = sf_prog_final(&P, rf1);
    sf_prog_prologue(&P, rf1);
    for (R_xlen_t s = 0; s < edge; s += P.B) {
      const int blen = (int)((edge - s < P.B) ? (edge - s) : P.B);
      sf_prog_block(&P, rf1, s, blen);
      memcpy(head + s, fin, (size_t)blen * sizeof(double));
    }
    const R_xlen_t tstart = n - edge;
    for (R_xlen_t s = 0; s < edge; s += P.B) {
      const int blen = (int)((edge - s < P.B) ? (edge - s) : P.B);
      sf_prog_block(&P, rf1, tstart + s, blen);
      memcpy(tail + s, fin, (size_t)blen * sizeof(double));
    }
  }

  const R_xlen_t ntask = (n + SF_TASK_ELEMS - 1) / SF_TASK_ELEMS;
  const int nthr = sf_prog_threads(n, ntask, INTEGER(s_threads)[0]);
  const size_t per = (size_t)P.nreg * (size_t)P.B;
  (void)per;   /* only read in the OpenMP branch */
  const size_t scr = (size_t)SF_TASK_ELEMS + (size_t)nw - 1;
  double *rf = sf_prog_regfile(&P, nthr);
  double *scratch = (double *)R_alloc(scr * (size_t)nthr, sizeof(double));

#ifdef _OPENMP
#pragma omp parallel num_threads(nthr)
#endif
  {
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    double *regs = rf + (size_t)tid * per;
    double *buf = scratch + (size_t)tid * scr;
    double *fin = sf_prog_final(&P, regs);
    sf_prog_prologue(&P, regs);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (R_xlen_t t = 0; t < ntask; t++) {
      const R_xlen_t lo = t * SF_TASK_ELEMS;
      R_xlen_t hi = lo + SF_TASK_ELEMS; if (hi > n) hi = n;
      const R_xlen_t base = lo - half;          /* buf[0] is logical index base */
      const R_xlen_t want = (hi + half) - base; /* values needed in buf */

      const R_xlen_t cs = (base < 0) ? 0 : base;
      R_xlen_t ce = hi + half; if (ce > n) ce = n;

      for (R_xlen_t s = cs; s < ce; s += P.B) {
        const int blen = (int)((ce - s < P.B) ? (ce - s) : P.B);
        sf_prog_block(&P, regs, s, blen);
        memcpy(buf + (s - base), fin, (size_t)blen * sizeof(double));
      }
      /* halo positions, resolved through the boundary rule */
      for (R_xlen_t q = 0; q < want; q++) {
        const R_xlen_t idx = base + q;
        if (idx >= cs && idx < ce) continue;
        double cval = NA_REAL;
        const R_xlen_t k = sf_wrap(idx, n, bnd, &cval);
        if (k < 0) { buf[q] = cval; continue; }
        if (k >= cs && k < ce)      buf[q] = buf[k - base];
        else if (k < edge)          buf[q] = head[k];
        else if (k >= n - edge)     buf[q] = tail[k - (n - edge)];
        else                        buf[q] = NA_REAL;  /* unreachable */
      }
      sf_prog_stencil(buf, hi - lo, w, nw, o + lo);
    }
  }

  UNPROTECT(1);
  return out;
}

/* -------------------------------------------------------- segmented reduce */

SEXP C_sf_segreduce(SEXP s_code, SEXP s_consts, SEXP s_vars, SEXP s_n,
                    SEXP s_nreg, SEXP s_nprologue, SEXP s_block, SEXP s_fast,
                    SEXP s_seg, SEXP s_red, SEXP s_narm, SEXP s_threads)
{
  sf_prog P;
  sf_prog_init(&P, s_code, s_consts, s_vars, s_n, s_nreg, s_nprologue,
               s_block, s_fast);
  const int red = INTEGER(s_red)[0];
  const int na_rm = (Rf_asLogical(s_narm) == TRUE);
  if (red <= 0 || red >= SF_RED_COUNT) Rf_error("simdfuse: bad reduction code");

  const R_xlen_t n = P.n;
  if (XLENGTH(s_seg) != n)
    Rf_error("simdfuse: `by` has length %lld but the expression is length %lld",
             (long long)XLENGTH(s_seg), (long long)n);
  const int *seg = INTEGER(s_seg);

  /* Count contiguous runs. An NA always starts a new run. */
  R_xlen_t nrun = 0;
  for (R_xlen_t i = 0; i < n; i++) {
    if (i == 0 || seg[i] == NA_INTEGER || seg[i - 1] == NA_INTEGER ||
        seg[i] != seg[i - 1])
      nrun++;
  }

  SEXP ids  = PROTECT(Rf_allocVector(INTSXP, nrun));
  SEXP vals = PROTECT(Rf_allocVector(REALSXP, nrun));
  SEXP lens = PROTECT(Rf_allocVector(INTSXP, nrun));
  R_xlen_t *start = (R_xlen_t *)R_alloc((size_t)nrun + 1, sizeof(R_xlen_t));
  {
    R_xlen_t r = 0;
    for (R_xlen_t i = 0; i < n; i++) {
      if (i == 0 || seg[i] == NA_INTEGER || seg[i - 1] == NA_INTEGER ||
          seg[i] != seg[i - 1]) {
        start[r] = i;
        INTEGER(ids)[r] = seg[i];
        r++;
      }
    }
    start[nrun] = n;
    for (R_xlen_t k = 0; k < nrun; k++)
      INTEGER(lens)[k] = (int)(start[k + 1] - start[k]);
  }

  /* Runs are independent, so this parallelises without any merge step and
   * stays deterministic: within a run, blocks combine in index order. */
  const int nthr = sf_prog_threads(n, nrun, INTEGER(s_threads)[0]);
  const size_t per = (size_t)P.nreg * (size_t)P.B;
  (void)per;   /* only read in the OpenMP branch */
  double *rf = sf_prog_regfile(&P, nthr);
  double *vp = REAL(vals);

#ifdef _OPENMP
#pragma omp parallel num_threads(nthr)
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
#pragma omp for schedule(dynamic, 64)
#endif
    for (R_xlen_t r = 0; r < nrun; r++) {
      const R_xlen_t lo = start[r], hi = start[r + 1];
      sf_acc A; sf_acc_init(&A, red);
      for (R_xlen_t s = lo; s < hi; s += P.B) {
        const int blen = (int)((hi - s < P.B) ? (hi - s) : P.B);
        sf_prog_block(&P, regs, s, blen);
        sf_prog_reduce(fin, blen, red, na_rm, blk);
        sf_acc_push(&A, red, blk);
      }
      vp[r] = sf_acc_value(&A, red, na_rm);
    }
  }

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, ids);
  SET_VECTOR_ELT(res, 1, vals);
  SET_VECTOR_ELT(res, 2, lens);
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(nm, 0, Rf_mkChar("id"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("value"));
  SET_STRING_ELT(nm, 2, Rf_mkChar("length"));
  Rf_setAttrib(res, R_NamesSymbol, nm);
  UNPROTECT(5);
  return res;
}
