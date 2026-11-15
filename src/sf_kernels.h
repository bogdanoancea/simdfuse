/* sf_kernels.h -- the block executor, compiled once per instruction set.
 *
 * Included by sf_kernels_<isa>.c after SF_ISA is defined. Everything here is
 * written against the macros in sf_simd.h, so a single source produces the
 * scalar, SSE2, AVX2, AVX-512 and NEON variants with no duplicated logic.
 *
 * Registers are plain double[B] arrays. That lets the rarely used elementary
 * functions fall back to a scalar libm loop over the same storage with no data
 * movement, which is what makes the bit-exact "accurate" math mode cheap to
 * support alongside the vectorized "fast" mode.
 *
 * Instruction layout: (op, a, b, c, dst), five ints.
 *
 * Editing this file requires `R CMD INSTALL --preclean`: make tracks only the
 * .c wrappers, not this body.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>
#include "sf_ops.h"
#include "sf_simd.h"
#include "sf_math.h"

#define SF_IARG(k) (code[ip*SF_INSTR_WIDTH + (k)])

/* Elements handled per iteration of the elementwise loops: one 64-byte cache
 * line. The inner loop has a compile-time trip count, so it unrolls fully --
 * eight lanes of work per iteration on every backend rather than one vector.
 * Without this the hand-written two-lane NEON kernels were *slower* than the
 * scalar kernel, because clang unrolls the plain C loop and we did not. */
#define SF_STEP 8


#define SF_BIN(VOP) do { \
    const double *a_ = ra, *b_ = rb; double *d_ = rd; \
    int j = 0; \
    for (; j + SF_STEP <= bpad; j += SF_STEP) \
      for (int u = 0; u < SF_STEP; u += SF_W) \
        V_STORE(d_ + j + u, VOP(V_LOAD(a_ + j + u), V_LOAD(b_ + j + u))); \
    for (; j < bpad; j += SF_W) \
      V_STORE(d_ + j, VOP(V_LOAD(a_ + j), V_LOAD(b_ + j))); \
  } while (0)

#define SF_UN(VOP) do { \
    const double *a_ = ra; double *d_ = rd; \
    int j = 0; \
    for (; j + SF_STEP <= bpad; j += SF_STEP) \
      for (int u = 0; u < SF_STEP; u += SF_W) \
        V_STORE(d_ + j + u, VOP(V_LOAD(a_ + j + u))); \
    for (; j < bpad; j += SF_W) \
      V_STORE(d_ + j, VOP(V_LOAD(a_ + j))); \
  } while (0)

#define SF_UN_LIBM(FN) do { \
    const double *a_ = ra; double *d_ = rd; \
    for (int j = 0; j < bpad; j++) d_[j] = FN(a_[j]); \
  } while (0)

#define SF_TER_LIBM(FN) do { \
    const double *a_ = ra, *b_ = rb; \
    const double *c_ = regs + (size_t)SF_IARG(3) * B; \
    double *d_ = rd; \
    for (int j = 0; j < bpad; j++) d_[j] = FN(a_[j], b_[j], c_[j]); \
  } while (0)

#define SF_BIN_LIBM(FN) do { \
    const double *a_ = ra, *b_ = rb; double *d_ = rd; \
    for (int j = 0; j < bpad; j++) d_[j] = FN(a_[j], b_[j]); \
  } while (0)

/* R comparison semantics: any NA/NaN operand yields NA, not FALSE. */
#define SF_CMP(MOP) do { \
    const double *a_ = ra, *b_ = rb; double *d_ = rd; \
    for (int j = 0; j < bpad; j += SF_W) { \
      sfv va = V_LOAD(a_ + j), vb = V_LOAD(b_ + j); \
      sfv res = M_TO_NUM(MOP(va, vb)); \
      res = V_SEL(M_OR(M_NAN(va), M_NAN(vb)), V_SET1(SF_NA), res); \
      V_STORE(d_ + j, res); \
    } \
  } while (0)

/* Distribution functions call R's own implementations rather than being
 * rebuilt from erfc(): R's pnorm uses the Cody rational approximations and
 * agrees with erfc only to within a few ulp, so reimplementing it would
 * quietly break bit-identity with base R. */
/* These take location and scale as real operands rather than being rebuilt
 * from a standardised argument. R's dnorm() folds the division by sigma into
 * a different expression order once |z| >= 5 (it splits z to keep z*z exact),
 * so dnorm((x-mu)/s)/s is *not* dnorm(x, mu, s) in the last bits. Passing all
 * three through is both exact and one instruction cheaper. */
static double sf_logis_cdf(double x, double m, double sg) { return plogis(x, m, sg, 1, 0); }
static double sf_logis_q(double p, double m, double sg)   { return qlogis(p, m, sg, 1, 0); }
static double sf_norm_cdf(double x, double m, double sg)  { return pnorm5(x, m, sg, 1, 0); }
static double sf_norm_pdf(double x, double m, double sg)  { return dnorm4(x, m, sg, 0); }
/* %% and %/% are built from the exact IEEE remainder rather than from
 * x1 - floor(x1/x2)*x2. The direct formula rounds twice, and whether those
 * two roundings collapse into one depends on whether the compiler contracted
 * the multiply-subtract into an FMA -- which made the AVX2 kernel disagree
 * with the SSE2 kernel by one ulp on the same data. fmod() is exact by
 * definition, so routing through it removes the dependence on code generation
 * entirely and yields the correctly rounded answer on every backend.
 * R's own large-divisor short circuit is kept, since it is observable. */
static double sf_mod(double x1, double x2) {
  if (x2 == 0.0) return R_NaN;
  if (ISNAN(x1) || ISNAN(x2)) return x1 + x2;   /* preserves R's NA payload */
  if (fabs(x2) * DBL_EPSILON > 1.0 && R_FINITE(x1) && fabs(x1) <= fabs(x2)) {
    return (fabs(x1) == fabs(x2)) ? 0.0
      : (((x1 < 0.0 && x2 > 0.0) || (x2 < 0.0 && x1 > 0.0)) ? x1 + x2 : x1);
  }
  if (!R_FINITE(x1)) return R_NaN;
  double r = fmod(x1, x2);                      /* exact, sign follows x1 */
  if (r == 0.0) return 0.0;
  if ((r < 0.0) != (x2 < 0.0)) {                /* R's result takes x2's sign */
    double s = r + x2;
    return (s == x2) ? 0.0 : s;
  }
  return r;
}

static double sf_idiv(double x1, double x2) {
  double q = x1 / x2;
  if (x2 == 0.0 || ISNAN(x1) || ISNAN(x2) || !R_FINITE(q)) return q;
  if (!R_FINITE(x2)) {
    /* A finite value over an infinite one is an infinitesimal: its floor is 0
     * when the signs agree and -1 when they differ. Base R agrees. */
    if (x1 == 0.0) return 0.0;
    return ((x1 < 0.0) != (x2 < 0.0)) ? -1.0 : 0.0;
  }
  double r = fmod(x1, x2);
  if (r != 0.0 && ((r < 0.0) != (x2 < 0.0))) r += x2;
  /* x1 - r is an exact multiple of x2, so the quotient is an integer and
   * rounding to nearest recovers it exactly. */
  return nearbyint((x1 - r) / x2);
}

/* Integer powers follow R_pow_di's multiply order so results are bit-identical
 * to base R. Non-finite bases are repaired by a scalar pass, which R_pow_di
 * also special-cases. */
SF_TARGET static void SF_NAME(sf_powi_block)(const double *a, double *d,
                                             int bpad, int blen, int n) {
  const int neg = n < 0, m = neg ? -n : n;
  /* Non-finite bases are the only ones that need R_pow's special cases, and
   * they are rare. Detecting them with one compare per vector, and running the
   * scalar repair pass only when one is actually present, matters a lot: the
   * unconditional pass was a branch per element over every block and made
   * x^2 expressions slower than base R on a two-lane machine. */
  sfv nonfinite = V_ZERO();
  const sfv huge = V_SET1(DBL_MAX);

  if (m == 0) {
    for (int j = 0; j < bpad; j++) d[j] = 1.0;      /* x^0 is 1, even for NaN */
    return;
  }
  if (m == 1 && !neg) {
    for (int j = 0; j < bpad; j++) d[j] = a[j];     /* pow(x, 1) == x exactly */
    return;
  }
  if (m == 2 && !neg) {                             /* much the commonest case */
    for (int j = 0; j < bpad; j += SF_W) {
      sfv x = V_LOAD(a + j);
      nonfinite = V_ADD(nonfinite, M_TO_NUM(M_NOT(M_LE(V_ABS(x), huge))));
      V_STORE(d + j, V_MUL(x, x));
    }
  } else {
    for (int j = 0; j < bpad; j += SF_W) {
      sfv x = V_LOAD(a + j);
      nonfinite = V_ADD(nonfinite, M_TO_NUM(M_NOT(M_LE(V_ABS(x), huge))));
      sfv acc = V_SET1(1.0);
      int k = m;
      for (;;) {
        if (k & 1) acc = V_MUL(acc, x);
        k >>= 1;
        if (k) x = V_MUL(x, x); else break;
      }
      if (neg) acc = V_DIV(V_SET1(1.0), acc);
      V_STORE(d + j, acc);
    }
  }

  double nf[SF_W];
  V_STORE(nf, nonfinite);
  double any = 0.0;
  for (int k = 0; k < SF_W; k++) any += nf[k];
  if (any > 0.0)
    for (int j = 0; j < blen; j++)
      if (!R_FINITE(a[j])) d[j] = R_pow(a[j], (double)n);
}

SF_TARGET void SF_NAME(sf_exec_block)(
    const int *code, int ncode, const double *consts,
    const void *const *vars, const int *var_scalar, const int *var_type,
    R_xlen_t start, int blen, double *regs, int nreg, int B, int math_fast)
{
  int bpad = ((blen + SF_W - 1) / SF_W) * SF_W;
  if (bpad > B) bpad = B;
  (void)nreg;

  for (int ip = 0; ip < ncode; ip++) {
    const int op = SF_IARG(0);
    const int a  = SF_IARG(1);
    const int b  = SF_IARG(2);
    const double *ra = regs + (size_t)a * B;
    const double *rb = regs + (size_t)b * B;
    double *rd = regs + (size_t)SF_IARG(4) * B;

    switch (op) {
    case SF_OP_LOADV: {
      double *d_ = rd;
      /* An integer or logical input is converted here rather than up front.
       * This path already copied the block into a register, so the conversion
       * replaces the memcpy instead of adding a pass: no second vector is
       * allocated and the data is touched exactly as often as before.
       * NA_INTEGER is INT_MIN and must become NA_real_, not the double
       * -2147483648, so that is.na() and base R agree. */
      if (var_scalar[a]) {
        double v;
        if (var_type[a] == SF_VT_DBL) {
          v = ((const double *)vars[a])[0];
        } else {
          int iv = ((const int *)vars[a])[0];
          v = (iv == NA_INTEGER) ? NA_REAL : (double)iv;
        }
        for (int j = 0; j < bpad; j++) d_[j] = v;
      } else if (var_type[a] == SF_VT_DBL) {
        memcpy(d_, (const double *)vars[a] + start, (size_t)blen * sizeof(double));
        for (int j = blen; j < bpad; j++) d_[j] = 0.0;
      } else {
        const int *s_ = (const int *)vars[a] + start;
        for (int j = 0; j < blen; j++)
          d_[j] = (s_[j] == NA_INTEGER) ? NA_REAL : (double)s_[j];
        for (int j = blen; j < bpad; j++) d_[j] = 0.0;
      }
      break;
    }
    case SF_OP_LOADC: {
      double v = consts[a]; double *d_ = rd;
      for (int j = 0; j < bpad; j++) d_[j] = v;
      break;
    }
    case SF_OP_ADD:  SF_BIN(V_ADD);  break;
    case SF_OP_SUB:  SF_BIN(V_SUB);  break;
    case SF_OP_MUL:  SF_BIN(V_MUL);  break;
    case SF_OP_DIV:  SF_BIN(V_DIV);  break;
    case SF_OP_NEG:  SF_UN(V_NEGV);  break;
    case SF_OP_MIN2: SF_BIN(V_MIN);  break;
    case SF_OP_MAX2: SF_BIN(V_MAX);  break;
    case SF_OP_SQRT: SF_UN(V_SQRT);  break;
    case SF_OP_ABS:  SF_UN(V_ABS);   break;
    case SF_OP_SIGN: SF_UN(V_SIGN);  break;
    case SF_OP_FLOOR:   SF_UN(V_FLOOR); break;
    case SF_OP_CEILING: SF_UN(V_CEIL);  break;
    case SF_OP_TRUNC:   SF_UN(V_TRUNC); break;
    case SF_OP_ROUND2:  SF_BIN_LIBM(fround); break;
    case SF_OP_FMA3: {
      /* a*b + c, emitted by the R-side compiler when it recognises the
       * pattern; one rounding instead of two, and one instruction instead
       * of two block passes. */
      const double *rc = regs + (size_t)SF_IARG(3) * B;
      const double *a_ = ra, *b_ = rb, *c_ = rc; double *d_ = rd;
      for (int j = 0; j < bpad; j += SF_W)
        V_STORE(d_ + j, V_FMA(V_LOAD(a_ + j), V_LOAD(b_ + j), V_LOAD(c_ + j)));
      break;
    }
    case SF_OP_EXP:
      if (math_fast & SF_FAST_EXP) SF_UN(V_EXP); else SF_UN_LIBM(exp);
      break;
    case SF_OP_LOG:
      if (math_fast & SF_FAST_LOG) SF_UN(V_LOG); else SF_UN_LIBM(log);
      break;
    case SF_OP_LOG1P: SF_UN_LIBM(log1p); break;
    case SF_OP_EXPM1: SF_UN_LIBM(expm1); break;
    case SF_OP_LOG2:  SF_UN_LIBM(log2);  break;
    case SF_OP_LOG10: SF_UN_LIBM(log10); break;
    case SF_OP_SIN:   SF_UN_LIBM(sin);   break;
    case SF_OP_COS:   SF_UN_LIBM(cos);   break;
    case SF_OP_TAN:   SF_UN_LIBM(tan);   break;
    case SF_OP_ASIN:  SF_UN_LIBM(asin);  break;
    case SF_OP_ACOS:  SF_UN_LIBM(acos);  break;
    case SF_OP_ATAN:  SF_UN_LIBM(atan);  break;
    case SF_OP_SINH:  SF_UN_LIBM(sinh);  break;
    case SF_OP_COSH:  SF_UN_LIBM(cosh);  break;
    case SF_OP_TANH:  SF_UN_LIBM(tanh);  break;
    case SF_OP_LGAMMA:  SF_UN_LIBM(lgammafn); break;
    case SF_OP_GAMMA:   SF_UN_LIBM(gammafn);  break;
    case SF_OP_DIGAMMA: SF_UN_LIBM(digamma);  break;
    case SF_OP_ERF:     SF_UN_LIBM(erf);      break;
    case SF_OP_ERFC:    SF_UN_LIBM(erfc);     break;
    case SF_OP_PLOGIS:  SF_TER_LIBM(sf_logis_cdf); break;
    case SF_OP_QLOGIS:  SF_TER_LIBM(sf_logis_q);   break;
    case SF_OP_PNORM:   SF_TER_LIBM(sf_norm_cdf);  break;
    case SF_OP_DNORM:   SF_TER_LIBM(sf_norm_pdf);  break;
    case SF_OP_ATAN2:   SF_BIN_LIBM(atan2);   break;
    case SF_OP_POW:     SF_BIN_LIBM(R_pow);   break;
    case SF_OP_IDIV:    SF_BIN_LIBM(sf_idiv); break;
    case SF_OP_MOD:     SF_BIN_LIBM(sf_mod);  break;
    case SF_OP_POWI:
      SF_NAME(sf_powi_block)(ra, rd, bpad, blen, (int)consts[b]);
      break;
    case SF_OP_LT: SF_CMP(M_LT); break;
    case SF_OP_LE: SF_CMP(M_LE); break;
    case SF_OP_GT: SF_CMP(M_GT); break;
    case SF_OP_GE: SF_CMP(M_GE); break;
    case SF_OP_EQ: SF_CMP(M_EQ); break;
    case SF_OP_NE: SF_CMP(M_NEQ); break;
    case SF_OP_AND: SF_BIN(V_LAND); break;
    case SF_OP_OR:  SF_BIN(V_LOR);  break;
    case SF_OP_NOT: SF_UN(V_LNOT);  break;
    case SF_OP_ISNA: {
      const double *a_ = ra; double *d_ = rd;
      for (int j = 0; j < bpad; j += SF_W)
        V_STORE(d_ + j, M_TO_NUM(M_NAN(V_LOAD(a_ + j))));
      break;
    }
    case SF_OP_ISFINITE: {
      const double *a_ = ra; double *d_ = rd;
      for (int j = 0; j < bpad; j++) d_[j] = R_FINITE(a_[j]) ? 1.0 : 0.0;
      break;
    }
    case SF_OP_ISINF: {
      const double *a_ = ra; double *d_ = rd;
      for (int j = 0; j < bpad; j++)
        d_[j] = (a_[j] == R_PosInf || a_[j] == R_NegInf) ? 1.0 : 0.0;
      break;
    }
    case SF_OP_SELECT: {
      const double *m_ = ra, *x_ = rb;
      const double *y_ = regs + (size_t)SF_IARG(3) * B;
      double *d_ = rd;
      for (int j = 0; j < bpad; j += SF_W) {
        sfv vm = V_LOAD(m_ + j);
        sfv res = V_SEL(M_NOT(M_EQ(vm, V_ZERO())), V_LOAD(x_ + j), V_LOAD(y_ + j));
        /* R's ifelse() yields NA wherever the test is NA or NaN, because the
         * comparison that produced the test already collapsed to NA. */
        res = V_SEL(M_NAN(vm), V_SET1(SF_NA), res);
        V_STORE(d_ + j, res);
      }
      break;
    }
    default: {
      /* Unreachable: the R compiler only emits codes from sf_op_table().
       * Fill with NA so any future bug surfaces as data, not a crash. */
      double *d_ = rd;
      for (int j = 0; j < bpad; j++) d_[j] = R_NaN;
      break;
    }
    }
  }
}

/* Reduce one finished block.
 * out[0] = value, out[1] = auxiliary count:
 *   SUM/MEAN  -> number of elements contributing
 *   MAX/MIN   -> 1 if an NA was seen and na.rm was FALSE
 *   ANY/ALL   -> number of NAs; out[0] counts TRUEs (ANY) or FALSEs (ALL)
 */
SF_TARGET void SF_NAME(sf_reduce_block)(const double *v, int blen, int red,
                                        int na_rm, double *out) {
  int vlen = (blen / SF_W) * SF_W;
  switch (red) {
  case SF_RED_SUM:
  case SF_RED_MEAN: {
    /* One accumulator, deliberately. Splitting this into four independent
     * chains was tried and measured: it is a 0.81x-0.87x regression on
     * AVX-512. The chains themselves are roughly free (0.99x), but making
     * them correct is not -- reordering the additions breaks base R's
     * index-order NA/NaN propagation, and the per-vector scan needed to
     * detect a block that must fall back costs 0.69x-0.89x on its own. The
     * single chain may well be latency-bound on a narrow vector unit with
     * memory bandwidth to spare, but the fix for that has to be cheaper than
     * this one before it is worth having. */
    sfv acc = V_ZERO(), cnt = V_ZERO();
    for (int j = 0; j < vlen; j += SF_W) {
      sfv x = V_LOAD(v + j);
      if (na_rm) {
        sfm ok = M_NOT(M_NAN(x));
        acc = V_ADD(acc, V_SEL(ok, x, V_ZERO()));
        cnt = V_ADD(cnt, M_TO_NUM(ok));
      } else {
        acc = V_ADD(acc, x);
      }
    }
    double buf[SF_W]; V_STORE(buf, acc);
    long double s = 0.0L;
    for (int k = 0; k < SF_W; k++) s += (long double)buf[k];
    double n = 0.0;
    if (na_rm) {
      double cb[SF_W]; V_STORE(cb, cnt);
      for (int k = 0; k < SF_W; k++) n += cb[k];
    }
    for (int j = vlen; j < blen; j++) {
      if (na_rm && ISNAN(v[j])) continue;
      s += (long double)v[j];
      if (na_rm) n += 1.0;
    }
    if (!na_rm) n = (double)blen;
    out[0] = (double)s; out[1] = n;
    break;
  }
  case SF_RED_PROD: {
    long double p = 1.0L; double n = 0.0;
    for (int j = 0; j < blen; j++) {
      if (na_rm && ISNAN(v[j])) continue;
      p *= (long double)v[j]; n += 1.0;
    }
    out[0] = (double)p; out[1] = n;
    break;
  }
  case SF_RED_MAX:
  case SF_RED_MIN: {
    const int isMax = (red == SF_RED_MAX);
    const double ident = isMax ? R_NegInf : R_PosInf;
    sfv acc = V_SET1(ident), nanacc = V_ZERO();
    for (int j = 0; j < vlen; j += SF_W) {
      sfv x = V_LOAD(v + j);
      sfm bad = M_NAN(x);
      nanacc = V_ADD(nanacc, M_TO_NUM(bad));
      x = V_SEL(bad, V_SET1(ident), x);
      acc = isMax ? V_MAXRAW(acc, x) : V_MINRAW(acc, x);
    }
    double buf[SF_W]; V_STORE(buf, acc);
    double r = ident;
    for (int k = 0; k < SF_W; k++)
      r = isMax ? (buf[k] > r ? buf[k] : r) : (buf[k] < r ? buf[k] : r);
    double nb[SF_W]; V_STORE(nb, nanacc);
    double nnan = 0.0;
    for (int k = 0; k < SF_W; k++) nnan += nb[k];
    for (int j = vlen; j < blen; j++) {
      double x = v[j];
      if (ISNAN(x)) { nnan += 1.0; continue; }
      r = isMax ? (x > r ? x : r) : (x < r ? x : r);
    }
    out[0] = r;
    out[1] = (!na_rm && nnan > 0.0) ? 1.0 : 0.0;
    break;
  }
  case SF_RED_ANY:
  case SF_RED_ALL: {
    double nhit = 0.0, nn = 0.0;
    for (int j = 0; j < blen; j++) {
      double x = v[j];
      if (ISNAN(x)) { nn += 1.0; continue; }
      if (red == SF_RED_ANY) { if (x != 0.0) nhit += 1.0; }
      else                   { if (x == 0.0) nhit += 1.0; }
    }
    out[0] = nhit; out[1] = nn;
    break;
  }
  case SF_RED_NOBS: {
    double n = 0.0;
    for (int j = 0; j < blen; j++) if (!ISNAN(v[j])) n += 1.0;
    out[0] = n; out[1] = 0.0;
    break;
  }
  default:
    out[0] = R_NaN; out[1] = 0.0;
  }
}

/* Apply stencil weights across a contiguous buffer.
 * v holds len + nw - 1 values; out[i] = sum_k w[k] * v[i + k].
 * Vectorised over i (the output index) rather than over k, so a short
 * window still gets full lane utilisation. */
SF_TARGET void SF_NAME(sf_stencil_apply)(const double *v, R_xlen_t len,
                                         const double *w, int nw, double *out) {
  R_xlen_t i = 0;
  for (; i + SF_W <= len; i += SF_W) {
    sfv acc = V_ZERO();
    for (int k = 0; k < nw; k++)
      acc = V_FMA(V_SET1(w[k]), V_LOAD(v + i + k), acc);
    V_STORE(out + i, acc);
  }
  for (; i < len; i++) {
    double acc = 0.0;
    for (int k = 0; k < nw; k++) acc += w[k] * v[i + k];
    out[i] = acc;
  }
}
