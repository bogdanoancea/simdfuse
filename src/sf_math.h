/* sf_math.h -- vectorized elementary functions for simdfuse's "fast" math mode.
 *
 * Include after sf_simd.h. Provides V_EXP and V_LOG as SIMD polynomial
 * evaluations. These are used only when the user selects math = "fast";
 * the default "accurate" mode routes every elementary function through the
 * platform libm per lane, which reproduces base R bit for bit.
 *
 * Accuracy is not asserted here -- it is measured against base R by
 * tests/testthat/test-accuracy.R and reported in the package vignette.
 */
#ifndef SF_MATH_H
#define SF_MATH_H

/* ---- exponent / mantissa extraction -------------------------------------- */

#if SF_ISA == SF_ISA_AVX512
#define V_GETEXP(x)  _mm512_getexp_pd(x)
#define V_GETMANT(x) _mm512_getmant_pd((x), _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)

#elif SF_ISA == SF_ISA_AVX2
SF_TARGET static __inline__ __attribute__((always_inline))
__m256d sf_getexp(__m256d x) {
  __m256i u = _mm256_castpd_si256(x);
  __m256i e = _mm256_srli_epi64(_mm256_slli_epi64(u, 1), 53);
  __m256d ed = _mm256_sub_pd(
      _mm256_castsi256_pd(_mm256_or_si256(e, _mm256_set1_epi64x(0x4330000000000000LL))),
      _mm256_set1_pd(4503599627370496.0));
  return _mm256_sub_pd(ed, _mm256_set1_pd(1023.0));
}
SF_TARGET static __inline__ __attribute__((always_inline))
__m256d sf_getmant(__m256d x) {
  __m256i u = _mm256_castpd_si256(x);
  return _mm256_castsi256_pd(_mm256_or_si256(
      _mm256_and_si256(u, _mm256_set1_epi64x(0x000FFFFFFFFFFFFFLL)),
      _mm256_set1_epi64x(0x3FF0000000000000LL)));
}
#define V_GETEXP(x)  sf_getexp(x)
#define V_GETMANT(x) sf_getmant(x)

#elif SF_ISA == SF_ISA_SSE2
SF_TARGET static __inline__ __attribute__((always_inline))
__m128d sf_getexp(__m128d x) {
  __m128i u = _mm_castpd_si128(x);
  __m128i e = _mm_srli_epi64(_mm_slli_epi64(u, 1), 53);
  __m128d ed = _mm_sub_pd(
      _mm_castsi128_pd(_mm_or_si128(e, _mm_set1_epi64x(0x4330000000000000LL))),
      _mm_set1_pd(4503599627370496.0));
  return _mm_sub_pd(ed, _mm_set1_pd(1023.0));
}
SF_TARGET static __inline__ __attribute__((always_inline))
__m128d sf_getmant(__m128d x) {
  __m128i u = _mm_castpd_si128(x);
  return _mm_castsi128_pd(_mm_or_si128(
      _mm_and_si128(u, _mm_set1_epi64x(0x000FFFFFFFFFFFFFLL)),
      _mm_set1_epi64x(0x3FF0000000000000LL)));
}
#define V_GETEXP(x)  sf_getexp(x)
#define V_GETMANT(x) sf_getmant(x)

#elif SF_ISA == SF_ISA_NEON
static __inline__ __attribute__((always_inline))
float64x2_t sf_getexp(float64x2_t x) {
  uint64x2_t u = vreinterpretq_u64_f64(x);
  uint64x2_t e = vshrq_n_u64(vshlq_n_u64(u, 1), 53);
  return vsubq_f64(vcvtq_f64_s64(vreinterpretq_s64_u64(e)), vdupq_n_f64(1023.0));
}
static __inline__ __attribute__((always_inline))
float64x2_t sf_getmant(float64x2_t x) {
  uint64x2_t u = vreinterpretq_u64_f64(x);
  u = vorrq_u64(vandq_u64(u, vdupq_n_u64(0x000FFFFFFFFFFFFFULL)),
                vdupq_n_u64(0x3FF0000000000000ULL));
  return vreinterpretq_f64_u64(u);
}
#define V_GETEXP(x)  sf_getexp(x)
#define V_GETMANT(x) sf_getmant(x)

#else /* scalar */
static __inline__ __attribute__((always_inline)) double sf_getexp(double x) {
  int e; (void)frexp(x, &e); return (double)(e - 1);
}
static __inline__ __attribute__((always_inline)) double sf_getmant(double x) {
  int e; double m = frexp(x, &e); return m * 2.0;   /* [1,2) */
}
#define V_GETEXP(x)  sf_getexp(x)
#define V_GETMANT(x) sf_getmant(x)
#endif

/* ---- constants ----------------------------------------------------------- */
#define SF_LOG2E   1.4426950408889634074
#define SF_LN2_HI  0.693147180369123816490
#define SF_LN2_LO  1.90821492927058770002e-10
#define SF_SQRT2   1.4142135623730951
#define SF_EXP_HI  709.782712893383996843
#define SF_EXP_LO  -745.133219101941108420
#define SF_DBL_MIN_NORM 2.2250738585072014e-308

/* ---- exp ----------------------------------------------------------------
 * Cody-Waite range reduction to |r| <= ln2/2, then a degree-13 Taylor
 * expansion of exp(r). The truncation term is r^14/14! < 5e-18 over the
 * reduced range, i.e. below one ulp of the result.
 * The 2^k scaling is split in half so that the full double range is
 * representable without building a denormal power of two.
 */
SF_TARGET static __inline__ __attribute__((always_inline))
sfv SF_NAME(sf_vexp)(sfv x) {
  sfv xc = V_MINRAW(V_MAXRAW(x, V_SET1(SF_EXP_LO)), V_SET1(SF_EXP_HI));
  sfv kd = V_RINT(V_MUL(xc, V_SET1(SF_LOG2E)));
  sfv r  = V_SUB(xc, V_MUL(kd, V_SET1(SF_LN2_HI)));
  r      = V_SUB(r,  V_MUL(kd, V_SET1(SF_LN2_LO)));

  sfv p = V_SET1(1.0 / 6227020800.0);                  /* 1/13! */
  p = V_FMA(p, r, V_SET1(1.0 / 479001600.0));          /* 1/12! */
  p = V_FMA(p, r, V_SET1(1.0 / 39916800.0));           /* 1/11! */
  p = V_FMA(p, r, V_SET1(1.0 / 3628800.0));            /* 1/10! */
  p = V_FMA(p, r, V_SET1(1.0 / 362880.0));             /* 1/9!  */
  p = V_FMA(p, r, V_SET1(1.0 / 40320.0));              /* 1/8!  */
  p = V_FMA(p, r, V_SET1(1.0 / 5040.0));               /* 1/7!  */
  p = V_FMA(p, r, V_SET1(1.0 / 720.0));                /* 1/6!  */
  p = V_FMA(p, r, V_SET1(1.0 / 120.0));                /* 1/5!  */
  p = V_FMA(p, r, V_SET1(1.0 / 24.0));                 /* 1/4!  */
  p = V_FMA(p, r, V_SET1(1.0 / 6.0));                  /* 1/3!  */
  p = V_FMA(p, r, V_SET1(0.5));                        /* 1/2!  */
  p = V_FMA(p, r, V_SET1(1.0));
  p = V_FMA(p, r, V_SET1(1.0));

  sfv kh = V_TRUNC(V_MUL(kd, V_SET1(0.5)));
  sfv res = V_SCALE2(V_SCALE2(p, kh), V_SUB(kd, kh));

  res = V_SEL(M_GT(x, V_SET1(SF_EXP_HI)), V_SET1(1.0 / 0.0), res);
  res = V_SEL(M_LT(x, V_SET1(SF_EXP_LO)), V_ZERO(), res);
  res = V_SEL(M_NAN(x), V_SET1(SF_NAN), res);
  return res;
}

/* ---- log ----------------------------------------------------------------
 * x = m * 2^e with m in [sqrt(1/2), sqrt(2)); log(m) via the atanh series
 * log((1+s)/(1-s)) = 2(s + s^3/3 + s^5/5 + ...) with s = (m-1)/(m+1),
 * |s| <= 0.1716 so s^2 <= 0.0295 and eleven terms put the truncation
 * below one ulp. Subnormal inputs are pre-scaled by 2^54.
 */
SF_TARGET static __inline__ __attribute__((always_inline))
sfv SF_NAME(sf_vlog)(sfv x) {
  sfm sub = M_LT(x, V_SET1(SF_DBL_MIN_NORM));
  sfv xs  = V_SEL(sub, V_MUL(x, V_SET1(18014398509481984.0)), x);  /* 2^54 */
  sfv corr = V_SEL(sub, V_SET1(-54.0), V_ZERO());

  sfv e = V_ADD(V_GETEXP(xs), corr);
  sfv m = V_GETMANT(xs);
  sfm big = M_GT(m, V_SET1(SF_SQRT2));
  m = V_SEL(big, V_MUL(m, V_SET1(0.5)), m);
  e = V_SEL(big, V_ADD(e, V_SET1(1.0)), e);

  sfv s  = V_DIV(V_SUB(m, V_SET1(1.0)), V_ADD(m, V_SET1(1.0)));
  sfv s2 = V_MUL(s, s);

  sfv p = V_SET1(1.0 / 23.0);
  p = V_FMA(p, s2, V_SET1(1.0 / 21.0));
  p = V_FMA(p, s2, V_SET1(1.0 / 19.0));
  p = V_FMA(p, s2, V_SET1(1.0 / 17.0));
  p = V_FMA(p, s2, V_SET1(1.0 / 15.0));
  p = V_FMA(p, s2, V_SET1(1.0 / 13.0));
  p = V_FMA(p, s2, V_SET1(1.0 / 11.0));
  p = V_FMA(p, s2, V_SET1(1.0 / 9.0));
  p = V_FMA(p, s2, V_SET1(1.0 / 7.0));
  p = V_FMA(p, s2, V_SET1(1.0 / 5.0));
  p = V_FMA(p, s2, V_SET1(1.0 / 3.0));
  p = V_FMA(p, s2, V_SET1(1.0));
  sfv logm = V_MUL(V_MUL(s, V_SET1(2.0)), p);

  sfv res = V_FMA(e, V_SET1(SF_LN2_HI), V_FMA(e, V_SET1(SF_LN2_LO), logm));

  res = V_SEL(M_EQ(x, V_ZERO()), V_SET1(-1.0 / 0.0), res);
  res = V_SEL(M_LT(x, V_ZERO()), V_SET1(SF_NAN), res);
  res = V_SEL(M_NAN(x), V_SET1(SF_NAN), res);
  res = V_SEL(M_EQ(x, V_SET1(1.0 / 0.0)), V_SET1(1.0 / 0.0), res);
  return res;
}

#define V_EXP(x) SF_NAME(sf_vexp)(x)
#define V_LOG(x) SF_NAME(sf_vlog)(x)

#endif /* SF_MATH_H */
