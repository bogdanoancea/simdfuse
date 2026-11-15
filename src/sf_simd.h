/* sf_simd.h -- ISA abstraction layer for simdfuse.
 *
 * This header is included once per translation unit, with SF_ISA defined to
 * select the instruction set. Every operation is a macro, so no call boundary
 * sits between the kernel body and the intrinsic; this matters because the
 * kernel functions carry a function-level `target` attribute and GCC/Clang
 * refuse to inline across mismatched target boundaries.
 *
 * Contract for every ISA block below:
 *   SF_TARGET        function attribute enabling the ISA (may be empty)
 *   SF_NAME(x)       appends the ISA suffix to a symbol name
 *   SF_W             lanes per vector (doubles)
 *   sfv              vector type
 *   sfm              mask type
 *   V_*              vector operations
 *   M_*              mask operations
 */
#ifndef SF_SIMD_H
#define SF_SIMD_H

#include <math.h>
#include <stdint.h>
#include <string.h>

#define SF_ISA_SCALAR 0
#define SF_ISA_SSE2   1
#define SF_ISA_AVX2   2
#define SF_ISA_AVX512 3
#define SF_ISA_NEON   4

#ifndef SF_ISA
#error "SF_ISA must be defined before including sf_simd.h"
#endif

/* ------------------------------------------------------------------ AVX-512 */
#if SF_ISA == SF_ISA_AVX512
#include <immintrin.h>
#define SF_TARGET __attribute__((target("avx512f")))
#define SF_NAME(x) x##_avx512
#define SF_ISANAME "avx512f"
#define SF_W 8
typedef __m512d sfv;
typedef __mmask8 sfm;
#define V_SET1(a)      _mm512_set1_pd(a)
#define V_ZERO()       _mm512_setzero_pd()
#define V_LOAD(p)      _mm512_loadu_pd(p)
#define V_STORE(p,v)   _mm512_storeu_pd((p),(v))
#define V_ADD(a,b)     _mm512_add_pd(a,b)
#define V_SUB(a,b)     _mm512_sub_pd(a,b)
#define V_MUL(a,b)     _mm512_mul_pd(a,b)
#define V_DIV(a,b)     _mm512_div_pd(a,b)
#define V_FMA(a,b,c)   _mm512_fmadd_pd(a,b,c)
#define V_SQRT(a)      _mm512_sqrt_pd(a)
#define V_ABS(a)       _mm512_abs_pd(a)
#define V_MINRAW(a,b)  _mm512_min_pd(a,b)
#define V_MAXRAW(a,b)  _mm512_max_pd(a,b)
#define V_FLOOR(a)     _mm512_roundscale_pd(a,0x01)
#define V_CEIL(a)      _mm512_roundscale_pd(a,0x02)
#define V_TRUNC(a)     _mm512_roundscale_pd(a,0x03)
#define V_RINT(a)      _mm512_roundscale_pd(a,0x00)
#define M_LT(a,b)      _mm512_cmp_pd_mask(a,b,_CMP_LT_OQ)
#define M_LE(a,b)      _mm512_cmp_pd_mask(a,b,_CMP_LE_OQ)
#define M_GT(a,b)      _mm512_cmp_pd_mask(a,b,_CMP_GT_OQ)
#define M_GE(a,b)      _mm512_cmp_pd_mask(a,b,_CMP_GE_OQ)
#define M_EQ(a,b)      _mm512_cmp_pd_mask(a,b,_CMP_EQ_OQ)
#define M_NEQ(a,b)     _mm512_cmp_pd_mask(a,b,_CMP_NEQ_OQ)
#define M_NAN(a)       _mm512_cmp_pd_mask(a,a,_CMP_UNORD_Q)
#define M_AND(a,b)     ((sfm)((a)&(b)))
#define M_OR(a,b)      ((sfm)((a)|(b)))
#define M_NOT(a)       ((sfm)(~(a)))
#define V_SEL(m,a,b)   _mm512_mask_blend_pd(m,b,a)   /* m ? a : b */
#define V_SCALE2(p,kd) _mm512_scalef_pd(p,kd)        /* p * 2^trunc(kd) */

/* --------------------------------------------------------------------- AVX2 */
#elif SF_ISA == SF_ISA_AVX2
#include <immintrin.h>
#define SF_TARGET __attribute__((target("avx2,fma")))
#define SF_NAME(x) x##_avx2
#define SF_ISANAME "avx2+fma"
#define SF_W 4
typedef __m256d sfv;
typedef __m256d sfm;
#define V_SET1(a)      _mm256_set1_pd(a)
#define V_ZERO()       _mm256_setzero_pd()
#define V_LOAD(p)      _mm256_loadu_pd(p)
#define V_STORE(p,v)   _mm256_storeu_pd((p),(v))
#define V_ADD(a,b)     _mm256_add_pd(a,b)
#define V_SUB(a,b)     _mm256_sub_pd(a,b)
#define V_MUL(a,b)     _mm256_mul_pd(a,b)
#define V_DIV(a,b)     _mm256_div_pd(a,b)
#define V_FMA(a,b,c)   _mm256_fmadd_pd(a,b,c)
#define V_SQRT(a)      _mm256_sqrt_pd(a)
#define V_ABS(a)       _mm256_andnot_pd(_mm256_set1_pd(-0.0),(a))
#define V_MINRAW(a,b)  _mm256_min_pd(a,b)
#define V_MAXRAW(a,b)  _mm256_max_pd(a,b)
#define V_FLOOR(a)     _mm256_round_pd(a,_MM_FROUND_TO_NEG_INF|_MM_FROUND_NO_EXC)
#define V_CEIL(a)      _mm256_round_pd(a,_MM_FROUND_TO_POS_INF|_MM_FROUND_NO_EXC)
#define V_TRUNC(a)     _mm256_round_pd(a,_MM_FROUND_TO_ZERO|_MM_FROUND_NO_EXC)
#define V_RINT(a)      _mm256_round_pd(a,_MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC)
#define M_LT(a,b)      _mm256_cmp_pd(a,b,_CMP_LT_OQ)
#define M_LE(a,b)      _mm256_cmp_pd(a,b,_CMP_LE_OQ)
#define M_GT(a,b)      _mm256_cmp_pd(a,b,_CMP_GT_OQ)
#define M_GE(a,b)      _mm256_cmp_pd(a,b,_CMP_GE_OQ)
#define M_EQ(a,b)      _mm256_cmp_pd(a,b,_CMP_EQ_OQ)
#define M_NEQ(a,b)     _mm256_cmp_pd(a,b,_CMP_NEQ_OQ)
#define M_NAN(a)       _mm256_cmp_pd(a,a,_CMP_UNORD_Q)
#define M_AND(a,b)     _mm256_and_pd(a,b)
#define M_OR(a,b)      _mm256_or_pd(a,b)
#define M_NOT(a)       _mm256_xor_pd(a,_mm256_castsi256_pd(_mm256_set1_epi64x(-1)))
#define V_SEL(m,a,b)   _mm256_blendv_pd(b,a,m)
SF_TARGET static __inline__ __attribute__((always_inline))
__m256d sf_scale2_avx2(__m256d p, __m256d kd) {
  __m128i k32 = _mm256_cvtpd_epi32(kd);                  /* 4 x int32 */
  __m256i k64 = _mm256_cvtepi32_epi64(k32);              /* 4 x int64 */
  k64 = _mm256_add_epi64(k64, _mm256_set1_epi64x(1023));
  k64 = _mm256_slli_epi64(k64, 52);
  return _mm256_mul_pd(p, _mm256_castsi256_pd(k64));
}
#define V_SCALE2(p,kd) sf_scale2_avx2(p,kd)

/* --------------------------------------------------------------------- SSE2 */
#elif SF_ISA == SF_ISA_SSE2
#include <emmintrin.h>
#define SF_TARGET __attribute__((target("sse2")))
#define SF_NAME(x) x##_sse2
#define SF_ISANAME "sse2"
#define SF_W 2
typedef __m128d sfv;
typedef __m128d sfm;
#define V_SET1(a)      _mm_set1_pd(a)
#define V_ZERO()       _mm_setzero_pd()
#define V_LOAD(p)      _mm_loadu_pd(p)
#define V_STORE(p,v)   _mm_storeu_pd((p),(v))
#define V_ADD(a,b)     _mm_add_pd(a,b)
#define V_SUB(a,b)     _mm_sub_pd(a,b)
#define V_MUL(a,b)     _mm_mul_pd(a,b)
#define V_DIV(a,b)     _mm_div_pd(a,b)
#define V_FMA(a,b,c)   _mm_add_pd(_mm_mul_pd(a,b),c)
#define V_SQRT(a)      _mm_sqrt_pd(a)
#define V_ABS(a)       _mm_andnot_pd(_mm_set1_pd(-0.0),(a))
#define V_MINRAW(a,b)  _mm_min_pd(a,b)
#define V_MAXRAW(a,b)  _mm_max_pd(a,b)
#define M_LT(a,b)      _mm_cmplt_pd(a,b)
#define M_LE(a,b)      _mm_cmple_pd(a,b)
#define M_GT(a,b)      _mm_cmpgt_pd(a,b)
#define M_GE(a,b)      _mm_cmpge_pd(a,b)
#define M_EQ(a,b)      _mm_cmpeq_pd(a,b)
#define M_NEQ(a,b)     _mm_cmpneq_pd(a,b)
#define M_NAN(a)       _mm_cmpunord_pd(a,a)
#define M_AND(a,b)     _mm_and_pd(a,b)
#define M_OR(a,b)      _mm_or_pd(a,b)
#define M_NOT(a)       _mm_xor_pd(a,_mm_castsi128_pd(_mm_set1_epi32(-1)))
#define V_SEL(m,a,b)   _mm_or_pd(_mm_and_pd(m,a),_mm_andnot_pd(m,b))
/* SSE2 has no round instruction: use the 2^52 magic-number technique. */
SF_TARGET static __inline__ __attribute__((always_inline))
__m128d sf_rint_sse2(__m128d a) {
  const __m128d magic = _mm_set1_pd(4503599627370496.0);   /* 2^52 */
  __m128d s = _mm_and_pd(_mm_set1_pd(-0.0), a);            /* sign */
  __m128d m = _mm_or_pd(magic, s);                         /* +-2^52 */
  __m128d r = _mm_sub_pd(_mm_add_pd(a, m), m);
  /* |a| >= 2^52 is already integral; the trick would lose it. */
  __m128d big = _mm_cmpge_pd(_mm_andnot_pd(_mm_set1_pd(-0.0), a), magic);
  return _mm_or_pd(_mm_and_pd(big, a), _mm_andnot_pd(big, r));
}
SF_TARGET static __inline__ __attribute__((always_inline))
__m128d sf_trunc_sse2(__m128d a) {
  __m128d r = sf_rint_sse2(a);
  __m128d over = _mm_cmpgt_pd(_mm_andnot_pd(_mm_set1_pd(-0.0), r),
                              _mm_andnot_pd(_mm_set1_pd(-0.0), a));
  __m128d one = _mm_or_pd(_mm_set1_pd(1.0), _mm_and_pd(_mm_set1_pd(-0.0), a));
  return _mm_sub_pd(r, _mm_and_pd(over, one));
}
#define V_RINT(a)  sf_rint_sse2(a)
#define V_TRUNC(a) sf_trunc_sse2(a)
SF_TARGET static __inline__ __attribute__((always_inline))
__m128d sf_floor_sse2(__m128d a) {
  __m128d t = sf_trunc_sse2(a);
  __m128d adj = _mm_cmplt_pd(a, t);
  return _mm_sub_pd(t, _mm_and_pd(adj, _mm_set1_pd(1.0)));
}
SF_TARGET static __inline__ __attribute__((always_inline))
__m128d sf_ceil_sse2(__m128d a) {
  __m128d t = sf_trunc_sse2(a);
  __m128d adj = _mm_cmpgt_pd(a, t);
  return _mm_add_pd(t, _mm_and_pd(adj, _mm_set1_pd(1.0)));
}
#define V_FLOOR(a) sf_floor_sse2(a)
#define V_CEIL(a)  sf_ceil_sse2(a)
SF_TARGET static __inline__ __attribute__((always_inline))
__m128d sf_scale2_sse2(__m128d p, __m128d kd) {
  __m128i k32 = _mm_cvtpd_epi32(kd);                       /* [k0,k1,0,0] */
  k32 = _mm_add_epi32(k32, _mm_set1_epi32(1023));
  /* widen to 64-bit lanes shifted left by 32, then by 20 more => <<52 */
  __m128i lo = _mm_unpacklo_epi32(_mm_setzero_si128(), k32);
  lo = _mm_slli_epi64(lo, 20);
  return _mm_mul_pd(p, _mm_castsi128_pd(lo));
}
#define V_SCALE2(p,kd) sf_scale2_sse2(p,kd)

/* --------------------------------------------------------------------- NEON */
#elif SF_ISA == SF_ISA_NEON
#include <arm_neon.h>
#define SF_TARGET
#define SF_NAME(x) x##_neon
#define SF_ISANAME "neon"
#define SF_W 2
typedef float64x2_t sfv;
typedef uint64x2_t  sfm;
#define V_SET1(a)      vdupq_n_f64(a)
#define V_ZERO()       vdupq_n_f64(0.0)
#define V_LOAD(p)      vld1q_f64(p)
#define V_STORE(p,v)   vst1q_f64((p),(v))
#define V_ADD(a,b)     vaddq_f64(a,b)
#define V_SUB(a,b)     vsubq_f64(a,b)
#define V_MUL(a,b)     vmulq_f64(a,b)
#define V_DIV(a,b)     vdivq_f64(a,b)
#define V_FMA(a,b,c)   vfmaq_f64(c,a,b)
#define V_SQRT(a)      vsqrtq_f64(a)
#define V_ABS(a)       vabsq_f64(a)
#define V_MINRAW(a,b)  vminq_f64(a,b)
#define V_MAXRAW(a,b)  vmaxq_f64(a,b)
#define V_FLOOR(a)     vrndmq_f64(a)
#define V_CEIL(a)      vrndpq_f64(a)
#define V_TRUNC(a)     vrndq_f64(a)
#define V_RINT(a)      vrndnq_f64(a)
#define M_LT(a,b)      vcltq_f64(a,b)
#define M_LE(a,b)      vcleq_f64(a,b)
#define M_GT(a,b)      vcgtq_f64(a,b)
#define M_GE(a,b)      vcgeq_f64(a,b)
#define M_EQ(a,b)      vceqq_f64(a,b)
#define M_NEQ(a,b)     vmvnq_u32_wrap(vceqq_f64(a,b))
#define M_NAN(a)       vmvnq_u32_wrap(vceqq_f64(a,a))
#define M_AND(a,b)     vandq_u64(a,b)
#define M_OR(a,b)      vorrq_u64(a,b)
#define M_NOT(a)       vmvnq_u32_wrap(a)
#define V_SEL(m,a,b)   vbslq_f64(m,a,b)
static __inline__ __attribute__((always_inline))
uint64x2_t vmvnq_u32_wrap(uint64x2_t a) {
  return vreinterpretq_u64_u32(vmvnq_u32(vreinterpretq_u32_u64(a)));
}
static __inline__ __attribute__((always_inline))
float64x2_t sf_scale2_neon(float64x2_t p, float64x2_t kd) {
  int64x2_t k = vcvtq_s64_f64(kd);
  k = vaddq_s64(k, vdupq_n_s64(1023));
  k = vshlq_n_s64(k, 52);
  return vmulq_f64(p, vreinterpretq_f64_s64(k));
}
#define V_SCALE2(p,kd) sf_scale2_neon(p,kd)

/* ------------------------------------------------------------------- scalar */
#elif SF_ISA == SF_ISA_SCALAR
#define SF_TARGET
#define SF_NAME(x) x##_scalar
#define SF_ISANAME "scalar"
#define SF_W 1
typedef double sfv;
typedef int    sfm;
#define V_SET1(a)      ((double)(a))
#define V_ZERO()       (0.0)
#define V_LOAD(p)      (*(const double*)(p))
#define V_STORE(p,v)   (*(double*)(p) = (v))
#define V_ADD(a,b)     ((a)+(b))
#define V_SUB(a,b)     ((a)-(b))
#define V_MUL(a,b)     ((a)*(b))
#define V_DIV(a,b)     ((a)/(b))
#define V_FMA(a,b,c)   ((a)*(b)+(c))
#define V_SQRT(a)      sqrt(a)
#define V_ABS(a)       fabs(a)
#define V_MINRAW(a,b)  ((a)<(b)?(a):(b))
#define V_MAXRAW(a,b)  ((a)>(b)?(a):(b))
#define V_FLOOR(a)     floor(a)
#define V_CEIL(a)      ceil(a)
#define V_TRUNC(a)     trunc(a)
#define V_RINT(a)      nearbyint(a)
#define M_LT(a,b)      ((a)<(b))
#define M_LE(a,b)      ((a)<=(b))
#define M_GT(a,b)      ((a)>(b))
#define M_GE(a,b)      ((a)>=(b))
#define M_EQ(a,b)      ((a)==(b))
#define M_NEQ(a,b)     (!((a)==(b)) && !(isnan(a)||isnan(b)))
#define M_NAN(a)       (isnan(a)!=0)
#define M_AND(a,b)     ((a)&&(b))
#define M_OR(a,b)      ((a)||(b))
#define M_NOT(a)       (!(a))
#define V_SEL(m,a,b)   ((m)?(a):(b))
#define V_SCALE2(p,kd) ldexp((p),(int)(kd))
#else
#error "unknown SF_ISA"
#endif

/* ---------------------------------------------------- derived, ISA-agnostic */

/* R semantics: min/max propagate missingness, unlike the hardware MINPD and
 * MAXPD. R's pmin/pmax do it by *adding* the operands when either is NaN,
 * which preserves the payload bits that distinguish NA_real_ from a plain
 * NaN; substituting a fresh NaN here would silently turn NA into NaN. The
 * operand order matches R's, so the same input wins the payload. */
#define V_MIN(a,b) V_SEL(M_OR(M_NAN(a),M_NAN(b)), V_ADD(a,b), V_MINRAW(a,b))
#define V_MAX(a,b) V_SEL(M_OR(M_NAN(a),M_NAN(b)), V_ADD(a,b), V_MAXRAW(a,b))

#define SF_NAN (0.0/0.0)

/* R distinguishes NA from NaN by the payload of a quiet NaN. Arithmetic that
 * genuinely produces an invalid result (0/0) yields NaN in R too, so SF_NAN is
 * right there. But comparisons, the logical operators and ifelse() yield R's
 * NA, so those use SF_NA. Expanded only inside translation units that have
 * already included R.h. */
#define SF_NA NA_REAL

/* Three-valued logic, matching R: FALSE dominates &, TRUE dominates |. */
#define V_LAND(a,b) \
  V_SEL(M_OR(M_AND(M_EQ(a,V_ZERO()),M_NOT(M_NAN(a))), \
             M_AND(M_EQ(b,V_ZERO()),M_NOT(M_NAN(b)))), \
        V_ZERO(), \
        V_SEL(M_OR(M_NAN(a),M_NAN(b)), V_SET1(SF_NA), V_SET1(1.0)))
#define V_LOR(a,b) \
  V_SEL(M_OR(M_AND(M_NOT(M_EQ(a,V_ZERO())),M_NOT(M_NAN(a))), \
             M_AND(M_NOT(M_EQ(b,V_ZERO())),M_NOT(M_NAN(b)))), \
        V_SET1(1.0), \
        V_SEL(M_OR(M_NAN(a),M_NAN(b)), V_SET1(SF_NA), V_ZERO()))
#define V_LNOT(a) \
  V_SEL(M_NAN(a), V_SET1(SF_NA), V_SEL(M_EQ(a,V_ZERO()), V_SET1(1.0), V_ZERO()))

/* R's sign() returns its argument unchanged when it is NA or NaN, which keeps
 * the two apart; a freshly built NaN would not. */
#define V_SIGN(a) \
  V_SEL(M_NAN(a), (a), \
    V_SEL(M_GT(a,V_ZERO()), V_SET1(1.0), \
      V_SEL(M_LT(a,V_ZERO()), V_SET1(-1.0), V_ZERO())))

#define M_TO_NUM(m) V_SEL(m, V_SET1(1.0), V_ZERO())

#define V_NEGV(a) V_SUB(V_ZERO(), (a))

#endif /* SF_SIMD_H */
