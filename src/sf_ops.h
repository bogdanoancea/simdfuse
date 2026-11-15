/* sf_ops.h -- the single authoritative opcode list.
 * The enum and the name table are generated from one X-macro so the C engine
 * and the R compiler can never drift apart: R builds its mapping at load time
 * from sf_op_table(), which is generated from this list.
 */
#ifndef SF_OPS_H
#define SF_OPS_H

#define SF_OP_LIST \
  X(LOADV) X(LOADC) \
  X(ADD) X(SUB) X(MUL) X(DIV) X(NEG) X(IDIV) X(MOD) \
  X(SQRT) X(EXP) X(LOG) X(ABS) X(POWI) X(POW) \
  X(MIN2) X(MAX2) \
  X(LT) X(LE) X(GT) X(GE) X(EQ) X(NE) \
  X(AND) X(OR) X(NOT) X(SELECT) \
  X(LOG1P) X(EXPM1) X(LOG2) X(LOG10) \
  X(SIN) X(COS) X(TAN) X(ASIN) X(ACOS) X(ATAN) X(ATAN2) \
  X(SINH) X(COSH) X(TANH) \
  X(FLOOR) X(CEILING) X(TRUNC) X(SIGN) \
  X(LGAMMA) X(GAMMA) X(DIGAMMA) X(ERF) X(ERFC) \
  X(ISNA) X(ISFINITE) X(ISINF) \
  X(PLOGIS) X(QLOGIS) X(PNORM) X(DNORM) X(FMA3) X(ROUND2)

enum {
#define X(nm) SF_OP_##nm,
  SF_OP_LIST
#undef X
  SF_OP_COUNT
};

/* reduction kinds */
enum {
  SF_RED_NONE = 0,
  SF_RED_SUM, SF_RED_MEAN, SF_RED_MAX, SF_RED_MIN, SF_RED_PROD,
  SF_RED_ANY, SF_RED_ALL, SF_RED_NOBS,
  SF_RED_COUNT
};

/* Fast-math is a bitmask rather than a flag: whether a vectorised polynomial
 * beats the platform libm depends on the lane width *and* on how good that
 * libm is, and those pull in different directions on different machines. A
 * two-lane NEON exp() beat Apple's libm by 1.7x while a two-lane SSE2 exp()
 * lost to glibc by 0.9x, so the choice has to be made per function and per
 * machine. sf_tune() measures it. */
/* How an input vector is stored. R's logical and integer vectors are both
 * int arrays and share NA_INTEGER == INT_MIN as their missing value, so they
 * convert identically and need no separate code. */
/* Upper bound on registers in one program, so the kernel's pointer table can
 * live on the stack. The compiler's use-counted allocator keeps the real
 * number small -- a seventeen-instruction log-likelihood uses six. */
#define SF_MAXREG 256

#define SF_VT_DBL 0
#define SF_VT_INT 1

#define SF_FAST_EXP 1
#define SF_FAST_LOG 2

/* one instruction = 5 ints: op, a, b, c, dst */
#define SF_INSTR_WIDTH 5

#endif /* SF_OPS_H */
