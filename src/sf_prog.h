/* sf_prog.h -- shared program representation and dispatch surface.
 *
 * Every entry point (elementwise evaluation, reductions, scan, stencil,
 * segmented reduce) runs the *same* compiled program through the *same*
 * block executor. That is what lets the skeletons fuse with an arbitrary
 * expression instead of being limited to a plain input vector.
 */
#ifndef SF_PROG_H
#define SF_PROG_H

#include "sf_omp.h"   /* must precede R.h: see sf_omp.h */

#include <R.h>
#include <Rinternals.h>
#include "sf_ops.h"

#define SF_TASK_ELEMS 65536

/* Tasks a worker must receive before it is worth creating.
 *
 * Joining a team is not free. Measured on a twelve-thread Raptor Lake-U
 * laptop, the team costs about 76 us at two threads and 227 us at twelve,
 * against roughly 46 us of work per task. The previous gate allowed one task
 * per thread, so just above the parallel threshold a four-thread team did one
 * task each and lost to the serial path on every expression a benchmark could
 * throw at it. Requiring several tasks per worker gives up a little
 * parallelism at small n in exchange for never being slower than one thread.
 *
 * The answer does not depend on this. Task boundaries are a function of n
 * alone, so the thread count cannot move a single bit of the result. */
#define SF_MIN_TASKS_PER_THREAD 4

/* Serial below this: fewer than two workers' worth of tasks exist. Kept in
 * step with the constant above so the two gates cannot disagree. */
#define SF_MIN_PARALLEL \
  ((R_xlen_t)SF_TASK_ELEMS * SF_MIN_TASKS_PER_THREAD * 2)

typedef struct {
  const int *code;
  int ncode;
  const double *consts;
  int nconst;
  /* Inputs are not coerced to double before evaluation: an integer or
   * logical vector is read in its own storage and converted during the block
   * load, which was already copying into a register anyway. So the pointers
   * are untyped and vtype says how to read each one. */
  const void *const *vars;
  const int *vtype;             /* SF_VT_DBL or SF_VT_INT, per input */
  const int *vscal;
  int nvar;
  int nreg;
  int nprologue;   /* leading instructions that run once per call, not per block */
  int B;
  int fast;
  int final_reg;
  R_xlen_t n;
} sf_prog;

/* Build and validate a program from the standard SEXP argument block.
 * Raises an R error on any malformed input, so callers may assume validity. */
void sf_prog_init(sf_prog *P, SEXP s_code, SEXP s_consts, SEXP s_vars,
                  SEXP s_n, SEXP s_nreg, SEXP s_nprologue, SEXP s_block,
                  SEXP s_fast);

/* Evaluate the block-invariant prologue into a thread's register file. Must be
 * called once per thread before that thread's first sf_prog_block(). */
void sf_prog_prologue(const sf_prog *P, double *regs);

/* Execute the program over [start, start+blen) into a thread's register file. */
void sf_prog_block(const sf_prog *P, double *regs, R_xlen_t start, int blen);

/* As sf_prog_block, but the final result of the block is written straight into
 * dest when the block needs no padding, saving the caller's copy. Returns 1 if
 * it did so, 0 if the caller must still copy out of sf_prog_final(). */
void sf_prog_block_to(const sf_prog *P, double *regs, R_xlen_t start, int blen,
                      double *dest);
int  sf_prog_dest_used(const sf_prog *P, int blen);

/* Reduce a finished block; see sf_kernels.h for the out[] contract. */
void sf_prog_reduce(const double *v, int blen, int red, int na_rm, double *out);

/* Pointer to the register holding the program's result. */
double *sf_prog_final(const sf_prog *P, double *regs);

/* Allocate an aligned per-thread register file via R_alloc. */
double *sf_prog_regfile(const sf_prog *P, int nthreads);

/* Resolve a usable thread count for a given workload. */
int sf_prog_threads(R_xlen_t n, R_xlen_t ntask, int requested);

/* Apply stencil weights to a contiguous buffer. `v` must hold len+nw-1
 * values (the window halo is the caller's responsibility); `out` receives
 * len results. ISA-specialised in sf_kernels.h. */
void sf_prog_stencil(const double *v, R_xlen_t len, const double *w, int nw,
                     double *out);

int sf_isa_lanes(void);

/* ---- shared reduction accumulator -------------------------------------
 * Used by both the fused reductions and the segmented reduce so the two
 * cannot disagree about NA handling or combine order.
 */
typedef struct { long double v; double aux, aux2; } sf_acc;

static inline void sf_acc_init(sf_acc *A, int red) {
  A->aux = 0.0; A->aux2 = 0.0;
  switch (red) {
  case SF_RED_PROD: A->v = 1.0L; break;
  case SF_RED_MAX:  A->v = (long double)R_NegInf; break;
  case SF_RED_MIN:  A->v = (long double)R_PosInf; break;
  default:          A->v = 0.0L; break;
  }
}

static inline void sf_acc_push(sf_acc *A, int red, const double *blk) {
  switch (red) {
  case SF_RED_SUM: case SF_RED_MEAN:
    A->v += (long double)blk[0]; A->aux += blk[1]; break;
  case SF_RED_NOBS:
    A->v += (long double)blk[0]; break;
  case SF_RED_PROD:
    A->v *= (long double)blk[0]; A->aux += blk[1]; break;
  case SF_RED_MAX:
    if ((long double)blk[0] > A->v) { A->v = (long double)blk[0]; }
    if (blk[1] != 0.0) { A->aux = 1.0; }
    break;
  case SF_RED_MIN:
    if ((long double)blk[0] < A->v) { A->v = (long double)blk[0]; }
    if (blk[1] != 0.0) { A->aux = 1.0; }
    break;
  case SF_RED_ANY: case SF_RED_ALL:
    A->aux += blk[0]; A->aux2 += blk[1]; break;
  default: break;
  }
}

/* Final scalar value of an accumulator, with R's NA conventions applied. */
static inline double sf_acc_value(const sf_acc *A, int red, int na_rm) {
  switch (red) {
  case SF_RED_SUM: case SF_RED_PROD: case SF_RED_NOBS:
    return (double)A->v;
  case SF_RED_MEAN:
    return (A->aux <= 0.0) ? R_NaN : (double)(A->v / (long double)A->aux);
  case SF_RED_MAX: case SF_RED_MIN:
    return (A->aux != 0.0) ? NA_REAL : (double)A->v;
  case SF_RED_ANY:
    if (A->aux > 0.0) return 1.0;
    if (A->aux2 > 0.0 && !na_rm) return NA_REAL;
    return 0.0;
  case SF_RED_ALL:
    if (A->aux > 0.0) return 0.0;
    if (A->aux2 > 0.0 && !na_rm) return NA_REAL;
    return 1.0;
  default: return R_NaN;
  }
}

#endif /* SF_PROG_H */
