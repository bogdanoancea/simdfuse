/* sf_dispatch.c -- instruction-set selection, program validation, registration.
 *
 * Dispatch happens once per call rather than once per element: a function
 * pointer pair is resolved at package load from the CPU's reported features,
 * and the whole evaluation runs through it. Users can pin a narrower backend
 * (sf_options(isa=)) so the test suite can sweep every code path on one machine.
 */

#include "sf_omp.h"   /* first: must precede any R header */
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <string.h>
#include <stdint.h>
#include "sf_prog.h"


typedef void (*sf_exec_fn)(const int*, int, const double*, const void *const*,
                           const int*, const int*, R_xlen_t, int, double*,
                           int, int, int, double*, int);
typedef void (*sf_reduce_fn)(const double*, int, int, int, double*);
typedef void (*sf_stencil_fn)(const double*, R_xlen_t, const double*, int, double*);

#define SF_DECL(suf) \
  extern void sf_exec_block_##suf(const int*, int, const double*, \
      const void *const*, const int*, const int*, R_xlen_t, int, double*, \
      int, int, int, double*, int); \
  extern void sf_reduce_block_##suf(const double*, int, int, int, double*); \
  extern void sf_stencil_apply_##suf(const double*, R_xlen_t, const double*, int, double*);

SF_DECL(scalar)
#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
SF_DECL(sse2)
SF_DECL(avx2)
SF_DECL(avx512)
#endif
#if defined(__aarch64__)
SF_DECL(neon)
#endif

static sf_exec_fn    g_exec    = NULL;
static sf_reduce_fn  g_reduce  = NULL;
static sf_stencil_fn g_stencil = NULL;
static const char  *g_isa    = "scalar";
static int          g_lanes  = 1;
static int          g_forced = 0;

#define SF_PICK(suf, nm, w) do { \
    g_exec = sf_exec_block_##suf; g_reduce = sf_reduce_block_##suf; \
    g_stencil = sf_stencil_apply_##suf; \
    g_isa = nm; g_lanes = w; } while (0)

static void sf_select_isa(const char *want) {
  int force = (want && strcmp(want, "auto") != 0);
  int unavailable = 0;
  g_forced = 0;
#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
  const int has512 = __builtin_cpu_supports("avx512f");
  const int has2   = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
  const int hasse  = __builtin_cpu_supports("sse2");
  if (force) {
    if (!strcmp(want, "avx512") && has512) { SF_PICK(avx512, "avx512f", 8); g_forced = 1; return; }
    if (!strcmp(want, "avx2")   && has2)   { SF_PICK(avx2, "avx2+fma", 4);  g_forced = 1; return; }
    if (!strcmp(want, "sse2")   && hasse)  { SF_PICK(sse2, "sse2", 2);      g_forced = 1; return; }
    if (!strcmp(want, "scalar"))           { SF_PICK(scalar, "scalar", 1);  g_forced = 1; return; }
    unavailable = 1;
  }
  if (has512)      SF_PICK(avx512, "avx512f", 8);
  else if (has2)   SF_PICK(avx2, "avx2+fma", 4);
  else if (hasse)  SF_PICK(sse2, "sse2", 2);
  else             SF_PICK(scalar, "scalar", 1);
#elif defined(__aarch64__)
  /* Same shape as the x86 branch above, and for the same reason: a request
   * that was honoured returns with g_forced set, while an unsatisfiable one
   * warns and falls through to the automatic choice with g_forced still 0.
   * Setting it unconditionally here made sf_info()$forced report TRUE after a
   * fallback, which is the opposite of what happened. */
  if (force) {
    if (!strcmp(want, "scalar")) { SF_PICK(scalar, "scalar", 1); g_forced = 1; return; }
    if (!strcmp(want, "neon"))   { SF_PICK(neon, "neon", 2);     g_forced = 1; return; }
    unavailable = 1;
  }
  SF_PICK(neon, "neon", 2);
#else
  if (force) {
    if (!strcmp(want, "scalar")) { SF_PICK(scalar, "scalar", 1); g_forced = 1; return; }
    unavailable = 1;
  }
  SF_PICK(scalar, "scalar", 1);
#endif

  /* Warn only now that a backend is fully selected. Rf_warning can unwind the
   * C frame -- options(warn = 2), or any tryCatch/expect_warning handler
   * established by the caller -- so warning before the fallback pick left the
   * previously selected backend active and skipped the pick entirely. */
  if (unavailable)
    Rf_warning("simdfuse: instruction set '%s' is not available here; using auto",
               want);
}

int sf_isa_lanes(void) { return g_lanes; }

/* ------------------------------------------------------------- validation */

static void sf_validate(const int *code, int ncode, int nreg, int nvar,
                        int nconst) {
  for (int ip = 0; ip < ncode; ip++) {
    const int *I = code + (size_t)ip * SF_INSTR_WIDTH;
    const int op = I[0];
    if (op < 0 || op >= SF_OP_COUNT)
      Rf_error("simdfuse: opcode %d out of range at instruction %d", op, ip + 1);
    if (I[4] < 0 || I[4] >= nreg)
      Rf_error("simdfuse: destination register out of range at instruction %d", ip + 1);
    if (op == SF_OP_LOADV) {
      if (I[1] < 0 || I[1] >= nvar)
        Rf_error("simdfuse: variable index out of range at instruction %d", ip + 1);
      continue;
    }
    if (op == SF_OP_LOADC) {
      if (I[1] < 0 || I[1] >= nconst)
        Rf_error("simdfuse: constant index out of range at instruction %d", ip + 1);
      continue;
    }
    if (I[1] < 0 || I[1] >= nreg)
      Rf_error("simdfuse: operand register out of range at instruction %d", ip + 1);
    if (op == SF_OP_POWI) {
      if (I[2] < 0 || I[2] >= nconst)
        Rf_error("simdfuse: constant index out of range at instruction %d", ip + 1);
      continue;
    }
    /* remaining ops read b, and SELECT/FMA3 also read c */
    if (I[2] < 0 || I[2] >= nreg) {
      /* unary ops are emitted with b == a, so an out-of-range b is a real bug */
      Rf_error("simdfuse: operand register out of range at instruction %d", ip + 1);
    }
    if ((op == SF_OP_SELECT || op == SF_OP_FMA3 || op == SF_OP_DNORM ||
         op == SF_OP_PNORM || op == SF_OP_PLOGIS || op == SF_OP_QLOGIS) &&
        (I[3] < 0 || I[3] >= nreg))
      Rf_error("simdfuse: third operand register out of range at instruction %d", ip + 1);
  }
}

void sf_prog_init(sf_prog *P, SEXP s_code, SEXP s_consts, SEXP s_vars,
                  SEXP s_n, SEXP s_nreg, SEXP s_nprologue, SEXP s_block,
                  SEXP s_fast) {
  P->code   = INTEGER(s_code);
  P->ncode  = (int)(XLENGTH(s_code) / SF_INSTR_WIDTH);
  P->consts = REAL(s_consts);
  P->nconst = (int)XLENGTH(s_consts);
  P->nvar   = (int)XLENGTH(s_vars);
  P->nreg   = INTEGER(s_nreg)[0];
  P->fast   = Rf_asInteger(s_fast);   /* bitmask, see sf_ops.h */
  P->n      = (R_xlen_t)REAL(s_n)[0];

  if (P->ncode < 1) Rf_error("simdfuse: empty program");
  if (P->nreg < 1)  Rf_error("simdfuse: no registers requested");

  int B = INTEGER(s_block)[0];
  if (B < g_lanes) B = g_lanes;
  B = ((B + g_lanes - 1) / g_lanes) * g_lanes;
  P->B = B;

  P->nprologue = INTEGER(s_nprologue)[0];
  if (P->nprologue < 0 || P->nprologue > P->ncode)
    Rf_error("simdfuse: prologue length %d out of range", P->nprologue);

  sf_validate(P->code, P->ncode, P->nreg, P->nvar, P->nconst);
  P->final_reg = P->code[(size_t)(P->ncode - 1) * SF_INSTR_WIDTH + 4];

  const int nv = P->nvar > 0 ? P->nvar : 1;
  const void **vars = (const void **)R_alloc((size_t)nv, sizeof(void *));
  int *vscal = (int *)R_alloc((size_t)nv, sizeof(int));
  int *vtype = (int *)R_alloc((size_t)nv, sizeof(int));
  for (int i = 0; i < P->nvar; i++) {
    SEXP v = VECTOR_ELT(s_vars, i);
    switch (TYPEOF(v)) {
    case REALSXP: vtype[i] = SF_VT_DBL; vars[i] = REAL_RO(v);    break;
    case INTSXP:
    case LGLSXP:  vtype[i] = SF_VT_INT; vars[i] = INTEGER_RO(v); break;
    default:
      Rf_error("simdfuse: input %d is not a double, integer or logical vector "
               "(internal error)", i + 1);
    }
    vscal[i] = (XLENGTH(v) == 1) ? 1 : 0;
    if (!vscal[i] && XLENGTH(v) != P->n)
      Rf_error("simdfuse: input %d has length %lld but the expression is length %lld",
               i + 1, (long long)XLENGTH(v), (long long)P->n);
  }
  P->vars = vars;
  P->vscal = vscal;
  P->vtype = vtype;
}

void sf_prog_prologue(const sf_prog *P, double *regs) {
  if (P->nprologue <= 0) return;
  /* Run at full block width so every lane of the invariant registers is set,
   * whatever length the final block turns out to be. */
  g_exec(P->code, P->nprologue, P->consts, P->vars, P->vscal, P->vtype, 0,
         P->B, regs, P->nreg, P->B, P->fast, NULL, -1);
}

/* dest, when not NULL, is where the final result of this block should land.
 * The kernel writes the last instruction straight into it and the caller skips
 * its copy, but only for a block that needs no padding -- see sf_kernels.h.
 * sf_prog_wrote_dest says whether that happened. */
void sf_prog_block_to(const sf_prog *P, double *regs, R_xlen_t start, int blen,
                      double *dest) {
  const int n = P->ncode - P->nprologue;
  if (n <= 0) return;      /* wholly invariant expression: prologue did it all */
  g_exec(P->code + (size_t)P->nprologue * SF_INSTR_WIDTH, n, P->consts,
         P->vars, P->vscal, P->vtype, start, blen, regs, P->nreg, P->B,
         P->fast, dest, P->final_reg);
  /* final_reg goes through whether or not dest does: the kernel needs it to
     know which register it must not leave aliased. */
}

void sf_prog_block(const sf_prog *P, double *regs, R_xlen_t start, int blen) {
  sf_prog_block_to(P, regs, start, blen, NULL);
}

/* Whether a block of this length will have been written to dest. Mirrors the
 * kernel's own condition: only a block that is already a whole number of lanes
 * can be written in place, because anything else would round up past the end
 * of the output. The lane count depends on the dispatched kernel, so ask the
 * dispatcher rather than assuming. */
int sf_prog_dest_used(const sf_prog *P, int blen) {
  const int w = sf_isa_lanes();
  if (P->ncode - P->nprologue <= 0) return 0;
  return blen == ((blen + w - 1) / w) * w;
}

void sf_prog_reduce(const double *v, int blen, int red, int na_rm, double *out) {
  g_reduce(v, blen, red, na_rm, out);
}

void sf_prog_stencil(const double *v, R_xlen_t len, const double *w, int nw,
                     double *out) {
  g_stencil(v, len, w, nw, out);
}

double *sf_prog_final(const sf_prog *P, double *regs) {
  return regs + (size_t)P->final_reg * P->B;
}

double *sf_prog_regfile(const sf_prog *P, int nthreads) {
  const size_t per = (size_t)P->nreg * (size_t)P->B;
  double *raw = (double *)R_alloc(per * (size_t)nthreads + 8, sizeof(double));
  uintptr_t a = ((uintptr_t)raw + 63u) & ~(uintptr_t)63u;
  return (double *)a;
}

int sf_prog_threads(R_xlen_t n, R_xlen_t ntask, int requested) {
#ifdef _OPENMP
  /* Each worker must receive SF_MIN_TASKS_PER_THREAD tasks to repay the cost
   * of joining the team; see the note in sf_prog.h. An explicit request is
   * honoured only up to that cap, just as it was already capped at the task
   * count -- a thread with no task to run was never anything but overhead.
   * Capping cannot change the result, only how long it takes to get it. */
  if (n < SF_MIN_PARALLEL || ntask < 2) return 1;
  const R_xlen_t cap = ntask / SF_MIN_TASKS_PER_THREAD;
  if (cap < 2) return 1;
  int mx = omp_get_max_threads();
  if (requested < 1) requested = mx;
  if (requested > mx) requested = mx;
  if ((R_xlen_t)requested > cap) requested = (int)cap;
  return requested < 1 ? 1 : requested;
#else
  (void)n; (void)ntask; (void)requested;
  return 1;
#endif
}

/* ------------------------------------------------------------ introspection */

SEXP C_sf_op_table(void) {
  static const char *nms[] = {
#define X(nm) #nm,
    SF_OP_LIST
#undef X
    NULL };
  SEXP out = PROTECT(Rf_allocVector(INTSXP, SF_OP_COUNT));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, SF_OP_COUNT));
  for (int i = 0; i < SF_OP_COUNT; i++) {
    INTEGER(out)[i] = i;
    SET_STRING_ELT(names, i, Rf_mkChar(nms[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

SEXP C_sf_red_table(void) {
  const char *nms[] = {"none","sum","mean","max","min","prod","any","all","nobs"};
  SEXP out = PROTECT(Rf_allocVector(INTSXP, SF_RED_COUNT));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, SF_RED_COUNT));
  for (int i = 0; i < SF_RED_COUNT; i++) {
    INTEGER(out)[i] = i;
    SET_STRING_ELT(names, i, Rf_mkChar(nms[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

SEXP C_sf_info(void) {
  const char *fields[] = {"isa","lanes","forced","available","openmp",
                          "max_threads","task_elems","instr_width"};
  const int nf = 8;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, nf));
  SEXP nm  = PROTECT(Rf_allocVector(STRSXP, nf));
  for (int i = 0; i < nf; i++) SET_STRING_ELT(nm, i, Rf_mkChar(fields[i]));

  SET_VECTOR_ELT(out, 0, Rf_mkString(g_isa));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(g_lanes));
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(g_forced));

  const char *avail[8]; int na = 0;
  avail[na++] = "scalar";
#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
  if (__builtin_cpu_supports("sse2")) avail[na++] = "sse2";
  if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) avail[na++] = "avx2";
  if (__builtin_cpu_supports("avx512f")) avail[na++] = "avx512";
#elif defined(__aarch64__)
  avail[na++] = "neon";
#endif
  SEXP av = PROTECT(Rf_allocVector(STRSXP, na));
  for (int i = 0; i < na; i++) SET_STRING_ELT(av, i, Rf_mkChar(avail[i]));
  SET_VECTOR_ELT(out, 3, av);
  UNPROTECT(1);

#ifdef _OPENMP
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(TRUE));
  SET_VECTOR_ELT(out, 5, Rf_ScalarInteger(omp_get_max_threads()));
#else
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(FALSE));
  SET_VECTOR_ELT(out, 5, Rf_ScalarInteger(1));
#endif
  SET_VECTOR_ELT(out, 6, Rf_ScalarInteger(SF_TASK_ELEMS));
  SET_VECTOR_ELT(out, 7, Rf_ScalarInteger(SF_INSTR_WIDTH));
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}

SEXP C_sf_set_isa(SEXP s_isa) {
  sf_select_isa(CHAR(STRING_ELT(s_isa, 0)));
  return Rf_mkString(g_isa);
}

/* ------------------------------------------------------------ registration */

extern SEXP C_sf_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_sf_scan(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_sf_stencil(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_sf_segreduce(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"C_sf_eval",      (DL_FUNC)&C_sf_eval,      11},
  {"C_sf_scan",      (DL_FUNC)&C_sf_scan,      10},
  {"C_sf_stencil",   (DL_FUNC)&C_sf_stencil,   11},
  {"C_sf_segreduce", (DL_FUNC)&C_sf_segreduce, 12},
  {"C_sf_op_table",  (DL_FUNC)&C_sf_op_table,   0},
  {"C_sf_red_table", (DL_FUNC)&C_sf_red_table,  0},
  {"C_sf_info",      (DL_FUNC)&C_sf_info,       0},
  {"C_sf_set_isa",   (DL_FUNC)&C_sf_set_isa,    1},
  {NULL, NULL, 0}
};

attribute_visible void R_init_simdfuse(DllInfo *dll) {
  sf_select_isa("auto");
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
