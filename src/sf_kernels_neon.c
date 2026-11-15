/* generated wrapper: compiles sf_kernels.h for the neon instruction set */
#if defined(__aarch64__)
#define SF_ISA SF_ISA_NEON
#include "sf_simd.h"
#include "sf_kernels.h"
#endif

/* ISO C has no empty translation unit. This file compiles to nothing on any
 * target that is not aarch64, which includes every x86 build. A file-scope
 * typedef emits no symbol and no code; it is here only so that -pedantic
 * has something to see. */
typedef int sf_kernels_neon_unused_t;
