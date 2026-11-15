/* generated wrapper: compiles sf_kernels.h for the sse2 instruction set */
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
#define SF_ISA SF_ISA_SSE2
#include "sf_simd.h"
#include "sf_kernels.h"
#endif

/* ISO C has no empty translation unit. This file compiles to nothing on any
 * target that is not x86, which includes every arm64 Mac. A file-scope
 * typedef emits no symbol and no code; it is here only so that -pedantic
 * has something to see. */
typedef int sf_kernels_sse2_unused_t;
