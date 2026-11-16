/* generated wrapper: compiles sf_kernels.h for the avx512 instruction set */
#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
#define SF_ISA SF_ISA_AVX512
#include "sf_simd.h"
#include "sf_kernels.h"
#endif

/* ISO C has no empty translation unit. This file compiles to nothing on any
 * target that is not x86, which includes every arm64 Mac. A file-scope
 * typedef emits no symbol and no code; it is here only so that -pedantic
 * has something to see. */
typedef int sf_kernels_avx512_unused_t;
