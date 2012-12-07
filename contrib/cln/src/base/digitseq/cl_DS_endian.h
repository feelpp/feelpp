// Digit sequence endianness

#ifndef _CL_DS_ENDIAN_H
#define _CL_DS_ENDIAN_H

#include "base/cl_gmpconfig.h"

// Set this to 1 for big-endian digit ordering in memory,
// set this to 0 for little-endian digit ordering in memory.
// We now support both.
#ifdef CL_USE_GMP
  // Use of gmp requires CL_DS_BIG_ENDIAN_P = 0.
  #define CL_DS_BIG_ENDIAN_P  0
#else
  // In general, the digit ordering has nearly no effect on speed.
  // We have used the big-endian ordering for a long time, so let's use
  // the little-endian ordering now for at least the same time :-)
  #define CL_DS_BIG_ENDIAN_P  0
#endif

#endif /* _CL_DS_ENDIAN_H */
