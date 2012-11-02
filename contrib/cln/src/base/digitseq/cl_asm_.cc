// Includes the CPU specific cl_asm_*.cc file.

#include "cl_config.h"
#include "base/digitseq/cl_DS_endian.h"

#ifndef NO_ASM

#if defined(__m68k__) && (intCsize==16)
  #include "base/digitseq/cl_asm_m68k_.cc"
#endif

#if defined(__sparc__) && !defined(__sparc64__)
  #include "base/digitseq/cl_asm_sparc_.cc"
#endif

#if defined(__sparc64__)
  #include "base/digitseq/cl_asm_sparc64_.cc"
#endif

#if defined(__i386__)
  #include "base/digitseq/cl_asm_i386_.cc"
#endif

#if defined(__mips__) && !defined(__mipsel__)
  #include "base/digitseq/cl_asm_mips_.cc"
#endif

#if defined(__mipsel__)
  #include "base/digitseq/cl_asm_mipsel_.cc"
#endif

#if defined(__hppa__)
  #include "base/digitseq/cl_asm_hppa_.cc"
#endif

#if defined(__arm__)
  #include "base/digitseq/cl_asm_arm_.cc"
#endif

#endif // ndef NO_ASM
