// Includes the CPU specific cl_asm_*.h file.

#include "cl_config.h"
#include "base/digitseq/cl_DS_endian.h"

#ifdef _MSC_VER
/* m$vc does not grok AT&T asm */
#undef NO_ASM_LOOPS
#define NO_ASM_LOOPS
#endif

#if !defined(NO_ASM) && !defined(NO_ASM_LOOPS)

#if defined(__m68k__) && (intCsize==16)
  #include "cl_asm_m68k.h"
#endif

#if defined(__sparc__) && !defined(__sparc64__)
  #include "cl_asm_sparc.h"
#endif

#if defined(__sparc64__)
  #include "cl_asm_sparc64.h"
#endif

#if defined(__i386__)
  #include "cl_asm_i386.h"
#endif

#if (defined(__mips__) || defined(__mipsel__)) && !defined(__mips64__) && (intDsize==32)
  #include "cl_asm_mips.h"
#endif

#if defined(__hppa__) && (intDsize==32)
  #include "cl_asm_hppa.h"
#endif

#if defined(__arm__)
  #include "cl_asm_arm.h"
#endif

#endif // ndef NO_ASM
