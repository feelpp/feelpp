// System dependent definitions

#ifndef _CL_SYSDEP_H
#define _CL_SYSDEP_H

// CPU and other
#include "cl_config.h"

// char_bitsize, short_bitsize, long_bitsize, long_long_bitsize
#include "cln/intparam.h"

// The CPU's endianness
#if defined(short_little_endian) || defined(int_little_endian) || defined(long_little_endian)
  // Z80, VAX, I80X86, DECALPHA, MIPSEL, ...:
  // Low byte at low address, high byte at high address
  #if defined(CL_CPU_BIG_ENDIAN_P)
    #error "Bogus CL_CPU_BIG_ENDIAN_P!"
  #endif
  #define CL_CPU_BIG_ENDIAN_P  0
#endif
#if defined(short_big_endian) || defined(int_big_endian) || defined(long_big_endian)
  // MC680X0, SPARC, HPPA, MIPSEB, M88000, RS6000, ...:
  // High byte at low address, low byte at high address
  #if defined(CL_CPU_BIG_ENDIAN_P)
    #error "Bogus CL_CPU_BIG_ENDIAN_P!"
  #endif
  #define CL_CPU_BIG_ENDIAN_P  1
#endif
#if !defined(CL_CPU_BIG_ENDIAN_P)
  #error "Bogus CL_CPU_BIG_ENDIAN_P!"
#endif

// Auswahl der Floating-Point-Fähigkeiten:
// FAST_DOUBLE sollte definiert werden, wenn ein Floating-Point-Coprozessor
// vorhanden ist, dessen `double'-Typ IEEE-Floating-Points mit 64 Bits sind.
// FAST_FLOAT sollte definiert werden, wenn ein Floating-Point-Coprozessor
// vorhanden ist, dessen `float'-Typ IEEE-Floating-Points mit 32 Bits sind,
// und der C++-Compiler auch `float'- und nicht `double'-Operationen generiert.
#if defined(__sparc__) || defined(__sparc64__) || defined(__hppa__) || defined(__m88k__) || defined(__rs6000__)
  #define FAST_DOUBLE
  #define FAST_FLOAT
#endif
#if defined(__i386__) && (defined(linux) || defined(__linux__) || defined(NeXT))
  // Linux hat einen funktionierenden Floating-Point-Coprozessor-Emulator.
  // NeXTstep läuft sowieso nur mit Floating-Point-Coprozessor.
  // Aber auf Intel-Pentium-Prozessoren ist die FPU fehlerhaft.
  #define FAST_DOUBLE
  #define FAST_FLOAT
#endif
#if defined(__arm__)
  // Bei Integers ist der Prozessor Little-Endian, bei Double-Floats Big-Endian!
  #undef FAST_DOUBLE
#endif

// Macros for internal use.
#include "base/cl_macros.h"

// Elementary types.
#include "cln/types.h"

// Dependencies among modules.
#include "cln/modules.h"

#endif /* _CL_SYSDEP_H */
