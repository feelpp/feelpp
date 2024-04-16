#ifndef CL_CONFIG_PUBLIC_H
#define CL_CONFIG_PUBLIC_H

#include "cln/host_cpu.h"
#include "cln/version.h"

/**
 * Alignment of a `void*'. CLN needs it to distinguish between pointers
 * and immediate values.
 */
#define cl_word_alignment @ALIGNOF_VOIDP@

/* 
 * Numbers in the heap are stored as "digit" (or "limb" in GMP speak)
 * sequences. A digit is an unsigned int with sizeof(void *)*CHAR_BIT bits.
 * It should be 8 or 16 or 32 or 64 bits. If CLN is sitting on top of GMP
 * it should match mp_limb_t
 */
#cmakedefine GMP_DEMANDS_UINTD_INT
#cmakedefine GMP_DEMANDS_UINTD_LONG
#cmakedefine GMP_DEMANDS_UINTD_LONG_LONG

#endif /* CL_CONFIG_PUBLIC_H */

