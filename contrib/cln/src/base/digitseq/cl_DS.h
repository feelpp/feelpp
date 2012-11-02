// Digit sequence arithmetic

#ifndef _CL_DS_H
#define _CL_DS_H

#include "cln/types.h"
#include "base/cl_gmpconfig.h"
#include "base/digit/cl_D.h"
#include "base/digitseq/cl_DS_endian.h"
#include "base/cl_alloca.h"

namespace cln {

// Digit Sequence (DS)
// a memory range with n digits (n an uintC),
// between two pointers MSDptr and LSDptr.
#if CL_DS_BIG_ENDIAN_P
//  MSDptr                  LSDptr = MSDptr+n
// | MSD ............. LSD |
// [short: MSDptr/n/LSDptr ]
// In C: uintD* MSDptr, uintC len, MSDptr[0] ... MSDptr[len-1] are the digits.
#else
//  LSDptr                  MSDptr = LSDptr+n
// | LSD ............. MSD |
// In C: uintD* LSDptr, uintC len, LSDptr[0] ... LSDptr[len-1] are the digits.
#endif
// If n = 0, this represents the number 0.
// If n > 0, the most significant bit (i.e. bit (intDsize-1) of
//           MSDptr[CL_DS_BIG_ENDIAN_P?0:-1]) is the sign bit. If the sign
//           bit were repeated infinitely often, one would get an
//           "infinite bit sequence".
//
// A Normalised Digit Sequence (NDS) is one for which the MSD is necessary,
// i.e. n = 0 or (n > 0 and the most significant intDsize+1 bits are not
// all the same).

// Unsigned Digit Sequence (UDS)
// like DS, but without sign.
//
// Normalized Unsigned Digit Sequence (NUDS):
// an UDS for which the MSD is necessary, i.e. n = 0 or
// (n > 0 and the most significant intDsize bits are not all zero).

// For the construction of constant DS, using "digit_header":
#define D1(byte0)  (uintD)(byte0)
#define D2(byte0,byte1)  (((uintD)(byte0)<<8)|(uintD)(byte1))
#define D4(byte0,byte1,byte2,byte3)  (((uintD)(byte0)<<24)|((uintD)(byte1)<<16)|((uintD)(byte2)<<8)|((uintD)(byte3)))
#define D8(byte0,byte1,byte2,byte3,byte4,byte5,byte6,byte7)  (((uintD)(byte0)<<56)|((uintD)(byte1)<<48)|((uintD)(byte2)<<40)|((uintD)(byte3)<<32)|((uintD)(byte4)<<24)|((uintD)(byte5)<<16)|((uintD)(byte6)<<8)|((uintD)(byte7)))

// i386-pc-solaris #defines DS.
// Grr...
#undef DS

struct DS {
	uintD* MSDptr;
	uintC len;
	uintD* LSDptr;
};

// Endianness independent access of digit sequences:
// mspref(MSDptr,i)    access a most significant digit
// lspref(LSDptr,i)    access a least significant digit
// msshrink(MSDptr)    shrinks the DS by throwing away the MSD
// msprefnext(MSDptr)  combines mspref(MSDptr,0) and msshrink(MSDptr)
// lsshrink(LSDptr)    shrinks the DS by throwing away the LSD
// lsprefnext(LSDptr)  combines lspref(LSDptr,0) and lsshrink(LSDptr)
// mspop     pointer operator corresponding to msshrink, arg is widened to an uintP
// lspop     pointer operator corresponding to lsshrink, arg is widened to an uintP
#if CL_DS_BIG_ENDIAN_P
  #define mspref(p,i)  (p)[i]
  #define lspref(p,i)  (p)[-(uintP)(i)-1]
  #define msshrink(p)  (p)++
  #define msprefnext(p)  (*(p)++)
  #define lsshrink(p)  (p)--
  #define lsprefnext(p)  (*--(p))
  #define mspop  +
  #define lspop  -
#else
  #define mspref(p,i)  (p)[-(uintP)(i)-1]
  #define lspref(p,i)  (p)[i]
  #define msshrink(p)  (p)--
  #define msprefnext(p)  (*--(p))
  #define lsshrink(p)  (p)++
  #define lsprefnext(p)  (*(p)++)
  #define mspop  -
  #define lspop  +
#endif

// Endianness independent macros for turning an array into a digit sequence.
// arrayMSDptr(array,length)  returns the MSDptr of array[0..length-1]
// arrayLSDptr(array,length)  returns the LSDptr of array[0..length-1]
#if CL_DS_BIG_ENDIAN_P
  #define arrayMSDptr(array,length)  &(array)[0]
  #define arrayLSDptr(array,length)  &(array)[length]
#else
  #define arrayMSDptr(array,length)  &(array)[length]
  #define arrayLSDptr(array,length)  &(array)[0]
#endif
#define arrayLSref(array,length,i)  lspref(arrayLSDptr(array,length),i)


// These functions on digit sequences are either inline C++ functions
// or external assembler functions (see files cl_asm_*).


// See which functions are defined as external functions.
#include "base/digitseq/cl_asm.h"


// Declare the external functions.

extern "C" {

#ifdef COPY_LOOPS

extern uintD* copy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count);

extern uintD* copy_loop_down (const uintD* sourceptr, uintD* destptr, uintC count);

#endif

#ifdef FILL_LOOPS

extern uintD* fill_loop_up (uintD* destptr, uintC count, uintD filler);

extern uintD* fill_loop_down (uintD* destptr, uintC count, uintD filler);

#endif

#ifdef CLEAR_LOOPS

extern uintD* clear_loop_up (uintD* destptr, uintC count);

extern uintD* clear_loop_down (uintD* destptr, uintC count);

#endif

#ifdef TEST_LOOPS

extern bool test_loop_up (const uintD* ptr, uintC count);

extern bool test_loop_down (const uintD* ptr, uintC count);

#endif

#if CL_DS_BIG_ENDIAN_P

#ifdef LOG_LOOPS

extern void or_loop_up (uintD* xptr, const uintD* yptr, uintC count);

extern void xor_loop_up (uintD* xptr, const uintD* yptr, uintC count);

extern void and_loop_up (uintD* xptr, const uintD* yptr, uintC count);

extern void eqv_loop_up (uintD* xptr, const uintD* yptr, uintC count);

extern void nand_loop_up (uintD* xptr, const uintD* yptr, uintC count);

extern void nor_loop_up (uintD* xptr, const uintD* yptr, uintC count);

extern void andc2_loop_up (uintD* xptr, const uintD* yptr, uintC count);

extern void orc2_loop_up (uintD* xptr, const uintD* yptr, uintC count);

extern void not_loop_up (uintD* xptr, uintC count);

#endif

#ifdef TEST_LOOPS

extern bool and_test_loop_up (const uintD* xptr, const uintD* yptr, uintC count);

extern cl_signean compare_loop_up (const uintD* xptr, const uintD* yptr, uintC count);

#endif

#ifdef ADDSUB_LOOPS

extern uintD add_loop_down (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count);

extern uintD addto_loop_down (const uintD* sourceptr, uintD* destptr, uintC count);

extern uintD inc_loop_down (uintD* ptr, uintC count);

extern uintD sub_loop_down (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count);

extern uintD subx_loop_down (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);

extern uintD subfrom_loop_down (const uintD* sourceptr, uintD* destptr, uintC count);

extern uintD dec_loop_down (uintD* ptr, uintC count);

extern uintD neg_loop_down (uintD* ptr, uintC count);

#endif

#ifdef SHIFT_LOOPS

extern uintD shift1left_loop_down (uintD* ptr, uintC count);

extern uintD shiftleft_loop_down (uintD* ptr, uintC count, uintC i, uintD carry);

extern uintD shiftleftcopy_loop_down (const uintD* sourceptr, uintD* destptr, uintC count, uintC i);

extern uintD shift1right_loop_up (uintD* ptr, uintC count, uintD carry);

extern uintD shiftright_loop_up (uintD* ptr, uintC count, uintC i);

extern uintD shiftrightsigned_loop_up (uintD* ptr, uintC count, uintC i);

extern uintD shiftrightcopy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry);

#endif

#ifdef MUL_LOOPS

extern uintD mulusmall_loop_down (uintD digit, uintD* ptr, uintC len, uintD newdigit);

extern void mulu_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len);

extern uintD muluadd_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len);

extern uintD mulusub_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len);

#endif

#ifdef DIV_LOOPS

extern uintD divu_loop_up (uintD digit, uintD* ptr, uintC len);

extern uintD divucopy_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len);

#endif

#else // !CL_DS_BIG_ENDIAN_P

#ifdef LOG_LOOPS

extern void or_loop_down (uintD* xptr, const uintD* yptr, uintC count);

extern void xor_loop_down (uintD* xptr, const uintD* yptr, uintC count);

extern void and_loop_down (uintD* xptr, const uintD* yptr, uintC count);

extern void eqv_loop_down (uintD* xptr, const uintD* yptr, uintC count);

extern void nand_loop_down (uintD* xptr, const uintD* yptr, uintC count);

extern void nor_loop_down (uintD* xptr, const uintD* yptr, uintC count);

extern void andc2_loop_down (uintD* xptr, const uintD* yptr, uintC count);

extern void orc2_loop_down (uintD* xptr, const uintD* yptr, uintC count);

extern void not_loop_down (uintD* xptr, uintC count);

#endif

#ifdef TEST_LOOPS

extern bool and_test_loop_down (const uintD* xptr, const uintD* yptr, uintC count);

extern cl_signean compare_loop_down (const uintD* xptr, const uintD* yptr, uintC count);

#endif

#ifdef ADDSUB_LOOPS

extern uintD add_loop_up (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count);

extern uintD addto_loop_up (const uintD* sourceptr, uintD* destptr, uintC count);

extern uintD inc_loop_up (uintD* ptr, uintC count);

extern uintD sub_loop_up (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count);

extern uintD subx_loop_up (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count, uintD carry);

extern uintD subfrom_loop_up (const uintD* sourceptr, uintD* destptr, uintC count);

extern uintD dec_loop_up (uintD* ptr, uintC count);

extern uintD neg_loop_up (uintD* ptr, uintC count);

#endif

#ifdef SHIFT_LOOPS

extern uintD shift1left_loop_up (uintD* ptr, uintC count);

extern uintD shiftleft_loop_up (uintD* ptr, uintC count, uintC i, uintD carry);

extern uintD shiftleftcopy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count, uintC i);

extern uintD shift1right_loop_down (uintD* ptr, uintC count, uintD carry);

extern uintD shiftright_loop_down (uintD* ptr, uintC count, uintC i);

extern uintD shiftrightsigned_loop_down (uintD* ptr, uintC count, uintC i);

extern uintD shiftrightcopy_loop_down (const uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry);

#endif

#ifdef MUL_LOOPS

extern uintD mulusmall_loop_up (uintD digit, uintD* ptr, uintC len, uintD newdigit);

extern void mulu_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len);

extern uintD muluadd_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len);

extern uintD mulusub_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len);

#endif

#ifdef DIV_LOOPS

extern uintD divu_loop_down (uintD digit, uintD* ptr, uintC len);

extern uintD divucopy_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len);

#endif

#endif // !CL_DS_BIG_ENDIAN_P

// Independently of CL_DS_BIG_ENDIAN_P:

#ifdef TEST_LOOPS

extern cl_signean compare_loop_up (const uintD* xptr, const uintD* yptr, uintC count);

#endif

#ifdef LOG_LOOPS

extern void xor_loop_up (uintD* xptr, const uintD* yptr, uintC count);

#endif

#ifdef SHIFT_LOOPS

extern uintD shiftleftcopy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count, uintC i);

extern void shiftxor_loop_up (uintD* xptr, const uintD* yptr, uintC count, uintC i);

#endif

} // "C" extern.


#if defined(CL_USE_GMP)

// Supersede the functions by wrappers around calls to gmp mpn,
// for those functions where gmp is believed to be faster.

}  // namespace cln

#include <gmp.h>
// Argh, gmp.h includes <stddef.h> which erases the definition of offsetof
// that we have provided in cl_offsetof.h. Restore it.
#include "base/cl_offsetof.h"

namespace cln {

#if 0 // not worth it, since gmp's mpn_cmp is not optimized
inline cl_signean compare_loop_down (const uintD* xptr, const uintD* yptr, uintC count)
{
	return mpn_cmp(xptr-count,yptr-count,count);
}
#endif

inline uintD add_loop_up (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count)
{
	if (count == 0)
		return 0;
	return mpn_add_n(destptr,sourceptr1,sourceptr2,count);
}

inline uintD addto_loop_up (const uintD* sourceptr, uintD* destptr, uintC count)
{
	if (count == 0)
		return 0;
	return mpn_add_n(destptr,destptr,sourceptr,count);
}

inline uintD inc_loop_up (uintD* ptr, uintC count)
{
	if (count == 0)
		return 1;
	return mpn_add_1(ptr,ptr,count,1);
}

inline uintD sub_loop_up (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count)
{
	if (count == 0)
		return 0;
	return mpn_sub_n(destptr,sourceptr1,sourceptr2,count);
}

inline uintD subx_loop_up (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count, uintD carry)
{
	if (count == 0)
		return carry;
	var uintD res_carry = mpn_sub_n(destptr,sourceptr1,sourceptr2,count);
	if (carry)
		res_carry |= mpn_sub_1(destptr,destptr,count,1);
	return res_carry;
}

inline uintD subfrom_loop_up (const uintD* sourceptr, uintD* destptr, uintC count)
{
	if (count == 0)
		return 0;
	return mpn_sub_n(destptr,destptr,sourceptr,count);
}

inline uintD dec_loop_up (uintD* ptr, uintC count)
{
	if (count == 0)
		return (uintD)(-1);
	return -mpn_sub_1(ptr,ptr,count,1);
}

#if !defined(ADDSUB_LOOPS)
// No equivalent for this in gmp. But we need this function, so write it in C.
inline uintD neg_loop_up (uintD* ptr, uintC count)
{
	// erstes Digit /=0 suchen:
	until (count==0) { if (!(*ptr == 0)) goto L1; ptr++; count--; }
	return 0;
    L1: // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
	*ptr = - *ptr; count--; // 1 Digit negieren
	dotimesC(count,count, { ptr++; *ptr = ~ *ptr; } ); // alle anderen Digits invertieren
	return (uintD)(-1);
}
#endif

#define ADDSUB_LOOPS

inline uintD shift1left_loop_up (uintD* ptr, uintC count)
{
	if (count == 0)
		return 0;
	return mpn_lshift(ptr,ptr,count,1);
}

inline uintD shiftleft_loop_up (uintD* ptr, uintC count, uintC i, uintD carry)
{
	if (count == 0)
		return carry;
	var uintD res_carry = mpn_lshift(ptr,ptr,count,i);
	ptr[0] |= carry;
	return res_carry;
}

inline uintD shiftleftcopy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count, uintC i)
{
	if (count == 0)
		return 0;
	return mpn_lshift(destptr,sourceptr,count,i);
}

inline uintD shift1right_loop_down (uintD* ptr, uintC count, uintD carry)
{
	if (count == 0)
		return carry;
	var uintD res_carry = mpn_rshift(ptr-count,ptr-count,count,1);
	if (carry)
		ptr[-1] |= bit(intDsize-1);
	return res_carry;
}

inline uintD shiftright_loop_down (uintD* ptr, uintC count, uintC i)
{
	if (count == 0)
		return 0;
	return mpn_rshift(ptr-count,ptr-count,count,i);
}

inline uintD shiftrightsigned_loop_down (uintD* ptr, uintC count, uintC i)
{
	var uintD carry = ((sintD)ptr[-1] >> (intDsize-1)) << (intDsize-i);
	var uintD res_carry = mpn_rshift(ptr-count,ptr-count,count,i);
	ptr[-1] |= carry;
	return res_carry;
}

inline uintD shiftrightcopy_loop_down (const uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry)
{
	carry = carry << (intDsize-i);
	if (count == 0)
		return carry;
	var uintD res_carry = mpn_rshift(destptr-count,sourceptr-count,count,i);
	destptr[-1] |= carry;
	return res_carry;
}

#define SHIFT_LOOPS

inline uintD mulusmall_loop_up (uintD digit, uintD* ptr, uintC len, uintD newdigit)
{
	if (len == 0)
		return newdigit;
	var uintD res_carry = mpn_mul_1(ptr,ptr,len,digit);
	res_carry += mpn_add_1(ptr,ptr,len,newdigit);
	return res_carry;
}

inline void mulu_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
{
	destptr[len] = (len==0 ? 0 : mpn_mul_1(destptr,sourceptr,len,digit));
}

inline uintD muluadd_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
{
	if (len == 0)
		return 0;
	return mpn_addmul_1(destptr,sourceptr,len,digit);
}

inline uintD mulusub_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
{
	if (len == 0)
		return 0;
	return mpn_submul_1(destptr,sourceptr,len,digit);
}

#define MUL_LOOPS

inline uintD divu_loop_up (uintD digit, uintD* ptr, uintC len)
{
	return mpn_divrem_1(ptr,0,ptr,len,digit);
}

inline uintD divu_loop_down (uintD digit, uintD* ptr, uintC len)
{
	return mpn_divrem_1(ptr-len,0,ptr-len,len,digit);
}

inline uintD divucopy_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
{
	return mpn_divrem_1(destptr,0,sourceptr,len,digit);
}

inline uintD divucopy_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
{
	return mpn_divrem_1(destptr-len,0,sourceptr-len,len,digit);
}

#define DIV_LOOPS

#endif // defined(CL_USE_GMP)


// Define the missing functions as inline functions.

#ifndef COPY_LOOPS

// Kopierschleife:
// destptr = copy_loop_up(sourceptr,destptr,count);
// kopiert count (uintC>=0) Digits aufwärts von sourceptr nach destptr
// und liefert das neue destptr.
  inline uintD* copy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count)
    { dotimesC(count,count, { *destptr++ = *sourceptr++; } );
      return destptr;
    }

// Kopierschleife:
// destptr = copy_loop_down(sourceptr,destptr,count);
// kopiert count (uintC>=0) Digits abwärts von sourceptr nach destptr
// und liefert das neue destptr.
  inline uintD* copy_loop_down (const uintD* sourceptr, uintD* destptr, uintC count)
    { dotimesC(count,count, { *--destptr = *--sourceptr; } );
      return destptr;
    }

#endif

#ifndef FILL_LOOPS

// Füllschleife:
// destptr = fill_loop_up(destptr,count,filler);
// kopiert count (uintC>=0) mal das Digit filler aufwärts nach destptr
// und liefert das neue destptr.
  inline uintD* fill_loop_up (uintD* destptr, uintC count, uintD filler)
    { dotimesC(count,count, { *destptr++ = filler; } );
      return destptr;
    }

// Füllschleife:
// destptr = fill_loop_down(destptr,count,filler);
// kopiert count (uintC>=0) mal das Digit filler abwärts nach destptr
// und liefert das neue destptr.
  inline uintD* fill_loop_down (uintD* destptr, uintC count, uintD filler)
    { dotimesC(count,count, { *--destptr = filler; } );
      return destptr;
    }

#endif

#ifndef CLEAR_LOOPS

// Lösch-Schleife:
// destptr = clear_loop_up(destptr,count);
// löscht count (uintC>=0) Digits aufwärts ab destptr
// und liefert das neue destptr.
  inline uintD* clear_loop_up (uintD* destptr, uintC count)
    { dotimesC(count,count, { *destptr++ = 0; } );
      return destptr;
    }

// Lösch-Schleife:
// destptr = clear_loop_down(destptr,count);
// löscht count (uintC>=0) Digits abwärts ab destptr
// und liefert das neue destptr.
  inline uintD* clear_loop_down (uintD* destptr, uintC count)
    { dotimesC(count,count, { *--destptr = 0; } );
      return destptr;
    }

#endif

#ifndef TEST_LOOPS

// Test-Schleife:
// test_loop_up(ptr,count)
// testet count (uintC>=0) Digits aufwärts ab ptr, ob darunter eines /=0 ist.
// Ergebnis /=0, falls ja.
  inline bool test_loop_up (const uintD* ptr, uintC count)
    { dotimesC(count,count, { if (*ptr++) return true; } );
      return false;
    }

// Test-Schleife:
// test_loop_down(ptr,count)
// testet count (uintC>=0) Digits abwärts ab ptr, ob darunter eines /=0 ist.
// Ergebnis /=0, falls ja.
  inline bool test_loop_down (const uintD* ptr, uintC count)
    { dotimesC(count,count, { if (*--ptr) return true; } );
      return false;
    }

#endif

#if CL_DS_BIG_ENDIAN_P

#ifndef LOG_LOOPS

// OR-Schleife:
// or_loop_up(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits aufwärts ab xptr und ab yptr
// mit Ziel ab xptr durch OR.
  inline void or_loop_up (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *xptr++ |= *yptr++; } ); }

// XOR-Schleife:
// xor_loop_up(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits aufwärts ab xptr und ab yptr
// mit Ziel ab xptr durch XOR.
  inline void xor_loop_up (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *xptr++ ^= *yptr++; } ); }

// AND-Schleife:
// and_loop_up(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits aufwärts ab xptr und ab yptr
// mit Ziel ab xptr durch AND.
  inline void and_loop_up (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *xptr++ &= *yptr++; } ); }

// EQV-Schleife:
// eqv_loop_up(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits aufwärts ab xptr und ab yptr
// mit Ziel ab xptr durch EQV (NOT XOR).
  inline void eqv_loop_up (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count,
        {var uintD temp = ~ (*xptr ^ *yptr++); *xptr++ = temp; }
        );
    }

// NAND-Schleife:
// nand_loop_up(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits aufwärts ab xptr und ab yptr
// mit Ziel ab xptr durch NAND (NOT AND).
  inline void nand_loop_up (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count,
        {var uintD temp = ~ (*xptr & *yptr++); *xptr++ = temp; }
        );
    }

// NOR-Schleife:
// nor_loop_up(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits aufwärts ab xptr und ab yptr
// mit Ziel ab xptr durch NOR (NOT OR).
  inline void nor_loop_up (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count,
        {var uintD temp = ~ (*xptr | *yptr++); *xptr++ = temp; }
        );
    }

// ANDC2-Schleife:
// andc2_loop_up(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits aufwärts ab xptr und ab yptr
// mit Ziel ab xptr durch ANDC2 (AND NOT).
  inline void andc2_loop_up (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *xptr++ &= ~(*yptr++); } ); }

// ORC2-Schleife:
// orc2_loop_up(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits aufwärts ab xptr und ab yptr
// mit Ziel ab xptr durch ORC2 (OR NOT).
  inline void orc2_loop_up (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *xptr++ |= ~(*yptr++); } ); }

// NOT-Schleife:
// not_loop_up(xptr,count);
// verknüpft count (uintC>0) Digits aufwärts ab xptr mit Ziel ab xptr
// durch NOT.
  inline void not_loop_up (uintD* xptr, uintC count)
    { dotimespC(count,count,
        {var uintD temp = ~ (*xptr); *xptr++ = temp; }
        );
    }

#endif

#ifndef TEST_LOOPS

// AND-Test-Schleife:
// and_test_loop_up(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits aufwärts ab xptr und ab yptr durch AND
// und testet, ob sich dabei ein Digit /=0 ergibt. Ergebnis true, falls ja.
  inline bool and_test_loop_up (const uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { if (*xptr++ & *yptr++) return true; } );
      return false;
    }

// Vergleichsschleife:
// result = compare_loop_up(xptr,yptr,count);
// vergleicht nacheinander xptr[0] mit yptr[0], xptr[1] mit yptr[1], usw.,
// insgesamt count Digits, und liefert 0 falls alle gleich sind,
// +1 falls zuerst ein xptr[i]>yptr[i] ist,
// -1 falls zuerst ein xptr[i]<yptr[i] ist.
  inline cl_signean compare_loop_up (const uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count,
        { if (!(*xptr++ == *yptr++))
            // verschiedene Digits gefunden
            return (*--xptr > *--yptr ? signean_plus : signean_minus);
        });
      return signean_null; // alle Digits gleich
    }

#endif

#ifndef ADDSUB_LOOPS

// Additionsschleife:
// übertrag = add_loop_down(sourceptr1,sourceptr2,destptr,count);
// addiert count (uintC>=0) Digits abwärts von sourceptr1, von sourceptr2
// abwärts nach destptr und liefert den Übertrag (0 oder /=0, was 1 bedeutet).
  inline uintD add_loop_down (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count)
    { var uintD source1;
      var uintD source2;
      if (!(count==0))
      do { source1 = *--sourceptr1;
           source2 = *--sourceptr2;
           *--destptr = source1 + source2;
           if (source1 > (uintD)(~source2)) goto carry_1;
           carry_0:
           count--;
         }
         until (count==0);
      return 0;
      do { source1 = *--sourceptr1;
           source2 = *--sourceptr2;
           *--destptr = source1 + source2 + 1;
           if (source1 < (uintD)(~source2)) goto carry_0;
           carry_1:
           count--;
         }
         until (count==0);
      return 1;
    }

// Additionsschleife:
// übertrag = addto_loop_down(sourceptr,destptr,count);
// addiert count (uintC>=0) Digits abwärts von sourceptr, von destptr
// abwärts nach destptr und liefert den Übertrag (0 oder /=0, was 1 bedeutet).
  inline uintD addto_loop_down (const uintD* sourceptr, uintD* destptr, uintC count)
    { var uintD source1;
      var uintD source2;
      if (!(count==0))
      do { source1 = *--sourceptr;
           source2 = *--destptr;
           *destptr = source1 + source2;
           if (source1 > (uintD)(~source2)) goto carry_1;
           carry_0:
           count--;
         }
         until (count==0);
      return 0;
      do { source1 = *--sourceptr;
           source2 = *--destptr;
           *destptr = source1 + source2 + 1;
           if (source1 < (uintD)(~source2)) goto carry_0;
           carry_1:
           count--;
         }
         until (count==0);
      return 1;
    }

// Incrementierschleife:
// übertrag = inc_loop_down(ptr,count);
// incrementiert count (uintC>=0) Digits abwärts von ptr, so lange bis kein
// Übertrag mehr auftritt und liefert den Übertrag (0 oder /=0, was 1 bedeutet).
  inline uintD inc_loop_down (uintD* ptr, uintC count)
    { dotimesC(count,count,
        { if (!( ++(*--ptr) == 0 )) return 0; } // kein weiterer Übertrag
        );
      return 1; // weiterer Übertrag
    }

// Subtraktionsschleife:
// übertrag = sub_loop_down(sourceptr1,sourceptr2,destptr,count);
// subtrahiert count (uintC>=0) Digits abwärts von sourceptr1, von sourceptr2
// abwärts nach destptr und liefert den Übertrag (0 oder /=0, was -1 bedeutet).
  inline uintD sub_loop_down (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count)
    { var uintD source1;
      var uintD source2;
      if (!(count==0))
      do { source1 = *--sourceptr1;
           source2 = *--sourceptr2;
           *--destptr = source1 - source2;
           if (source1 < source2) goto carry_1;
           carry_0:
           count--;
         }
         until (count==0);
      return 0;
      do { source1 = *--sourceptr1;
           source2 = *--sourceptr2;
           *--destptr = source1 - source2 - 1;
           if (source1 > source2) goto carry_0;
           carry_1:
           count--;
         }
         until (count==0);
      return (uintD)(-1);
    }

// Subtraktionsschleife:
// übertrag = subx_loop_down(sourceptr1,sourceptr2,destptr,count,carry);
// subtrahiert count (uintC>=0) Digits abwärts von sourceptr1 und addiert
// einen Carry (0 oder -1), von sourceptr2 abwärts nach destptr und
// liefert den Übertrag (0 oder /=0, was -1 bedeutet).
  inline uintD subx_loop_down (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count, uintD carry)
    { var uintD source1;
      var uintD source2;
      if (carry==0)
        { if (!(count==0))
            do { source1 = *--sourceptr1;
                 source2 = *--sourceptr2;
                 *--destptr = source1 - source2;
                 if (source1 < source2) goto carry_1;
                 carry_0:
                 count--;
               }
               until (count==0);
          return 0;
        }
        else
        { if (!(count==0))
            do { source1 = *--sourceptr1;
                 source2 = *--sourceptr2;
                 *--destptr = source1 - source2 - 1;
                 if (source1 > source2) goto carry_0;
                 carry_1:
                 count--;
               }
               until (count==0);
          return (uintD)(-1);
    }   }

// Subtraktionsschleife:
// übertrag = subfrom_loop_down(sourceptr,destptr,count);
// subtrahiert count (uintC>=0) Digits abwärts von sourceptr, von destptr
// abwärts nach destptr (dest := dest - source)
// und liefert den Übertrag (0 oder /=0, was -1 bedeutet).
  inline uintD subfrom_loop_down (const uintD* sourceptr, uintD* destptr, uintC count)
    { var uintD source1;
      var uintD source2;
      if (!(count==0))
      do { source1 = *--destptr;
           source2 = *--sourceptr;
           *destptr = source1 - source2;
           if (source1 < source2) goto carry_1;
           carry_0:
           count--;
         }
         until (count==0);
      return 0;
      do { source1 = *--destptr;
           source2 = *--sourceptr;
           *destptr = source1 - source2 - 1;
           if (source1 > source2) goto carry_0;
           carry_1:
           count--;
         }
         until (count==0);
      return (uintD)(-1);
    }

// Decrementierschleife:
// übertrag = dec_loop_down(ptr,count);
// decrementiert count (uintC>=0) Digits abwärts von ptr, so lange bis kein
// Übertrag mehr auftritt und liefert den Übertrag (0 oder -1).
  inline uintD dec_loop_down (uintD* ptr, uintC count)
    { dotimesC(count,count,
        { if (!( (*--ptr)-- == 0 )) return 0; } // kein weiterer Übertrag
        );
      return (uintD)(-1); // weiterer Übertrag
    }

// Negierschleife:
// übertrag = neg_loop_down(ptr,count);
// negiert count (uintC>=0) Digits abwärts von ptr,
// und liefert den Übertrag (0 oder -1).
  inline uintD neg_loop_down (uintD* ptr, uintC count)
    { // erstes Digit /=0 suchen:
      until (count==0) { if (!(*--ptr == 0)) goto L1; count--; }
      return 0;
      L1: // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
      *ptr = - *ptr; count--; // 1 Digit negieren
      dotimesC(count,count, { --ptr; *ptr = ~ *ptr; } ); // alle anderen Digits invertieren
      return (uintD)(-1);
    }

#endif

#ifndef SHIFT_LOOPS

// Schiebeschleife um 1 Bit nach links:
// übertrag = shift1left_loop_down(ptr,count);
// schiebt count (uintC>=0) Digits abwärts von ptr um 1 Bit nach links,
// und liefert den Übertrag (0 oder /=0, was 1 bedeutet).
  #if HAVE_DD
  inline uintD shift1left_loop_down (uintD* ptr, uintC count)
    { var uintDD accu = 0;
      dotimesC(count,count,
        { accu = ((uintDD)(*--ptr)<<1)+accu; *ptr = lowD(accu);
          accu = (uintDD)(highD(accu));
        });
      return (uintD)accu;
    }
  #else
  inline uintD shift1left_loop_down (uintD* ptr, uintC count)
    { var uintD carry = 0;
      dotimesC(count,count,
        { var uintD accu = *--ptr;
          *ptr = (accu<<1) | carry;
          carry = accu>>(intDsize-1);
        });
      return carry;
    }
  #endif

// Schiebeschleife um i Bits nach links:
// übertrag = shiftleft_loop_down(ptr,count,i,übertrag_init);
// schiebt count (uintC>=0) Digits abwärts von ptr um i Bits (0<i<intDsize)
// nach links, schiebt dabei die i Bits aus übertrag_init rechts rein,
// und liefert den Übertrag (was links rauskommt, >=0, <2^i).
  #if HAVE_DD
  inline uintD shiftleft_loop_down (uintD* ptr, uintC count, uintC i, uintD carry)
    { var uintDD accu = (uintDD)carry;
      dotimesC(count,count,
        { accu = ((uintDD)(*--ptr)<<i)+accu; *ptr = lowD(accu);
          accu = (uintDD)(highD(accu));
        });
      return (uintD)accu;
    }
  #else
  inline uintD shiftleft_loop_down (uintD* ptr, uintC count, uintC i, uintD carry)
    { var uintC j = intDsize-i;
      dotimesC(count,count,
        { var uintD accu = *--ptr;
          *ptr = (accu<<i) | carry;
          carry = accu>>j;
        });
      return carry;
    }
  #endif

// Schiebe- und Kopierschleife um i Bits nach links:
// übertrag = shiftleftcopy_loop_down(sourceptr,destptr,count,i);
// kopiert count (uintC>=0) Digits abwärts von sourceptr nach destptr
// und schiebt sie dabei um i Bits (0<i<intDsize) nach links,
// wobei ganz rechts mit i Nullbits aufgefüllt wird,
// und liefert den Übertrag (was links rauskommt, >=0, <2^i).
  #if HAVE_DD
  inline uintD shiftleftcopy_loop_down (const uintD* sourceptr, uintD* destptr, uintC count, uintC i)
    { var uintDD accu = 0;
      dotimesC(count,count,
        { accu = ((uintDD)(*--sourceptr)<<i)+accu; *--destptr = lowD(accu);
          accu = (uintDD)(highD(accu));
        });
      return (uintD)accu;
    }
  #else
  inline uintD shiftleftcopy_loop_down (const uintD* sourceptr, uintD* destptr, uintC count, uintC i)
    { var uintC j = intDsize-i;
      var uintD carry = 0;
      dotimesC(count,count,
        { var uintD accu = *--sourceptr;
          *--destptr = (accu<<i) | carry;
          carry = accu>>j;
        });
      return carry;
    }
  #endif

// Schiebeschleife um 1 Bit nach rechts:
// übertrag = shift1right_loop_up(ptr,count,übertrag_init);
// schiebt count (uintC>=0) Digits aufwärts von ptr um 1 Bit nach rechts,
// wobei links das Bit übertrag_init (sollte =0 oder =-1 sein) hineingeschoben
// wird, und liefert den Übertrag (0 oder /=0, was 1 bedeutet).
  #if HAVE_DD
  inline uintD shift1right_loop_up (uintD* ptr, uintC count, uintD carry)
    { var uintDD accu = (sintDD)(sintD)carry & ((uintDD)1 << (2*intDsize-1)); // 0 oder bit(2*intDsize-1)
      dotimesC(count,count,
        { accu = (highlowDD_0(*ptr)>>1)+accu; *ptr++ = highD(accu);
          accu = highlowDD_0(lowD(accu));
        });
      return highD(accu);
    }
  #else
  inline uintD shift1right_loop_up (uintD* ptr, uintC count, uintD carry)
    { carry = carry << (intDsize-1); // carry zu einem einzigen Bit machen
      dotimesC(count,count,
        { var uintD accu = *ptr;
          *ptr++ = (accu >> 1) | carry;
          carry = accu << (intDsize-1);
        });
      return carry;
    }
  #endif

// Schiebeschleife um i Bits nach rechts:
// übertrag = shiftright_loop_up(ptr,count,i);
// schiebt count (uintC>=0) Digits aufwärts von ptr um i Bits (0<i<intDsize)
// nach rechts, wobei links Nullen eingeschoben werden,
// und liefert den Übertrag (was rechts rauskommt, als Bits intDsize-1..intDsize-i).
  #if HAVE_DD
  inline uintD shiftright_loop_up (uintD* ptr, uintC count, uintC i)
    { var uintDD accu = 0;
      dotimesC(count,count,
        { // Die oberen i Bits von (uintD)accu bilden hier den Übertrag.
          accu = highlowDD_0(lowD(accu));
          // Die oberen i Bits von (uintDD)accu bilden hier den Übertrag.
          accu = (highlowDD_0(*ptr)>>i)+accu; *ptr++ = highD(accu);
        });
      return lowD(accu);
    }
  #else
  inline uintD shiftright_loop_up (uintD* ptr, uintC count, uintC i)
    { var uintC j = intDsize-i;
      var uintD carry = 0;
      dotimesC(count,count,
        { var uintD accu = *ptr;
          *ptr++ = (accu >> i) | carry;
          carry = accu << j;
        });
      return carry;
    }
  #endif

// Schiebeschleife um i Bits nach rechts:
// übertrag = shiftrightsigned_loop_up(ptr,count,i);
// schiebt count (uintC>0) Digits aufwärts von ptr um i Bits (0<i<intDsize)
// nach rechts, wobei links das MSBit ver-i-facht wird,
// und liefert den Übertrag (was rechts rauskommt, als Bits intDsize-1..intDsize-i).
  #if HAVE_DD
  inline uintD shiftrightsigned_loop_up (uintD* ptr, uintC count, uintC i)
    { var uintDD accu = // Übertrag mit i Vorzeichenbits initialisieren
                           highlowDD_0(sign_of_sintD((sintD)(*ptr)))>>i;
      dotimespC(count,count,
        { // Die oberen i Bits von (uintD)accu bilden hier den Übertrag.
          accu = highlowDD_0(lowD(accu));
          // Die oberen i Bits von (uintDD)accu bilden hier den Übertrag.
          accu = (highlowDD_0(*ptr)>>i)+accu; *ptr++ = highD(accu);
        });
      return lowD(accu);
    }
  #else
  inline uintD shiftrightsigned_loop_up (uintD* ptr, uintC count, uintC i)
    { var uintC j = intDsize-i;
      var uintD carry;
      { var uintD accu = *ptr;
        *ptr++ = (sintD)accu >> i;
        carry = accu << j;
        count--;
      }
      dotimesC(count,count,
        { var uintD accu = *ptr;
          *ptr++ = (accu >> i) | carry;
          carry = accu << j;
        });
      return carry;
    }
  #endif

// Schiebe- und Kopier-Schleife um i Bits nach rechts:
// übertrag = shiftrightcopy_loop_up(sourceptr,destptr,count,i,carry);
// kopiert count (uintC>=0) Digits aufwärts von sourceptr nach destptr
// und schiebt sie dabei um i Bits (0<i<intDsize) nach rechts, wobei carry
// (sozusagen als sourceptr[-1]) die i Bits ganz links bestimmt,
// und liefert den Übertrag (was rechts rauskommt, als Bits intDsize-1..intDsize-i).
  #if HAVE_DD
  inline uintD shiftrightcopy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry)
    { var uintDD accu = // Übertrag mit carry initialisieren
                           highlowDD_0(carry)>>i;
      dotimesC(count,count,
        { // Die oberen i Bits von (uintD)accu bilden hier den Übertrag.
          accu = highlowDD_0(lowD(accu));
          // Die oberen i Bits von (uintDD)accu bilden hier den Übertrag.
          accu = (highlowDD_0(*sourceptr++)>>i)+accu; *destptr++ = highD(accu);
        });
      return lowD(accu);
    }
  #else
  inline uintD shiftrightcopy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry)
    { var uintC j = intDsize-i;
      carry = carry << j;
      dotimesC(count,count,
        { var uintD accu = *sourceptr++;
          *destptr++ = (accu >> i) | carry;
          carry = accu << j;
        });
      return carry;
    }
  #endif

#endif

#ifndef MUL_LOOPS

// Multiplikations-Einfachschleife:
// Multipliziert eine UDS mit einem kleinen Digit und addiert ein kleines Digit.
// mulusmall_loop_down(digit,ptr,len,newdigit)
// multipliziert die UDS  ptr[-len..-1]  mit digit (>=2, <=36),
// addiert dabei newdigit (>=0, <digit) zur letzten Ziffer,
// und liefert den Carry (>=0, <digit).
  #if HAVE_DD
  inline uintD mulusmall_loop_down (uintD digit, uintD* ptr, uintC len, uintD newdigit)
    { var uintDD carry = newdigit;
      dotimesC(len,len,
        { // Hier ist 0 <= carry < digit.
          carry = carry + muluD(digit,*--ptr);
          // Hier ist 0 <= carry < 2^intDsize*digit.
          *ptr = lowD(carry);
          carry = (uintDD)highD(carry); // carry := floor(carry/2^intDsize) < digit
        });
      return lowD(carry);
    }
  #else
  inline uintD mulusmall_loop_down (uintD digit, uintD* ptr, uintC len, uintD newdigit)
    { var uintD carry = newdigit;
      dotimesC(len,len,
        { // Hier ist 0 <= carry < digit.
          var uintD hi;
          var uintD lo;
          muluD(digit,*--ptr,hi=,lo=);
          // Hier ist 0 <= 2^intDsize*hi + lo + carry < 2^intDsize*digit.
          lo += carry; if (lo < carry) { hi += 1; }
          *ptr = lo;
          carry = hi;
        });
      return carry;
    }
  #endif

// Multiplikations-Einfachschleife:
// Multipliziert eine UDS mit einem Digit und legt das Ergebnis in einer
// zweiten UDS ab.
// mulu_loop_down(digit,sourceptr,destptr,len);
// multipliziert die UDS  sourceptr[-len..-1]  (len>0)
// mit dem einzelnen  digit
// und legt das Ergebnis in der UDS  destptr[-len-1..-1]  ab.
  #if HAVE_DD
  inline void mulu_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintDD carry = 0;
      dotimespC(len,len,
        { // Hier ist carry=digit=0 oder 0 <= carry < digit.
          carry = carry + muluD(digit,*--sourceptr);
          // Hier ist carry=digit=0 oder 0 <= carry < 2^intDsize*digit.
          *--destptr = lowD(carry);
          carry = (uintDD)highD(carry); // carry := floor(carry/2^intDsize) < digit
        });
      *--destptr = lowD(carry);
    }
  #else
  inline void mulu_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintD carry = 0;
      dotimespC(len,len,
        { // Hier ist carry=digit=0 oder 0 <= carry < digit.
          var uintD hi;
          var uintD lo;
          muluD(digit,*--sourceptr,hi=,lo=);
          // Hier ist 0 <= 2^intDsize*hi + lo + carry < 2^intDsize*digit oder hi=lo=carry=digit=0.
          lo += carry; if (lo < carry) { hi += 1; }
          *--destptr = lo;
          carry = hi;
        });
      *--destptr = carry;
    }
  #endif

// Multiplikations-Einfachschleife mit Akkumulation:
// Multipliziert eine UDS mit einem Digit und addiert das Ergebnis zu einer
// zweiten UDS auf.
// muluadd_loop_down(digit,sourceptr,destptr,len);
// multipliziert die UDS  sourceptr[-len..-1]  (len>0)
// mit dem einzelnen digit, legt das Ergebnis in der UDS  destptr[-len..-1]
// ab und liefert den weiteren Übertrag.
  #if HAVE_DD
  inline uintD muluadd_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintDD carry = 0;
      if (!(digit==0))
        { dotimespC(len,len,
            { // Hier ist 0 <= carry <= digit.
              carry = carry + muluD(digit,*--sourceptr) + (uintDD)*--destptr;
              // Hier ist 0 <= carry <= 2^intDsize*digit + 2^intDsize-1.
              *destptr = lowD(carry);
              carry = (uintDD)highD(carry); // carry := floor(carry/2^intDsize) <= digit
            });
        }
      return lowD(carry);
    }
  #else
  inline uintD muluadd_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintD carry = 0;
      if (!(digit==0))
        { dotimespC(len,len,
            { // Hier ist 0 <= carry <= digit.
              var uintD hi;
              var uintD lo;
              muluD(digit,*--sourceptr,hi=,lo=);
              // Hier ist 0 <= 2^intDsize*hi + lo + carry + *--destptr <= 2^intDsize*digit+2^intDsize-1.
              lo += carry; if (lo < carry) { hi += 1; }
              carry = *--destptr;
              lo += carry; if (lo < carry) { hi += 1; }
              *destptr = lo;
              carry = hi;
            });
        }
      return carry;
    }
  #endif

// Multiplikations-Einfachschleife mit Diminution:
// Multipliziert eine UDS mit einem Digit und subtrahiert das Ergebnis von
// einer zweiten UDS.
// mulusub_loop_down(digit,sourceptr,destptr,len);
// multipliziert die UDS  sourceptr[-len..-1]  (len>0)  mit dem einzelnen
// digit, subtrahiert das Ergebnis von der UDS  destptr[-len..-1]  und liefert
// den weiteren Übertrag (>=0, evtl. von destptr[-len-1] zu subtrahieren).
  #if HAVE_DD
  inline uintD mulusub_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintDD carry = 0;
      if (!(digit==0))
        { dotimespC(len,len,
            { // Hier ist 0 <= carry <= digit.
              carry = carry + muluD(digit,*--sourceptr) + (uintD)(~(*--destptr));
              // Hier ist 0 <= carry <= 2^intDsize*digit + 2^intDsize-1.
              *destptr = ~lowD(carry);
              carry = (uintDD)highD(carry); // carry := floor(carry/2^intDsize) <= digit
              // Hier ist 0 <= carry <= digit.
            });
          return lowD(carry);
        }
        else
        return 0; // nichts zu subtrahieren -> kein Übertrag
    }
  #else
  inline uintD mulusub_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintD carry = 0;
      if (!(digit==0))
        { dotimespC(len,len,
            { // Hier ist 0 <= carry <= digit.
              var uintD hi;
              var uintD lo;
              muluD(digit,*--sourceptr,hi=,lo=);
              // Hier ist 0 <= 2^intDsize*hi + lo + carry + ~(*--destptr) <= 2^intDsize*digit+2^intDsize-1.
              lo += carry; if (lo < carry) { hi += 1; }
              carry = *--destptr;
              *destptr = carry - lo; if (carry < lo) { hi += 1; }
              carry = hi;
            });
          return carry;
        }
        else
        return 0; // nichts zu subtrahieren -> kein Übertrag
    }
  #endif

#endif

#ifndef DIV_LOOPS

// Divisions-Einfachschleife:
// Dividiert eine UDS durch ein Digit.
// divu_loop_up(digit,ptr,len)
// dividiert die UDS  ptr[0..len-1] durch digit,
// legt das Ergebnis in derselben UDS ab, und liefert den Rest (>=0, <digit).
  #if HAVE_DD
  inline uintD divu_loop_up (uintD digit, uintD* ptr, uintC len)
    { var uintD rest = 0;
      dotimesC(len,len,
        { divuD(highlowDD(rest,*ptr),digit,*ptr =, rest =); ptr++; }
        );
      return rest;
    }
  #else
  inline uintD divu_loop_up (uintD digit, uintD* ptr, uintC len)
    { var uintD rest = 0;
      dotimesC(len,len,
        { divuD(rest,*ptr,digit,*ptr =, rest =); ptr++; }
        );
      return rest;
    }
  #endif

// Divisions-Einfachschleife:
// Dividiert eine UDS durch ein Digit und legt das Ergebnis in einer
// zweiten UDS ab.
// divucopy_loop_up(digit,sourceptr,destptr,len)
// dividiert die UDS  sourceptr[0..len-1]  durch digit,
// legt das Ergebnis in der UDS  destptr[0..len-1]  ab,
// und liefert den Rest (>=0, <digit).
  #if HAVE_DD
  inline uintD divucopy_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintD rest = 0;
      dotimesC(len,len,
        { divuD(highlowDD(rest,*sourceptr++),digit,*destptr++ =, rest =); }
        );
      return rest;
    }
  #else
  inline uintD divucopy_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintD rest = 0;
      dotimesC(len,len,
        { divuD(rest,*sourceptr++,digit,*destptr++ =, rest =); }
        );
      return rest;
    }
  #endif

#endif

#else // !CL_DS_BIG_ENDIAN_P

#ifndef LOG_LOOPS

// OR-Schleife:
// or_loop_down(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits abwärts ab xptr und ab yptr
// mit Ziel ab xptr durch OR.
  inline void or_loop_down (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *--xptr |= *--yptr; } ); }

// XOR-Schleife:
// xor_loop_down(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits abwärts ab xptr und ab yptr
// mit Ziel ab xptr durch XOR.
  inline void xor_loop_down (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *--xptr ^= *--yptr; } ); }

// AND-Schleife:
// and_loop_down(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits abwärts ab xptr und ab yptr
// mit Ziel ab xptr durch AND.
  inline void and_loop_down (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *--xptr &= *--yptr; } ); }

// EQV-Schleife:
// eqv_loop_down(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits abwärts ab xptr und ab yptr
// mit Ziel ab xptr durch EQV (NOT XOR).
  inline void eqv_loop_down (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count,
        {var uintD temp = ~ (*--xptr ^ *--yptr); *xptr = temp; }
        );
    }

// NAND-Schleife:
// nand_loop_down(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits abwärts ab xptr und ab yptr
// mit Ziel ab xptr durch NAND (NOT AND).
  inline void nand_loop_down (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count,
        {var uintD temp = ~ (*--xptr & *--yptr); *xptr = temp; }
        );
    }

// NOR-Schleife:
// nor_loop_down(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits abwärts ab xptr und ab yptr
// mit Ziel ab xptr durch NOR (NOT OR).
  inline void nor_loop_down (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count,
        {var uintD temp = ~ (*--xptr | *--yptr); *xptr = temp; }
        );
    }

// ANDC2-Schleife:
// andc2_loop_down(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits abwärts ab xptr und ab yptr
// mit Ziel ab xptr durch ANDC2 (AND NOT).
  inline void andc2_loop_down (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *--xptr &= ~(*--yptr); } ); }

// ORC2-Schleife:
// orc2_loop_down(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits abwärts ab xptr und ab yptr
// mit Ziel ab xptr durch ORC2 (OR NOT).
  inline void orc2_loop_down (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *--xptr |= ~(*--yptr); } ); }

// NOT-Schleife:
// not_loop_down(xptr,count);
// verknüpft count (uintC>0) Digits abwärts ab xptr mit Ziel ab xptr
// durch NOT.
  inline void not_loop_down (uintD* xptr, uintC count)
    { dotimespC(count,count,
        {var uintD temp = ~ (*--xptr); *xptr = temp; }
        );
    }

#endif

#ifndef TEST_LOOPS

// AND-Test-Schleife:
// and_test_loop_down(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits abwärts ab xptr und ab yptr durch AND
// und testet, ob sich dabei ein Digit /=0 ergibt. Ergebnis true, falls ja.
  inline bool and_test_loop_down (const uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { if (*--xptr & *--yptr) return true; } );
      return false;
    }

// Vergleichsschleife:
// result = compare_loop_down(xptr,yptr,count);
// vergleicht nacheinander xptr[-1] mit yptr[-1], xptr[-2] mit yptr[-2], usw.,
// insgesamt count Digits, und liefert 0 falls alle gleich sind,
// +1 falls zuerst ein xptr[i]>yptr[i] ist,
// -1 falls zuerst ein xptr[i]<yptr[i] ist.
  inline cl_signean compare_loop_down (const uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count,
        { if (!(*--xptr == *--yptr))
            // verschiedene Digits gefunden
            return (*xptr > *yptr ? signean_plus : signean_minus);
        });
      return signean_null; // alle Digits gleich
    }

#endif

#ifndef ADDSUB_LOOPS

// Additionsschleife:
// übertrag = add_loop_up(sourceptr1,sourceptr2,destptr,count);
// addiert count (uintC>=0) Digits aufwärts von sourceptr1, von sourceptr2
// aufwärts nach destptr und liefert den Übertrag (0 oder /=0, was 1 bedeutet).
  inline uintD add_loop_up (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count)
    { var uintD source1;
      var uintD source2;
      if (!(count==0))
      do { source1 = *sourceptr1++;
           source2 = *sourceptr2++;
           *destptr++ = source1 + source2;
           if (source1 > (uintD)(~source2)) goto carry_1;
           carry_0:
           count--;
         }
         until (count==0);
      return 0;
      do { source1 = *sourceptr1++;
           source2 = *sourceptr2++;
           *destptr++ = source1 + source2 + 1;
           if (source1 < (uintD)(~source2)) goto carry_0;
           carry_1:
           count--;
         }
         until (count==0);
      return 1;
    }

// Additionsschleife:
// übertrag = addto_loop_up(sourceptr,destptr,count);
// addiert count (uintC>=0) Digits aufwärts von sourceptr, von destptr
// aufwärts nach destptr und liefert den Übertrag (0 oder /=0, was 1 bedeutet).
  inline uintD addto_loop_up (const uintD* sourceptr, uintD* destptr, uintC count)
    { var uintD source1;
      var uintD source2;
      if (!(count==0))
      do { source1 = *sourceptr++;
           source2 = *destptr;
           *destptr++ = source1 + source2;
           if (source1 > (uintD)(~source2)) goto carry_1;
           carry_0:
           count--;
         }
         until (count==0);
      return 0;
      do { source1 = *sourceptr++;
           source2 = *destptr;
           *destptr++ = source1 + source2 + 1;
           if (source1 < (uintD)(~source2)) goto carry_0;
           carry_1:
           count--;
         }
         until (count==0);
      return 1;
    }

// Incrementierschleife:
// übertrag = inc_loop_up(ptr,count);
// incrementiert count (uintC>=0) Digits aufwärts von ptr, so lange bis kein
// Übertrag mehr auftritt und liefert den Übertrag (0 oder /=0, was 1 bedeutet).
  inline uintD inc_loop_up (uintD* ptr, uintC count)
    { dotimesC(count,count,
        { if (!( ++(*ptr++) == 0 )) return 0; } // kein weiterer Übertrag
        );
      return 1; // weiterer Übertrag
    }

// Subtraktionsschleife:
// übertrag = sub_loop_up(sourceptr1,sourceptr2,destptr,count);
// subtrahiert count (uintC>=0) Digits aufwärts von sourceptr1, von sourceptr2
// aufwärts nach destptr und liefert den Übertrag (0 oder /=0, was -1 bedeutet).
  inline uintD sub_loop_up (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count)
    { var uintD source1;
      var uintD source2;
      if (!(count==0))
      do { source1 = *sourceptr1++;
           source2 = *sourceptr2++;
           *destptr++ = source1 - source2;
           if (source1 < source2) goto carry_1;
           carry_0:
           count--;
         }
         until (count==0);
      return 0;
      do { source1 = *sourceptr1++;
           source2 = *sourceptr2++;
           *destptr++ = source1 - source2 - 1;
           if (source1 > source2) goto carry_0;
           carry_1:
           count--;
         }
         until (count==0);
      return (uintD)(-1);
    }

// Subtraktionsschleife:
// übertrag = subx_loop_up(sourceptr1,sourceptr2,destptr,count,carry);
// subtrahiert count (uintC>=0) Digits aufwärts von sourceptr1 und addiert
// einen Carry (0 oder -1), von sourceptr2 aufwärts nach destptr und
// liefert den Übertrag (0 oder /=0, was -1 bedeutet).
  inline uintD subx_loop_up (const uintD* sourceptr1, const uintD* sourceptr2, uintD* destptr, uintC count, uintD carry)
    { var uintD source1;
      var uintD source2;
      if (carry==0)
        { if (!(count==0))
            do { source1 = *sourceptr1++;
                 source2 = *sourceptr2++;
                 *destptr++ = source1 - source2;
                 if (source1 < source2) goto carry_1;
                 carry_0:
                 count--;
               }
               until (count==0);
          return 0;
        }
        else
        { if (!(count==0))
            do { source1 = *sourceptr1++;
                 source2 = *sourceptr2++;
                 *destptr++ = source1 - source2 - 1;
                 if (source1 > source2) goto carry_0;
                 carry_1:
                 count--;
               }
               until (count==0);
          return (uintD)(-1);
    }   }

// Subtraktionsschleife:
// übertrag = subfrom_loop_up(sourceptr,destptr,count);
// subtrahiert count (uintC>=0) Digits aufwärts von sourceptr, von destptr
// aufwärts nach destptr (dest := dest - source)
// und liefert den Übertrag (0 oder /=0, was -1 bedeutet).
  inline uintD subfrom_loop_up (const uintD* sourceptr, uintD* destptr, uintC count)
    { var uintD source1;
      var uintD source2;
      if (!(count==0))
      do { source1 = *destptr;
           source2 = *sourceptr++;
           *destptr++ = source1 - source2;
           if (source1 < source2) goto carry_1;
           carry_0:
           count--;
         }
         until (count==0);
      return 0;
      do { source1 = *destptr;
           source2 = *sourceptr++;
           *destptr++ = source1 - source2 - 1;
           if (source1 > source2) goto carry_0;
           carry_1:
           count--;
         }
         until (count==0);
      return (uintD)(-1);
    }

// Decrementierschleife:
// übertrag = dec_loop_up(ptr,count);
// decrementiert count (uintC>=0) Digits aufwärts von ptr, so lange bis kein
// Übertrag mehr auftritt und liefert den Übertrag (0 oder -1).
  inline uintD dec_loop_up (uintD* ptr, uintC count)
    { dotimesC(count,count,
        { if (!( (*ptr++)-- == 0 )) return 0; } // kein weiterer Übertrag
        );
      return (uintD)(-1); // weiterer Übertrag
    }

// Negierschleife:
// übertrag = neg_loop_up(ptr,count);
// negiert count (uintC>=0) Digits aufwärts von ptr,
// und liefert den Übertrag (0 oder -1).
  inline uintD neg_loop_up (uintD* ptr, uintC count)
    { // erstes Digit /=0 suchen:
      until (count==0) { if (!(*ptr == 0)) goto L1; ptr++; count--; }
      return 0;
      L1: // erstes Digit /=0 gefunden, ab jetzt gibt's Carrys
      *ptr = - *ptr; count--; // 1 Digit negieren
      dotimesC(count,count, { ptr++; *ptr = ~ *ptr; } ); // alle anderen Digits invertieren
      return (uintD)(-1);
    }

#endif

#ifndef SHIFT_LOOPS

// Schiebeschleife um 1 Bit nach links:
// übertrag = shift1left_loop_up(ptr,count);
// schiebt count (uintC>=0) Digits aufwärts von ptr um 1 Bit nach links,
// und liefert den Übertrag (0 oder /=0, was 1 bedeutet).
  #if HAVE_DD
  inline uintD shift1left_loop_up (uintD* ptr, uintC count)
    { var uintDD accu = 0;
      dotimesC(count,count,
        { accu = ((uintDD)(*ptr)<<1)+accu; *ptr++ = lowD(accu);
          accu = (uintDD)(highD(accu));
        });
      return (uintD)accu;
    }
  #else
  inline uintD shift1left_loop_up (uintD* ptr, uintC count)
    { var uintD carry = 0;
      dotimesC(count,count,
        { var uintD accu = *ptr;
          *ptr++ = (accu<<1) | carry;
          carry = accu>>(intDsize-1);
        });
      return carry;
    }
  #endif

// Schiebeschleife um i Bits nach links:
// übertrag = shiftleft_loop_up(ptr,count,i,übertrag_init);
// schiebt count (uintC>=0) Digits aufwärts von ptr um i Bits (0<i<intDsize)
// nach links, schiebt dabei die i Bits aus übertrag_init rechts rein,
// und liefert den Übertrag (was links rauskommt, >=0, <2^i).
  #if HAVE_DD
  inline uintD shiftleft_loop_up (uintD* ptr, uintC count, uintC i, uintD carry)
    { var uintDD accu = (uintDD)carry;
      dotimesC(count,count,
        { accu = ((uintDD)(*ptr)<<i)+accu; *ptr++ = lowD(accu);
          accu = (uintDD)(highD(accu));
        });
      return (uintD)accu;
    }
  #else
  inline uintD shiftleft_loop_up (uintD* ptr, uintC count, uintC i, uintD carry)
    { var uintC j = intDsize-i;
      dotimesC(count,count,
        { var uintD accu = *ptr;
          *ptr++ = (accu<<i) | carry;
          carry = accu>>j;
        });
      return carry;
    }
  #endif

// Schiebe- und Kopierschleife um i Bits nach links:
// übertrag = shiftleftcopy_loop_up(sourceptr,destptr,count,i);
// kopiert count (uintC>=0) Digits aufwärts von sourceptr nach destptr
// und schiebt sie dabei um i Bits (0<i<intDsize) nach links,
// wobei ganz rechts mit i Nullbits aufgefüllt wird,
// und liefert den Übertrag (was links rauskommt, >=0, <2^i).
  #if HAVE_DD
  inline uintD shiftleftcopy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count, uintC i)
    { var uintDD accu = 0;
      dotimesC(count,count,
        { accu = ((uintDD)(*sourceptr++)<<i)+accu; *destptr++ = lowD(accu);
          accu = (uintDD)(highD(accu));
        });
      return (uintD)accu;
    }
  #else
  inline uintD shiftleftcopy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count, uintC i)
    { var uintC j = intDsize-i;
      var uintD carry = 0;
      dotimesC(count,count,
        { var uintD accu = *sourceptr++;
          *destptr++ = (accu<<i) | carry;
          carry = accu>>j;
        });
      return carry;
    }
  #endif

// Schiebeschleife um 1 Bit nach rechts:
// übertrag = shift1right_loop_down(ptr,count,übertrag_init);
// schiebt count (uintC>=0) Digits abwärts von ptr um 1 Bit nach rechts,
// wobei links das Bit übertrag_init (sollte =0 oder =-1 sein) hineingeschoben
// wird, und liefert den Übertrag (0 oder /=0, was 1 bedeutet).
  #if HAVE_DD
  inline uintD shift1right_loop_down (uintD* ptr, uintC count, uintD carry)
    { var uintDD accu = (sintDD)(sintD)carry & ((uintDD)1 << (2*intDsize-1)); // 0 oder bit(2*intDsize-1)
      dotimesC(count,count,
        { accu = (highlowDD_0(*--ptr)>>1)+accu; *ptr = highD(accu);
          accu = highlowDD_0(lowD(accu));
        });
      return highD(accu);
    }
  #else
  inline uintD shift1right_loop_down (uintD* ptr, uintC count, uintD carry)
    { carry = carry << (intDsize-1); // carry zu einem einzigen Bit machen
      dotimesC(count,count,
        { var uintD accu = *--ptr;
          *ptr = (accu >> 1) | carry;
          carry = accu << (intDsize-1);
        });
      return carry;
    }
  #endif

// Schiebeschleife um i Bits nach rechts:
// übertrag = shiftright_loop_down(ptr,count,i);
// schiebt count (uintC>=0) Digits abwärts von ptr um i Bits (0<i<intDsize)
// nach rechts, wobei links Nullen eingeschoben werden,
// und liefert den Übertrag (was rechts rauskommt, als Bits intDsize-1..intDsize-i).
  #if HAVE_DD
  inline uintD shiftright_loop_down (uintD* ptr, uintC count, uintC i)
    { var uintDD accu = 0;
      dotimesC(count,count,
        { // Die oberen i Bits von (uintD)accu bilden hier den Übertrag.
          accu = highlowDD_0(lowD(accu));
          // Die oberen i Bits von (uintDD)accu bilden hier den Übertrag.
          accu = (highlowDD_0(*--ptr)>>i)+accu; *ptr = highD(accu);
        });
      return lowD(accu);
    }
  #else
  inline uintD shiftright_loop_down (uintD* ptr, uintC count, uintC i)
    { var uintC j = intDsize-i;
      var uintD carry = 0;
      dotimesC(count,count,
        { var uintD accu = *--ptr;
          *ptr = (accu >> i) | carry;
          carry = accu << j;
        });
      return carry;
    }
  #endif

// Schiebeschleife um i Bits nach rechts:
// übertrag = shiftrightsigned_loop_down(ptr,count,i);
// schiebt count (uintC>0) Digits abwärts von ptr um i Bits (0<i<intDsize)
// nach rechts, wobei links das MSBit ver-i-facht wird,
// und liefert den Übertrag (was rechts rauskommt, als Bits intDsize-1..intDsize-i).
  #if HAVE_DD
  inline uintD shiftrightsigned_loop_down (uintD* ptr, uintC count, uintC i)
    { var uintDD accu = // Übertrag mit i Vorzeichenbits initialisieren
                           highlowDD_0(sign_of_sintD((sintD)(ptr[-1])))>>i;
      dotimespC(count,count,
        { // Die oberen i Bits von (uintD)accu bilden hier den Übertrag.
          accu = highlowDD_0(lowD(accu));
          // Die oberen i Bits von (uintDD)accu bilden hier den Übertrag.
          accu = (highlowDD_0(*--ptr)>>i)+accu; *ptr = highD(accu);
        });
      return lowD(accu);
    }
  #else
  inline uintD shiftrightsigned_loop_down (uintD* ptr, uintC count, uintC i)
    { var uintC j = intDsize-i;
      var uintD carry;
      { var uintD accu = *--ptr;
        *ptr = (sintD)accu >> i;
        carry = accu << j;
        count--;
      }
      dotimesC(count,count,
        { var uintD accu = *--ptr;
          *ptr = (accu >> i) | carry;
          carry = accu << j;
        });
      return carry;
    }
  #endif

// Schiebe- und Kopier-Schleife um i Bits nach rechts:
// übertrag = shiftrightcopy_loop_down(sourceptr,destptr,count,i,carry);
// kopiert count (uintC>=0) Digits abwärts von sourceptr nach destptr
// und schiebt sie dabei um i Bits (0<i<intDsize) nach rechts, wobei carry
// (sozusagen als sourceptr[0]) die i Bits ganz links bestimmt,
// und liefert den Übertrag (was rechts rauskommt, als Bits intDsize-1..intDsize-i).
  #if HAVE_DD
  inline uintD shiftrightcopy_loop_down (const uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry)
    { var uintDD accu = // Übertrag mit carry initialisieren
                           highlowDD_0(carry)>>i;
      dotimesC(count,count,
        { // Die oberen i Bits von (uintD)accu bilden hier den Übertrag.
          accu = highlowDD_0(lowD(accu));
          // Die oberen i Bits von (uintDD)accu bilden hier den Übertrag.
          accu = (highlowDD_0(*--sourceptr)>>i)+accu; *--destptr = highD(accu);
        });
      return lowD(accu);
    }
  #else
  inline uintD shiftrightcopy_loop_down (const uintD* sourceptr, uintD* destptr, uintC count, uintC i, uintD carry)
    { var uintC j = intDsize-i;
      carry = carry << j;
      dotimesC(count,count,
        { var uintD accu = *--sourceptr;
          *--destptr = (accu >> i) | carry;
          carry = accu << j;
        });
      return carry;
    }
  #endif

#endif

#ifndef MUL_LOOPS

// Multiplikations-Einfachschleife:
// Multipliziert eine UDS mit einem kleinen Digit und addiert ein kleines Digit.
// mulusmall_loop_up(digit,ptr,len,newdigit)
// multipliziert die UDS  ptr[0..len-1]  mit digit (>=2, <=36),
// addiert dabei newdigit (>=0, <digit) zur letzten Ziffer,
// und liefert den Carry (>=0, <digit).
  #if HAVE_DD
  inline uintD mulusmall_loop_up (uintD digit, uintD* ptr, uintC len, uintD newdigit)
    { var uintDD carry = newdigit;
      dotimesC(len,len,
        { // Hier ist 0 <= carry < digit.
          carry = carry + muluD(digit,*ptr);
          // Hier ist 0 <= carry < 2^intDsize*digit.
          *ptr++ = lowD(carry);
          carry = (uintDD)highD(carry); // carry := floor(carry/2^intDsize) < digit
        });
      return lowD(carry);
    }
  #else
  inline uintD mulusmall_loop_up (uintD digit, uintD* ptr, uintC len, uintD newdigit)
    { var uintD carry = newdigit;
      dotimesC(len,len,
        { // Hier ist 0 <= carry < digit.
          var uintD hi;
          var uintD lo;
          muluD(digit,*ptr,hi=,lo=);
          // Hier ist 0 <= 2^intDsize*hi + lo + carry < 2^intDsize*digit.
          lo += carry; if (lo < carry) { hi += 1; }
          *ptr++ = lo;
          carry = hi;
        });
      return carry;
    }
  #endif

// Multiplikations-Einfachschleife:
// Multipliziert eine UDS mit einem Digit und legt das Ergebnis in einer
// zweiten UDS ab.
// mulu_loop_up(digit,sourceptr,destptr,len);
// multipliziert die UDS  sourceptr[0..len-1]  (len>0)
// mit dem einzelnen  digit
// und legt das Ergebnis in der UDS  destptr[0..len]  ab.
  #if HAVE_DD
  inline void mulu_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintDD carry = 0;
      dotimespC(len,len,
        { // Hier ist carry=digit=0 oder 0 <= carry < digit.
          carry = carry + muluD(digit,*sourceptr++);
          // Hier ist carry=digit=0 oder 0 <= carry < 2^intDsize*digit.
          *destptr++ = lowD(carry);
          carry = (uintDD)highD(carry); // carry := floor(carry/2^intDsize) < digit
        });
      *destptr++ = lowD(carry);
    }
  #else
  inline void mulu_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintD carry = 0;
      dotimespC(len,len,
        { // Hier ist carry=digit=0 oder 0 <= carry < digit.
          var uintD hi;
          var uintD lo;
          muluD(digit,*sourceptr++,hi=,lo=);
          // Hier ist 0 <= 2^intDsize*hi + lo + carry < 2^intDsize*digit oder hi=lo=carry=digit=0.
          lo += carry; if (lo < carry) { hi += 1; }
          *destptr++ = lo;
          carry = hi;
        });
      *destptr++ = carry;
    }
  #endif

// Multiplikations-Einfachschleife mit Akkumulation:
// Multipliziert eine UDS mit einem Digit und addiert das Ergebnis zu einer
// zweiten UDS auf.
// muluadd_loop_up(digit,sourceptr,destptr,len);
// multipliziert die UDS  sourceptr[0..len-1]  (len>0)
// mit dem einzelnen digit, legt das Ergebnis in der UDS  destptr[0..len-1]
// ab und liefert den weiteren Übertrag.
  #if HAVE_DD
  inline uintD muluadd_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintDD carry = 0;
      if (!(digit==0))
        { dotimespC(len,len,
            { // Hier ist 0 <= carry <= digit.
              carry = carry + muluD(digit,*sourceptr++) + (uintDD)*destptr;
              // Hier ist 0 <= carry <= 2^intDsize*digit + 2^intDsize-1.
              *destptr++ = lowD(carry);
              carry = (uintDD)highD(carry); // carry := floor(carry/2^intDsize) <= digit
            });
        }
      return lowD(carry);
    }
  #else
  inline uintD muluadd_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintD carry = 0;
      if (!(digit==0))
        { dotimespC(len,len,
            { // Hier ist 0 <= carry <= digit.
              var uintD hi;
              var uintD lo;
              muluD(digit,*sourceptr++,hi=,lo=);
              // Hier ist 0 <= 2^intDsize*hi + lo + carry + *destptr <= 2^intDsize*digit+2^intDsize-1.
              lo += carry; if (lo < carry) { hi += 1; }
              carry = *destptr;
              lo += carry; if (lo < carry) { hi += 1; }
              *destptr++ = lo;
              carry = hi;
            });
        }
      return carry;
    }
  #endif

// Multiplikations-Einfachschleife mit Diminution:
// Multipliziert eine UDS mit einem Digit und subtrahiert das Ergebnis von
// einer zweiten UDS.
// mulusub_loop_up(digit,sourceptr,destptr,len);
// multipliziert die UDS  sourceptr[0..len-1]  (len>0)  mit dem einzelnen
// digit, subtrahiert das Ergebnis von der UDS  destptr[0..len-1]  und liefert
// den weiteren Übertrag (>=0, evtl. von destptr[len] zu subtrahieren).
  #if HAVE_DD
  inline uintD mulusub_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintDD carry = 0;
      if (!(digit==0))
        { dotimespC(len,len,
            { // Hier ist 0 <= carry <= digit.
              carry = carry + muluD(digit,*sourceptr++) + (uintD)(~(*destptr));
              // Hier ist 0 <= carry <= 2^intDsize*digit + 2^intDsize-1.
              *destptr++ = ~lowD(carry);
              carry = (uintDD)highD(carry); // carry := floor(carry/2^intDsize) <= digit
              // Hier ist 0 <= carry <= digit.
            });
          return lowD(carry);
        }
        else
        return 0; // nichts zu subtrahieren -> kein Übertrag
    }
  #else
  inline uintD mulusub_loop_up (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintD carry = 0;
      if (!(digit==0))
        { dotimespC(len,len,
            { // Hier ist 0 <= carry <= digit.
              var uintD hi;
              var uintD lo;
              muluD(digit,*sourceptr++,hi=,lo=);
              // Hier ist 0 <= 2^intDsize*hi + lo + carry + ~(*destptr) <= 2^intDsize*digit+2^intDsize-1.
              lo += carry; if (lo < carry) { hi += 1; }
              carry = *destptr;
              *destptr++ = carry - lo; if (carry < lo) { hi += 1; }
              carry = hi;
            });
          return carry;
        }
        else
        return 0; // nichts zu subtrahieren -> kein Übertrag
    }
  #endif

#endif

#ifndef DIV_LOOPS

// Divisions-Einfachschleife:
// Dividiert eine UDS durch ein Digit.
// divu_loop_down(digit,ptr,len)
// dividiert die UDS  ptr[-len..-1] durch digit,
// legt das Ergebnis in derselben UDS ab, und liefert den Rest (>=0, <digit).
  #if HAVE_DD
  inline uintD divu_loop_down (uintD digit, uintD* ptr, uintC len)
    { var uintD rest = 0;
      dotimesC(len,len,
        { --ptr; divuD(highlowDD(rest,*ptr),digit,*ptr =, rest =); }
        );
      return rest;
    }
  #else
  inline uintD divu_loop_down (uintD digit, uintD* ptr, uintC len)
    { var uintD rest = 0;
      dotimesC(len,len,
        { --ptr; divuD(rest,*ptr,digit,*ptr =, rest =); }
        );
      return rest;
    }
  #endif

// Divisions-Einfachschleife:
// Dividiert eine UDS durch ein Digit und legt das Ergebnis in einer
// zweiten UDS ab.
// divucopy_loop_down(digit,sourceptr,destptr,len)
// dividiert die UDS  sourceptr[-len..-1]  durch digit,
// legt das Ergebnis in der UDS  destptr[-len..-1]  ab,
// und liefert den Rest (>=0, <digit).
  #if HAVE_DD
  inline uintD divucopy_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintD rest = 0;
      dotimesC(len,len,
        { divuD(highlowDD(rest,*--sourceptr),digit,*--destptr =, rest =); }
        );
      return rest;
    }
  #else
  inline uintD divucopy_loop_down (uintD digit, const uintD* sourceptr, uintD* destptr, uintC len)
    { var uintD rest = 0;
      dotimesC(len,len,
        { divuD(rest,*--sourceptr,digit,*--destptr =, rest =); }
        );
      return rest;
    }
  #endif

#endif

#endif // !CL_DS_BIG_ENDIAN_P

#if !defined(TEST_LOOPS) && !CL_DS_BIG_ENDIAN_P

// Vergleichsschleife:
// result = compare_loop_up(xptr,yptr,count);
// vergleicht nacheinander xptr[0] mit yptr[0], xptr[1] mit yptr[1], usw.,
// insgesamt count Digits, und liefert 0 falls alle gleich sind,
// +1 falls zuerst ein xptr[i]>yptr[i] ist,
// -1 falls zuerst ein xptr[i]<yptr[i] ist.
  inline cl_signean compare_loop_up (const uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count,
        { if (!(*xptr++ == *yptr++))
            // verschiedene Digits gefunden
            return (*--xptr > *--yptr ? signean_plus : signean_minus);
        });
      return signean_null; // alle Digits gleich
    }

#endif

#if !defined(LOG_LOOPS) && !CL_DS_BIG_ENDIAN_P

// XOR-Schleife:
// xor_loop_up(xptr,yptr,count);
// verknüpft count (uintC>=0) Digits aufwärts ab xptr und ab yptr
// mit Ziel ab xptr durch XOR.
  inline void xor_loop_up (uintD* xptr, const uintD* yptr, uintC count)
    { dotimesC(count,count, { *xptr++ ^= *yptr++; } ); }

#endif

#if !defined(SHIFT_LOOPS) && CL_DS_BIG_ENDIAN_P

// Schiebe- und Kopierschleife um i Bits nach links:
// übertrag = shiftleftcopy_loop_up(sourceptr,destptr,count,i);
// kopiert count (uintC>=0) Digits aufwärts von sourceptr nach destptr
// und schiebt sie dabei um i Bits (0<i<intDsize) nach links,
// wobei ganz rechts mit i Nullbits aufgefüllt wird,
// und liefert den Übertrag (was links rauskommt, >=0, <2^i).
  #if HAVE_DD
  inline uintD shiftleftcopy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count, uintC i)
    { var uintDD accu = 0;
      dotimesC(count,count,
        { accu = ((uintDD)(*sourceptr++)<<i)+accu; *destptr++ = lowD(accu);
          accu = (uintDD)(highD(accu));
        });
      return (uintD)accu;
    }
  #else
  inline uintD shiftleftcopy_loop_up (const uintD* sourceptr, uintD* destptr, uintC count, uintC i)
    { var uintC j = intDsize-i;
      var uintD carry = 0;
      dotimesC(count,count,
        { var uintD accu = *sourceptr++;
          *destptr++ = (accu<<i) | carry;
          carry = accu>>j;
        });
      return carry;
    }
  #endif

#endif

#if !defined(SHIFT_LOOPS)

// Schiebe- und XOR-Schleife:
// shiftxor_loop_up(xptr,yptr,count,i);
// verknüpft count+1 Digits aufwärts ab xptr mit count Digits aufwärts ab yptr,
// um i Bits verschoben, durch XOR. (count uintC>=0, 0<i<intDsize)
  #if HAVE_DD
  inline void shiftxor_loop_up (uintD* xptr, const uintD* yptr, uintC count, uintC i)
    { if (count > 0)
        { var uintD carry = xptr[0];
          dotimespC(count,count,
            { var uintDD accu = highlowDD(xptr[1],carry);
              accu = ((uintDD)(*yptr++)<<i) ^ accu;
              *xptr++ = lowD(accu);
              carry = highD(accu);
            });
          *xptr = carry;
    }   }
  #else
  inline void shiftxor_loop_up (uintD* xptr, const uintD* yptr, uintC count, uintC i)
    { if (count > 0)
        { var uintC j = intDsize-i;
          var uintD carry = *xptr;
          dotimespC(count,count,
            { var uintD accu = *yptr++;
              *xptr++ = (accu<<i) ^ carry;
              carry = (accu>>j) ^ *xptr;
            });
          *xptr = carry;
    }   }
  #endif

#endif


// Endianness independent names for these functions.
#if CL_DS_BIG_ENDIAN_P
  #define copy_loop_msp              copy_loop_up
  #define copy_loop_lsp              copy_loop_down
  #define fill_loop_msp              fill_loop_up
  #define fill_loop_lsp              fill_loop_down
  #define clear_loop_msp             clear_loop_up
  #define clear_loop_lsp             clear_loop_down
  #define test_loop_msp              test_loop_up
  #define or_loop_msp                or_loop_up
  #define xor_loop_msp               xor_loop_up
  #define and_loop_msp               and_loop_up
  #define eqv_loop_msp               eqv_loop_up
  #define nand_loop_msp              nand_loop_up
  #define nor_loop_msp               nor_loop_up
  #define andc2_loop_msp             andc2_loop_up
  #define orc2_loop_msp              orc2_loop_up
  #define not_loop_msp               not_loop_up
  #define and_test_loop_msp          and_test_loop_up
  #define compare_loop_msp           compare_loop_up
  #define add_loop_lsp               add_loop_down
  #define addto_loop_lsp             addto_loop_down
  #define inc_loop_lsp               inc_loop_down
  #define sub_loop_lsp               sub_loop_down
  #define subx_loop_lsp              subx_loop_down
  #define subfrom_loop_lsp           subfrom_loop_down
  #define dec_loop_lsp               dec_loop_down
  #define neg_loop_lsp               neg_loop_down
  #define shift1left_loop_lsp        shift1left_loop_down
  #define shiftleft_loop_lsp         shiftleft_loop_down
  #define shiftleftcopy_loop_lsp     shiftleftcopy_loop_down
  #define shift1right_loop_msp       shift1right_loop_up
  #define shiftright_loop_msp        shiftright_loop_up
  #define shiftrightsigned_loop_msp  shiftrightsigned_loop_up
  #define shiftrightcopy_loop_msp    shiftrightcopy_loop_up
  #define mulusmall_loop_lsp         mulusmall_loop_down
  #define mulu_loop_lsp              mulu_loop_down
  #define muluadd_loop_lsp           muluadd_loop_down
  #define mulusub_loop_lsp           mulusub_loop_down
  #define divu_loop_msp              divu_loop_up
  #define divucopy_loop_msp          divucopy_loop_up
#else
  #define copy_loop_msp              copy_loop_down
  #define copy_loop_lsp              copy_loop_up
  #define fill_loop_msp              fill_loop_down
  #define fill_loop_lsp              fill_loop_up
  #define clear_loop_msp             clear_loop_down
  #define clear_loop_lsp             clear_loop_up
  #define test_loop_msp              test_loop_down
  #define or_loop_msp                or_loop_down
  #define xor_loop_msp               xor_loop_down
  #define and_loop_msp               and_loop_down
  #define eqv_loop_msp               eqv_loop_down
  #define nand_loop_msp              nand_loop_down
  #define nor_loop_msp               nor_loop_down
  #define andc2_loop_msp             andc2_loop_down
  #define orc2_loop_msp              orc2_loop_down
  #define not_loop_msp               not_loop_down
  #define and_test_loop_msp          and_test_loop_down
  #define compare_loop_msp           compare_loop_down
  #define add_loop_lsp               add_loop_up
  #define addto_loop_lsp             addto_loop_up
  #define inc_loop_lsp               inc_loop_up
  #define sub_loop_lsp               sub_loop_up
  #define subx_loop_lsp              subx_loop_up
  #define subfrom_loop_lsp           subfrom_loop_up
  #define dec_loop_lsp               dec_loop_up
  #define neg_loop_lsp               neg_loop_up
  #define shift1left_loop_lsp        shift1left_loop_up
  #define shiftleft_loop_lsp         shiftleft_loop_up
  #define shiftleftcopy_loop_lsp     shiftleftcopy_loop_up
  #define shift1right_loop_msp       shift1right_loop_down
  #define shiftright_loop_msp        shiftright_loop_down
  #define shiftrightsigned_loop_msp  shiftrightsigned_loop_down
  #define shiftrightcopy_loop_msp    shiftrightcopy_loop_down
  #define mulusmall_loop_lsp         mulusmall_loop_up
  #define mulu_loop_lsp              mulu_loop_up
  #define muluadd_loop_lsp           muluadd_loop_up
  #define mulusub_loop_lsp           mulusub_loop_up
  #define divu_loop_msp              divu_loop_down
  #define divucopy_loop_msp          divucopy_loop_down
#endif

// Endianness independent loops where the direction doesn't matter.
#if CL_DS_BIG_ENDIAN_P
  #define DS_clear_loop(MSDptr,len,LSDptr)  (void)clear_loop_up(MSDptr,len)
  #define DS_test_loop(MSDptr,len,LSDptr)  test_loop_up(MSDptr,len)
#else
  #define DS_clear_loop(MSDptr,len,LSDptr)  (void)clear_loop_up(LSDptr,len)
  #define DS_test_loop(MSDptr,len,LSDptr)  test_loop_up(LSDptr,len)
#endif


// Umwandlungsroutinen Digit-Sequence-Teil <--> Longword:

// get_32_Dptr(ptr)
//   holt die nächsten 32 Bits aus den 32/intDsize Digits ab ptr.
// set_32_Dptr(ptr,wert);
//   speichert den Wert wert (32 Bits) in die 32/intDsize Digits ab ptr.
// get_max32_Dptr(count,ptr)
//   holt die nächsten count Bits aus den ceiling(count/intDsize) Digits ab ptr.
// set_max32_Dptr(count,ptr,wert)
//   speichert wert (count Bits) in die ceiling(count/intDsize) Digits ab ptr.
// Jeweils ptr eine Variable vom Typ uintD*,
//         wert eine Variable vom Typ uint32,
//         count eine Variable oder constant-expression mit Wert >=0, <=32.
  #if (intDsize==32)
    inline uint32 get_32_Dptr (const uintD* ptr)
    {
	return mspref(ptr,0);
    }
    inline void set_32_Dptr (uintD* ptr, uint32 wert)
    {
	mspref(ptr,0) = wert;
    }
    inline uint32 get_max32_Dptr (uintC count, const uintD* ptr)
    {
	return count==0 ? 0 :
                          mspref(ptr,0);
    }
    inline void set_max32_Dptr (uintC count, uintD* ptr, uint32 wert)
    {
	if (count==0) return;
	mspref(ptr,0) = wert; return;
    }
  #endif
  #if (intDsize==16)
    inline uint32 get_32_Dptr (const uintD* ptr)
    {
	return ((uint32)mspref(ptr,0)<<16) | (uint32)mspref(ptr,1);
    }
    inline void set_32_Dptr (uintD* ptr, uint32 wert)
    {
	mspref(ptr,0) = (uintD)(wert>>16); mspref(ptr,1) = (uintD)wert;
    }
    inline uint32 get_max32_Dptr (uintC count, const uintD* ptr)
    {
	return count==0 ? 0 :
	       count<=16 ? mspref(ptr,0) :
	                   ((uint32)mspref(ptr,0)<<16) | (uint32)mspref(ptr,1);
    }
    inline void set_max32_Dptr (uintC count, uintD* ptr, uint32 wert)
    {
	if (count==0) return;
	if (count<=16) { mspref(ptr,0) = (uintD)wert; return; }
	mspref(ptr,0) = (uintD)(wert>>16); mspref(ptr,1) = (uintD)wert; return;
    }
  #endif
  #if (intDsize==8)
    inline uint32 get_32_Dptr (const uintD* ptr)
    {
	return ((((((uint32)mspref(ptr,0) <<8) | (uint32)mspref(ptr,1)) <<8) | (uint32)mspref(ptr,2)) <<8) | (uint32)mspref(ptr,3);
    }
    inline void set_32_Dptr (uintD* ptr, uint32 wert)
    {
	mspref(ptr,0) = (uintD)(wert>>24); mspref(ptr,1) = (uintD)(wert>>16); mspref(ptr,2) = (uintD)(wert>>8); mspref(ptr,3) = (uintD)wert;
    }
    inline uint32 get_max32_Dptr (uintC count, const uintD* ptr)
    {
	return count==0 ? 0 :
	       count<=8 ? mspref(ptr,0) :
	       count<=16 ? ((uint32)mspref(ptr,0)<<8) | (uint32)mspref(ptr,1) :
	       count<=24 ? ((((uint32)mspref(ptr,0)<<8) | (uint32)mspref(ptr,1))<<8) | (uint32)mspref(ptr,2) :
	                   ((((((uint32)mspref(ptr,0)<<8) | (uint32)mspref(ptr,1))<<8) | (uint32)mspref(ptr,2))<<8) | (uint32)mspref(ptr,3);
    }
    inline void set_max32_Dptr (uintC count, uintD* ptr, uint32 wert)
    {
	if (count==0) return;
	if (count<=8) { mspref(ptr,0) = (uintD)wert; return; }
	if (count<=16) { mspref(ptr,0) = (uintD)(wert>>8); mspref(ptr,1) = (uintD)wert; return; }
	if (count<=24) { mspref(ptr,0) = (uintD)(wert>>16); mspref(ptr,1) = (uintD)(wert>>8); mspref(ptr,2) = (uintD)wert; return; }
	mspref(ptr,0) = (uintD)(wert>>24); mspref(ptr,1) = (uintD)(wert>>16); mspref(ptr,2) = (uintD)(wert>>8); mspref(ptr,3) = (uintD)wert; return;
    }
  #endif

#if (cl_word_size==64)
// get_64_Dptr(ptr)
//   holt die nächsten 64 Bits aus den 64/intDsize Digits ab ptr.
// set_64_Dptr(ptr,wert);
//   speichert den Wert wert (64 Bits) in die 64/intDsize Digits ab ptr.
// get_max64_Dptr(count,ptr)
//   holt die nächsten count Bits aus den ceiling(count/intDsize) Digits ab ptr.
// set_max64_Dptr(count,ptr,wert)
//   speichert wert (count Bits) in die ceiling(count/intDsize) Digits ab ptr.
// Jeweils ptr eine Variable vom Typ uintD*,
//         wert eine Variable vom Typ uint64,
//         count eine Variable oder constant-expression mit Wert >=0, <=64.
  #if (intDsize==64)
    inline uint64 get_64_Dptr (const uintD* ptr)
    {
	return mspref(ptr,0);
    }
    inline void set_64_Dptr (uintD* ptr, uint64 wert)
    {
	mspref(ptr,0) = wert;
    }
    inline uint64 get_max64_Dptr (uintC count, const uintD* ptr)
    {
	return count==0 ? 0 : mspref(ptr,0);
    }
    inline void set_max64_Dptr (uintC count, uintD* ptr, uint64 wert)
    {
	if (count==0) return;
	mspref(ptr,0) = wert; return;
    }
  #else // (intDsize<=32)
    inline uint64 get_64_Dptr (const uintD* ptr)
    {
	return ((uint64)get_32_Dptr(ptr) << 32) | (uint64)get_32_Dptr(ptr mspop 32/intDsize);
    }
    inline void set_64_Dptr (uintD* ptr, uint64 wert)
    {
	set_32_Dptr(ptr,(uint32)(wert>>32));
	set_32_Dptr(ptr mspop 32/intDsize,(uint32)wert);
    }
    inline uint64 get_max64_Dptr (uintC count, const uintD* ptr)
    {
	return count==0 ? 0 :
               count<=32 ? (uint64)get_max32_Dptr(count,ptr) :
	                   ((uint64)get_max32_Dptr(count-32,ptr) << 32) | (uint64)get_32_Dptr(ptr mspop ceiling(count-32,intDsize));
    }
    inline void set_max64_Dptr (uintC count, uintD* ptr, uint64 wert)
    {
	if (count==0) return;
	if (count<=32) { set_max32_Dptr(count,ptr,(uint32)wert); return; }
	set_max32_Dptr(count-32,ptr,(uint32)(wert>>32));
	set_32_Dptr(ptr mspop ceiling(count-32,intDsize),(uint32)wert); return;
    }
  #endif
#endif

// get_uint1D_Dptr(ptr)  holt 1 Digit (unsigned) ab ptr
// get_uint2D_Dptr(ptr)  holt 2 Digits (unsigned) ab ptr
// get_uint3D_Dptr(ptr)  holt 3 Digits (unsigned) ab ptr
// get_uint4D_Dptr(ptr)  holt 4 Digits (unsigned) ab ptr
// get_sint1D_Dptr(ptr)  holt 1 Digit (signed) ab ptr
// get_sint2D_Dptr(ptr)  holt 2 Digits (signed) ab ptr
// get_sint3D_Dptr(ptr)  holt 3 Digits (signed) ab ptr
// get_sint4D_Dptr(ptr)  holt 4 Digits (signed) ab ptr
// Jeweils ptr eine Variable vom Typ uintD*.
// NB: Bei intDsize==64 sind diese Funktionen nur sehr bedingt tauglich.
  inline uint32 get_uint1D_Dptr (const uintD* ptr)
  {
	return lspref(ptr,0);
  }
  inline sint32 get_sint1D_Dptr (const uintD* ptr)
  {
	return (sint32)(sintD)lspref(ptr,0);
  }
  #if (intDsize < 32)
  inline uint32 get_uint2D_Dptr (const uintD* ptr)
  {
	return ((uint32)lspref(ptr,1) << intDsize) | (uint32)lspref(ptr,0);
  }
  inline sint32 get_sint2D_Dptr (const uintD* ptr)
  {
	return ((uint32)(sint32)(sintD)lspref(ptr,1) << intDsize) | (uint32)lspref(ptr,0);
  }
  #else
  #define get_uint2D_Dptr(ptr)  get_uint1D_Dptr(ptr)
  #define get_sint2D_Dptr(ptr)  (sint32)get_uint2D_Dptr(ptr)
  #endif
  #if (intDsize < 16)
  inline uint32 get_uint3D_Dptr (const uintD* ptr)
  {
	return ((((uint32)lspref(ptr,2) << intDsize) | (uint32)lspref(ptr,1)) << intDsize) | (uint32)lspref(ptr,0);
  }
  inline sint32 get_sint3D_Dptr (const uintD* ptr)
  {
	return ((((uint32)(sint32)(sintD)lspref(ptr,2) << intDsize) | (uint32)lspref(ptr,1)) << intDsize) | (uint32)lspref(ptr,0);
  }
  inline uint32 get_uint4D_Dptr (const uintD* ptr)
  {
	return ((((((uint32)lspref(ptr,3) << intDsize) | (uint32)lspref(ptr,2)) << intDsize) | (uint32)lspref(ptr,1)) << intDsize) | (uint32)lspref(ptr,0);
  }
  inline sint32 get_sint4D_Dptr (const uintD* ptr)
  {
	return ((((((uint32)(sint32)(sintD)lspref(ptr,3) << intDsize) | (uint32)lspref(ptr,2)) << intDsize) | (uint32)lspref(ptr,1)) << intDsize) | (uint32)lspref(ptr,0);
  }
  #else
  #define get_uint3D_Dptr(ptr)  get_uint2D_Dptr(ptr)
  #define get_sint3D_Dptr(ptr)  (sint32)get_uint3D_Dptr(ptr)
  #define get_uint4D_Dptr(ptr)  get_uint2D_Dptr(ptr)
  #define get_sint4D_Dptr(ptr)  (sint32)get_uint4D_Dptr(ptr)
  #endif


// NUM_STACK ist eine Art Zahlen-Stack-Pointer.
// Verwendung:
//   {CL_ALLOCA_STACK;
//    ...
//    num_stack_alloc(...);
//    ...
//    num_stack_array(...);
//    ...
//   }
// CL_ALLOCA_STACK rettet den aktuellen Wert von NUM_STACK.
// Dann darf beliebig oft mit num_stack_alloc/num_stack_array Platz auf dem
// Zahlen-Stack belegt werden.
// Beim Ende des Blocks wird NUM_STACK wieder auf den vorigen Wert gesetzt,
// und der Platz gilt als wieder freigegeben.
// In jeder C-Funktion sollte CL_ALLOCA_STACK nur einmal aufgerufen werden.
// Wegen eines GCC-Bugs sollten Funktionen, die diese Macros benutzen,
// nicht inline deklariert sein.

// num_stack_array(need, low_addr = , high_addr = );
// num_stack_small_array(need, low_addr = , high_addr = );
// belegt need Digits auf dem Zahlen-Stack und legt die untere Grenze des
// allozierten Bereichs in low_addr und die obere Grenze in high_addr ab.
// Jedes von beiden ist optional.

// num_stack_alloc(need, MSDptr = , LSDptr = );
// num_stack_small_alloc(need, MSDptr = , LSDptr = );
// belegt need Digits auf dem Zahlen-Stack und legt den MSDptr und den
// LSDptr ab. Jedes von beiden ist optional.

// num_stack_alloc_1(need, MSDptr = , LSDptr = );
// num_stack_small_alloc_1(need, MSDptr = , LSDptr = );
// wie num_stack_alloc, nur daß unterhalb von MSDptr noch ein Digit Platz
// zusätzlich belegt wird.

#define num_stack_array(need,low_zuweisung,high_zuweisung)  \
  {var uintC __need = (uintC)(need);				\
   var uintD* __array = cl_alloc_array(uintD,__need);		\
   unused (low_zuweisung &__array[0]); unused (high_zuweisung &__array[__need]); \
  }
#define num_stack_small_array(need,low_zuweisung,high_zuweisung)  \
  {var uintC __need = (uintC)(need);				\
   var uintD* __array = cl_small_alloc_array(uintD,__need);	\
   unused (low_zuweisung &__array[0]); unused (high_zuweisung &__array[__need]); \
  }
#if CL_DS_BIG_ENDIAN_P
  #define num_stack_alloc(need,MSDptr_zuweisung,LSDptr_zuweisung)  \
    num_stack_array(need,MSDptr_zuweisung,LSDptr_zuweisung)
  #define num_stack_small_alloc(need,MSDptr_zuweisung,LSDptr_zuweisung)  \
    num_stack_small_array(need,MSDptr_zuweisung,LSDptr_zuweisung)
  #define num_stack_alloc_1(need,MSDptr_zuweisung,LSDptr_zuweisung)  \
    num_stack_array((uintC)(need)+1,MSDptr_zuweisung 1 + ,LSDptr_zuweisung)
  #define num_stack_small_alloc_1(need,MSDptr_zuweisung,LSDptr_zuweisung)  \
    num_stack_small_array((uintC)(need)+1,MSDptr_zuweisung 1 + ,LSDptr_zuweisung)
#else
  #define num_stack_alloc(need,MSDptr_zuweisung,LSDptr_zuweisung)  \
    num_stack_array(need,LSDptr_zuweisung,MSDptr_zuweisung)
  #define num_stack_small_alloc(need,MSDptr_zuweisung,LSDptr_zuweisung)  \
    num_stack_small_array(need,LSDptr_zuweisung,MSDptr_zuweisung)
  #define num_stack_alloc_1(need,MSDptr_zuweisung,LSDptr_zuweisung)  \
    num_stack_array((uintC)(need)+1,LSDptr_zuweisung,MSDptr_zuweisung -1 + )
  #define num_stack_small_alloc_1(need,MSDptr_zuweisung,LSDptr_zuweisung)  \
    num_stack_small_array((uintC)(need)+1,LSDptr_zuweisung,MSDptr_zuweisung -1 + )
#endif


// Macro: In der DS MSDptr/len/LSDptr wird eine 1 unterhalb des Pointers ptr
// addiert. Unterhalb von MSDptr muß 1 Digit Platz sein.
// Dabei ist  ptr - MSDptr = count  und  0 < count <= len .
// Eventuell wird MSDptr erniedrigt und len erhöht.
  #define DS_1_plus(ptr,count)  \
    {var uintD* ptr_from_DS_1_plus = (ptr);				\
     var uintC count_from_DS_1_plus = (count);				\
     loop { if (--count_from_DS_1_plus==0) /* Zähler erniedrigen      */\
              { /* Beim Most Significant Digit angelangt              */\
                lsprefnext(ptr_from_DS_1_plus) += 1;			\
                /* jetzt ist ptr_from_DS_1_plus = MSDptr              */\
                if (mspref(ptr_from_DS_1_plus,0) == (uintD)bit(intDsize-1)) \
                  { /* 7FFF + 1 muß zu 00008000 werden:               */\
                    lsprefnext(MSDptr) = 0;				\
                    len++;						\
                  }							\
                break;							\
              }								\
            if (!((lsprefnext(ptr_from_DS_1_plus) += 1) == 0)) /* weiterincrementieren */\
              break; /* kein weiterer Übertrag -> Schleife abbrechen  */\
    }     }

// Macro: In der DS MSDptr/len/LSDptr wird eine 1 unterhalb des Pointers ptr
// subtrahiert. Unterhalb von MSDptr muß 1 Digit Platz sein.
// Dabei ist  ptr - MSDptr = count  und  0 < count <= len .
// Eventuell wird MSDptr erniedrigt und len erhöht.
  #define DS_minus1_plus(ptr,count)  \
    {var uintD* ptr_from_DS_minus1_plus = (ptr);			\
     var uintC count_from_DS_minus1_plus = (count);			\
     loop { if (--count_from_DS_minus1_plus==0) /* Zähler erniedrigen */\
              { /* Beim Most Significant Digit angelangt              */\
                lsprefnext(ptr_from_DS_minus1_plus) -= 1;		\
                /* jetzt ist ptr_from_DS_minus1_plus = MSDptr         */\
                if (mspref(ptr_from_DS_minus1_plus,0) == (uintD)bit(intDsize-1)-1) \
                  { /* 8000 - 1 muß zu FFFF7FFF werden:               */\
                    lsprefnext(MSDptr) = (uintD)(-1);			\
                    len++;						\
                  }							\
                break;							\
              }								\
            if (!((sintD)(lsprefnext(ptr_from_DS_minus1_plus) -= 1) == -1)) /* weiterdecrementieren */\
              break; /* kein weiterer Übertrag -> Schleife abbrechen  */\
    }     }


// Multiplikations-Doppelschleife:
// Multipliziert zwei UDS und legt das Ergebnis in einer dritten UDS ab.
// cl_UDS_mul(sourceptr1,len1,sourceptr2,len2,destptr);
// multipliziert die UDS  sourceptr1[-len1..-1]  (len1>0)
//           mit der UDS  sourceptr2[-len1..-1]  (len2>0)
// und legt das Ergebnis in der UDS  destptr[-len..-1]  (len=len1+len2) ab.
// Unterhalb von destptr werden len Digits Platz benötigt.
extern void cl_UDS_mul (const uintD* sourceptr1, uintC len1,
                        const uintD* sourceptr2, uintC len2,
                        uintD* destptr);
// Spezialfall sourceptr1 == sourceptr2 && len1 == len2.
extern void cl_UDS_mul_square (const uintD* sourceptr, uintC len,
                               uintD* destptr);

// Multipliziert zwei Unsigned-Digit-sequences.
// UDS_UDS_mul_UDS(len1,LSDptr1, len2,LSDptr2, MSDptr=,len=,LSDptr=);
// multipliziert die UDS ../len1/LSDptr1 und ../len2/LSDptr2.
// Dabei sollte len1>0 und len2>0 sein.
// Ergebnis ist die UDS MSDptr/len/LSDptr, mit len=len1+len2, im Stack.
// Dabei wird num_stack erniedrigt.
  #define UDS_UDS_mul_UDS(len1,LSDptr1,len2,LSDptr2, MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung)  \
    var uintC CONCAT(len_from_UDSmul_,__LINE__) = (uintC)(len1) + (uintC)(len2); \
    var uintD* CONCAT(LSDptr_from_UDSmul_,__LINE__);				\
    unused (len_zuweisung CONCAT(len_from_UDSmul_,__LINE__));			\
    num_stack_alloc(CONCAT(len_from_UDSmul_,__LINE__),MSDptr_zuweisung,LSDptr_zuweisung CONCAT(LSDptr_from_UDSmul_,__LINE__) =); \
    cl_UDS_mul((LSDptr1),(len1),(LSDptr2),(len2),CONCAT(LSDptr_from_UDSmul_,__LINE__));

// Multipliziert zwei Digit-sequences.
// DS_DS_mul_DS(MSDptr1,len1,LSDptr1, MSDptr2,len2,LSDptr2, MSDptr=,len=,LSDptr=);
// multipliziert die DS MSDptr1/len1/LSDptr1 und MSDptr2/len2/LSDptr2.
// Dabei sollte len1>0 und len2>0 sein, und beide DS sollten /= 0 sein.
// Alles sollten Variablen sein!
// Ergebnis ist die DS MSDptr/len/LSDptr, mit len=len1+len2, im Stack.
// Dabei wird num_stack erniedrigt.
  // Methode:
  // Erst unsigned multiplizieren. Dann bis zu zwei Subtraktionen.
  // Sei b=2^intDsize, k=len1, l=len2, n=DS1, m=DS2.
  // Gesucht ist n * m.
  // Wir errechnen erst das unsigned-product p (mod b^(k+l)).
  // n>0, m>0: p = n*m,             n*m = p
  // n<0, m>0: p = (n+b^k)*m,       n*m + b^(k+l) = p - b^k * m (mod b^(k+l)).
  // n>0, m<0: p = n*(m+b^l),       n*m + b^(k+l) = p - b^l * n (mod b^(k+l)).
  // n<0, m<0: p = (n+b^k)*(m+b^l),
  //           n*m = p - b^k * (m+b^l) - b^l * (n+b^k) (mod b^(k+l)).
  #define DS_DS_mul_DS(MSDptr1,len1,LSDptr1,MSDptr2,len2,LSDptr2, MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung)  \
    var uintD* MSDptr0;							\
    var uintD* LSDptr0;							\
    var uintC len_from_DSmal = (uintC)(len1) + (uintC)(len2);		\
    unused (len_zuweisung len_from_DSmal);				\
    num_stack_alloc(len_from_DSmal,MSDptr_zuweisung MSDptr0 =,LSDptr_zuweisung LSDptr0 =); \
    var uintD MSD1_from_DSmal = mspref(MSDptr1,0);			\
    var uintD MSD2_from_DSmal = mspref(MSDptr2,0);			\
    var uintC len1_from_DSmal = (len1);					\
    var uintC len2_from_DSmal = (len2);					\
    if (MSD1_from_DSmal==0) { msprefnext(MSDptr0) = 0; len1_from_DSmal--; } \
    if (MSD2_from_DSmal==0) { msprefnext(MSDptr0) = 0; len2_from_DSmal--; } \
    cl_UDS_mul((LSDptr1),len1_from_DSmal,(LSDptr2),len2_from_DSmal,LSDptr0); \
    if ((sintD)MSD1_from_DSmal < 0) /* n<0 ?                          */\
      /* muß m bzw. m+b^l subtrahieren, um k Digits verschoben:       */\
      { subfrom_loop_lsp(LSDptr2,LSDptr0 lspop len1,len2); }		\
    if ((sintD)MSD2_from_DSmal < 0) /* m<0 ?                          */\
      /* muß n bzw. n+b^k subtrahieren, um l Digits verschoben:       */\
      { subfrom_loop_lsp(LSDptr1,LSDptr0 lspop len2,len1); }


// Dividiert zwei Unsigned Digit sequences durcheinander.
// UDS_divide(a_MSDptr,a_len,a_LSDptr, b_MSDptr,b_len,b_LSDptr, &q,&r);
// Die UDS a = a_MSDptr/a_len/a_LSDptr (a>=0) wird durch
// die UDS b = b_MSDptr/b_len/b_LSDptr (b>=0) dividiert:
// a = q * b + r mit 0 <= r < b. Bei b=0 Error.
// q der Quotient, r der Rest.
// q = q_MSDptr/q_len/q_LSDptr, r = r_MSDptr/r_len/r_LSDptr beides
// Normalized Unsigned Digit sequences.
// Vorsicht: q_LSDptr <= r_MSDptr,
//           Vorzeichenerweiterung von r kann q zerstören!
//           Vorzeichenerweiterung von q ist erlaubt.
// a und b werden nicht modifiziert.
// num_stack wird erniedrigt.
  #define UDS_divide(a_MSDptr,a_len,a_LSDptr,b_MSDptr,b_len,b_LSDptr,q_,r_)  \
    /* Platz fürs Ergebnis machen. Brauche maximal a_len+1 Digits.    */\
    var uintC _a_len = (a_len);						\
    var uintD* roomptr; num_stack_alloc_1(_a_len+1,roomptr=,);		\
    cl_UDS_divide(a_MSDptr,_a_len,a_LSDptr,b_MSDptr,b_len,b_LSDptr,roomptr,q_,r_);
  extern void cl_UDS_divide (const uintD* a_MSDptr, uintC a_len, const uintD* a_LSDptr,
                             const uintD* b_MSDptr, uintC b_len, const uintD* b_LSDptr,
                             uintD* roomptr, DS* q_, DS* r_);


// Bildet zu einer Unsigned Digit sequence a die Wurzel
// (genauer: Gaußklammer aus Wurzel aus a).
// UDS_sqrt(a_MSDptr,a_len,a_LSDptr, &b, squarep=)
// > a_MSDptr/a_len/a_LSDptr: eine UDS
// < NUDS b: Gaußklammer der Wurzel aus a
// < squarep: true falls a = b^2, false falls b^2 < a < (b+1)^2.
// a wird nicht modifiziert.
// Vorzeichenerweiterung von b ist erlaubt.
// num_stack wird erniedrigt.
  #define UDS_sqrt(a_MSDptr,a_len,a_LSDptr,b_,squarep_zuweisung)  \
    { /* ceiling(a_len,2) Digits Platz fürs Ergebnis machen:          */\
      var uintC _a_len = (a_len);					\
      num_stack_alloc_1(ceiling(_a_len,2),(b_)->MSDptr=,);		\
      squarep_zuweisung cl_UDS_sqrt(a_MSDptr,_a_len,a_LSDptr,b_);	\
    }
  extern bool cl_UDS_sqrt (const uintD* a_MSDptr, uintC a_len, const uintD* a_LSDptr, DS* b_);


// Auxiliary function for approximately computing 1/x
// using Newton iteration.
  extern void cl_UDS_recip (const uintD* a_MSDptr, uintC a_len,
                            uintD* b_MSDptr, uintC b_len);

// Auxiliary function for approximately computing 1/sqrt(x)
// using Newton iteration.
  extern void cl_UDS_recipsqrt (const uintD* a_MSDptr, uintC a_len,
                                uintD* b_MSDptr, uintC b_len);

}  // namespace cln

#endif /* _CL_DS_H */
