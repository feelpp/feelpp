// cl_I internals

#ifndef _CL_I_H
#define _CL_I_H

#include "cln/number.h"
#include "cln/integer.h"
#include "base/cl_macros.h"
#include "cln/malloc.h"
#include "cln/exception.h"
#include "base/cl_offsetof.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

// maximal needed length of a digit sequence for a fixnum
  #define FN_maxlength  ceiling(cl_value_len,intDsize)
// maximal needed length (without sign) of a digit sequence for a fixnum
  #define pFN_maxlength  ceiling(cl_value_len-1,intDsize)
// minimum length of a digit sequence for a bignum
  #define bn_minlength  ceiling((cl_value_len+1),intDsize)
  // Because bignums with n < ceiling((cl_value_len+1)/intDsize) digits are
  // integers with at most n*intDsize < cl_value_len+1 bits, so they fit into
  // a fixnum, including the sign bit.
  // 1 <= bn_minlength <= 5.
// We have pFN_maxlength <= FN_maxlength <= bn_minlength.

// Private fixnum constructor.
inline cl_I::cl_I (struct cl_fixnum * null, cl_uint w)
	: cl_RA ((cl_private_thing) w) { unused null; }
inline const cl_I cl_I_from_word (cl_uint w)
{
	return cl_I((struct cl_fixnum *) 0, w);
}

inline cl_uint cl_FN_word (const cl_I& x)
{
	return x.word;
}


// Bignums.

struct cl_heap_bignum : cl_heap {
	uintC length;		// length (in digits)
	uintD data[1];		// number in two's complement representation
};

inline cl_heap_bignum* TheBignum (cl_heap_bignum* p)
	{ return p; }
inline cl_heap_bignum* TheBignum (const cl_number& obj)
	{ return (cl_heap_bignum*)(obj.pointer); }

inline cl_heap_bignum* allocate_bignum (uintC length)
{
	cl_heap_bignum* p = (cl_heap_bignum*) malloc_hook(offsetofa(cl_heap_bignum,data)+sizeof(uintD)*length);
	p->refcount = 1;
	p->type = &cl_class_bignum;
	p->length = length;
	return p;
}

// Private constructor.
// ptr should be the result of some allocate_bignum() call.
inline cl_I::cl_I (cl_heap_bignum* ptr)
	: cl_RA ((cl_private_thing) ptr) {}

// Both work, but the first definition results in less compiler-generated
// temporaries.
#if 1
  #define Bignum  cl_heap_bignum*
#else
  #define Bignum  cl_I
#endif


// Integers in general.

// Type tests.
inline bool integerp (const cl_I& x)
	{ unused x; return true; }
inline bool fixnump (const cl_I& x)
	{ return !x.pointer_p(); }
inline bool bignump (const cl_I& x)
	{ return x.pointer_p(); }

// Sign test:

// (MINUSP x) == (< x 0)
inline bool minusp (const cl_I& x)
{
	if (fixnump(x))
		// This assumes cl_value_shift + cl_value_len == cl_pointer_size.
		return (cl_sint) x.word < 0;
	else
		return (sintD)mspref(arrayMSDptr(TheBignum(x)->data,TheBignum(x)->length),0) < 0;
}

// (ZEROP x) == (= x 0)
inline bool zerop (const cl_I& x)
{
	return x.word == cl_combine(cl_FN_tag,0);
}

// (EQ x y) == (= x y), assuming y a fixnum
inline bool eq (const cl_I& x, sint32 y)
{
	return x.word == cl_combine(cl_FN_tag,y);
}


// Umwandlungsroutinen Integer <--> Longword:

// Wandelt Fixnum >=0 in Unsigned Longword um.
// FN_to_UV(obj)
// > obj: ein Fixnum >=0
// < ergebnis: der Wert des Fixnum als intVsize-Bit-Zahl.
inline uintV FN_to_UV (const cl_I& x)
{
	// This assumes cl_value_shift + cl_value_len == cl_pointer_size.
	return (cl_uint)(x.word) >> cl_value_shift;
}

// Wandelt Fixnum in Longword um.
// FN_to_V(obj)
// > obj: ein Fixnum
// < ergebnis: der Wert des Fixnum als intVsize-Bit-Zahl.
inline sintV FN_to_V (const cl_I& x)
{
	// This assumes cl_value_shift + cl_value_len == cl_pointer_size.
	return (cl_sint)(x.word) >> cl_value_shift;
}

// FN_V_zerop(x,x_) stellt fest, ob x = 0 ist.
// Dabei ist x ein Fixnum und x_ = FN_to_V(x).
  #define FN_V_zerop(x,x_)  (x_==0)

// FN_V_minusp(x,x_) stellt fest, ob x < 0 ist.
// Dabei ist x ein Fixnum und x_ = FN_to_V(x).
  #define FN_V_minusp(x,x_)  (x_<0)

#ifdef intQsize

// Wandelt Fixnum in Quadword um.
// FN_to_Q(obj)
// > obj: ein Fixnum
// < ergebnis: der Wert des Fixnum als 64-Bit-Zahl.
inline sint64 FN_to_Q (const cl_I& x)
{
	// This assumes cl_value_shift + cl_value_len == cl_pointer_size.
	return (cl_sint)(x.word) >> cl_value_shift;
}

// Wandelt Integer >=0 in Unsigned Quadword um.
// cl_I_to_UQ(obj)
// > obj: ein Objekt, sollte ein Integer >=0, <2^64 sein
// < ergebnis: der Wert des Integer als 64-Bit-Zahl.
  extern uint64 cl_I_to_UQ (const cl_I& obj);

// Wandelt Integer in Signed Quadword um.
// cl_I_to_Q(obj)
// > obj: ein Objekt, sollte ein Integer >=-2^63, <2^63 sein
// < ergebnis: der Wert des Integer als 64-Bit-Zahl.
  extern sint64 cl_I_to_Q (const cl_I& obj);

#endif

// Wandelt Longword in Fixnum um.
// L_to_FN(wert)
// > wert: Wert des Fixnums, ein signed 32-Bit-Integer
//         >= -2^(cl_value_len-1), < 2^(cl_value_len-1)
// < ergebnis: Fixnum mit diesem Wert.
  inline const cl_I L_to_FN (sint32 wert)
  {
	return cl_I_from_word(cl_combine(cl_FN_tag,wert));
  }

// Wandelt Longword in Integer um.
// L_to_I(wert)
// > wert: Wert des Integers, ein signed 32-Bit-Integer.
// < ergebnis: Integer mit diesem Wert.
#if (cl_value_len >= 32)
  #define L_to_I(wert)  L_to_FN(wert)
#else
  extern cl_private_thing cl_I_constructor_from_L (sint32 wert);
  inline const cl_I L_to_I (sint32 wert)
  {
	return cl_I(cl_I_constructor_from_L(wert));
  }
#endif

// Wandelt Unsigned Longword in Fixnum >=0 um.
// UL_to_FN(wert)
// > wert: Wert des Integers, ein unsigned 32-Bit-Integer < 2^(cl_value_len-1).
// < ergebnis: Fixnum mit diesem Wert.
  inline const cl_I UL_to_FN (uint32 wert)
  {
	return cl_I_from_word(cl_combine(cl_FN_tag,wert));
  }

// Wandelt Unsigned Longword in Integer >=0 um.
// UL_to_I(wert)
// > wert: Wert des Integers, ein unsigned 32-Bit-Integer.
// < ergebnis: Integer mit diesem Wert.
#if (cl_value_len > 32)
  #define UL_to_I(wert)  UL_to_FN(wert)
#else
  extern cl_private_thing cl_I_constructor_from_UL (uint32 wert);
  inline const cl_I UL_to_I (uint32 wert)
  {
	return cl_I(cl_I_constructor_from_UL(wert));
  }
#endif

#ifdef intQsize

// Wandelt Quadword in Integer um.
// Q_to_I(wert)
// > wert: Wert des Integers, ein signed 64-Bit-Integer.
// < ergebnis: Integer mit diesem Wert.
  extern cl_private_thing cl_I_constructor_from_Q (sint64 wert);
  inline const cl_I Q_to_I (sint64 wert)
  {
	return cl_I(cl_I_constructor_from_Q(wert));
  }

// Wandelt Unsigned Quadword in Integer >=0 um.
// UQ_to_I(wert)
// > wert: Wert des Integers, ein unsigned 64-Bit-Integer.
// < ergebnis: Integer mit diesem Wert.
  extern cl_private_thing cl_I_constructor_from_UQ (uint64 wert);
  inline const cl_I UQ_to_I (uint64 wert)
  {
	return cl_I(cl_I_constructor_from_UQ(wert));
  }

  extern cl_private_thing cl_I_constructor_from_Q2 (sint64 wert_hi, uint64 wert_lo );
  inline const cl_I Q2_to_I( sint64 wert_hi, uint64 wert_lo)
  {
	return cl_I(cl_I_constructor_from_Q2(wert_hi, wert_lo));
  }
#endif

// Wandelt Doppel-Longword in Integer um.
// L2_to_I(wert_hi,wert_lo)
// > wert_hi|wert_lo: Wert des Integers, ein signed 64-Bit-Integer.
// < ergebnis: Integer mit diesem Wert.
#if (cl_word_size==64)
  inline cl_private_thing cl_I_constructor_from_L2 (sint32 wert_hi, uint32 wert_lo)
  {
	return cl_I_constructor_from_Q(((sint64)wert_hi<<32) | (sint64)wert_lo);
  }
#else
  extern cl_private_thing cl_I_constructor_from_L2 (sint32 wert_hi, uint32 wert_lo);
#endif
  inline const cl_I L2_to_I (sint32 wert_hi, uint32 wert_lo)
  {
	return cl_I(cl_I_constructor_from_L2(wert_hi,wert_lo));
  }

// Wandelt Unsigned Doppel-Longword in Integer um.
// UL2_to_I(wert_hi,wert_lo)
// > wert_hi|wert_lo: Wert des Integers, ein unsigned 64-Bit-Integer.
// < ergebnis: Integer mit diesem Wert.
#if (cl_word_size==64)
  inline cl_private_thing cl_I_constructor_from_UL2 (uint32 wert_hi, uint32 wert_lo)
  {
	return cl_I_constructor_from_UQ(((uint64)wert_hi<<32) | (uint64)wert_lo);
  }
#else
  extern cl_private_thing cl_I_constructor_from_UL2 (uint32 wert_hi, uint32 wert_lo);
#endif
  inline const cl_I UL2_to_I (uint32 wert_hi, uint32 wert_lo)
  {
	return cl_I(cl_I_constructor_from_UL2(wert_hi,wert_lo));
  }

// Wandelt sintV in Integer um.
// V_to_I(wert)
// > wert: Wert des Integers, ein sintV.
// < ergebnis: Integer mit diesem Wert.
#if (intVsize<=32)
  #define V_to_I(wert)  L_to_I(wert)
#else
  #define V_to_I(wert)  Q_to_I(wert)
#endif

// Wandelt uintV in Integer >=0 um.
// UV_to_I(wert)
// > wert: Wert des Integers, ein uintV.
// < ergebnis: Integer mit diesem Wert.
#if (intVsize<=32)
  #define UV_to_I(wert)  UL_to_I(wert)
#else
  #define UV_to_I(wert)  UQ_to_I(wert)
#endif

// Wandelt sintE in Integer um.
// E_to_I(wert)
// > wert: Wert des Integers, ein sintE.
// < ergebnis: Integer mit diesem Wert.
#if (intEsize<=32)
  #define E_to_I(wert)  L_to_I(wert)
#else
  #define E_to_I(wert)  Q_to_I(wert)
#endif

// Wandelt uintE in Integer >=0 um.
// UE_to_I(wert)
// > wert: Wert des Integers, ein uintE.
// < ergebnis: Integer mit diesem Wert.
#if (intEsize<=32)
  #define UE_to_I(wert)  UL_to_I(wert)
#else
  #define UE_to_I(wert)  UQ_to_I(wert)
#endif

// Wandelt uintD in Integer >=0 um.
// UD_to_I(wert)
// > wert: Wert des Integers, ein uintD.
// < ergebnis: Integer mit diesem Wert.
#if (intDsize < cl_value_len)
  #define UD_to_I(wert)  UL_to_FN(wert)
#elif (intDsize<=32)
  #define UD_to_I(wert)  UL_to_I(wert)
#elif (intDsize==64)
  #define UD_to_I(wert)  UQ_to_I(wert)
#endif

// Liefert die Differenz x-y zweier Unsigned Longwords x,y als Integer.
// minus(x,y)
inline const cl_I minus (uintL x, uintL y)
{
#if (cl_word_size==64)
	return Q_to_I((sintQ)(uintQ)x - (sintQ)(uintQ)y);
#else
	return L2_to_I( (x<y ? -1 : 0), x-y );
#endif
}

#ifdef intQsize

inline const cl_I minus (uintQ x, uintQ y)
{
	return Q2_to_I( (x<y ? -1 : 0), x-y );
}

#endif

// Umwandlungsroutinen Digit sequence <--> Longword:

#if (intDsize<=32)

// Holt die nächsten pFN_maxlength Digits in ein uintV.
inline uintV pFN_maxlength_digits_at (const uintD* ptr)
{
#if (pFN_maxlength==1)
	return (uintV)lspref(ptr,0);
#elif (pFN_maxlength==2)
	return ((uintV)lspref(ptr,1)<<intDsize) | (uintV)lspref(ptr,0);
#elif (pFN_maxlength==3)
	return ((((uintV)lspref(ptr,2)<<intDsize) | (uintV)lspref(ptr,1))<<intDsize) | (uintV)lspref(ptr,0);
#elif (pFN_maxlength==4)
	return ((((((uintV)lspref(ptr,3)<<intDsize) | (uintV)lspref(ptr,2))<<intDsize) | (uintV)lspref(ptr,1))<<intDsize) | (uintV)lspref(ptr,0);
#elif (pFN_maxlength==5)
	return ((((((((uintV)lspref(ptr,4)<<intDsize) | (uintV)lspref(ptr,3))<<intDsize) | (uintV)lspref(ptr,2))<<intDsize) | (uintV)lspref(ptr,1))<<intDsize) | (uintV)lspref(ptr,0);
#elif (pFN_maxlength==6)
	return ((((((((((uintV)lspref(ptr,5)<<intDsize) | (uintV)lspref(ptr,4))<<intDsize) | (uintV)lspref(ptr,3))<<intDsize) | (uintV)lspref(ptr,2))<<intDsize) | (uintV)lspref(ptr,1))<<intDsize) | (uintV)lspref(ptr,0);
#elif (pFN_maxlength==7)
	return ((((((((((((uintV)lspref(ptr,6)<<intDsize) | (uintV)lspref(ptr,5))<<intDsize) | (uintV)lspref(ptr,4))<<intDsize) | (uintV)lspref(ptr,3))<<intDsize) | (uintV)lspref(ptr,2))<<intDsize) | (uintV)lspref(ptr,1))<<intDsize) | (uintV)lspref(ptr,0);
#elif (pFN_maxlength==8)
	return ((((((((((((((uintV)lspref(ptr,7)<<intDsize) | (uintV)lspref(ptr,6))<<intDsize) | (uintV)lspref(ptr,5))<<intDsize) | (uintV)lspref(ptr,4))<<intDsize) | (uintV)lspref(ptr,3))<<intDsize) | (uintV)lspref(ptr,2))<<intDsize) | (uintV)lspref(ptr,1))<<intDsize) | (uintV)lspref(ptr,0);
#endif
}

// Schreibt ein uintV in die nächsten pFN_maxlength Digits.
inline void set_pFN_maxlength_digits_at (uintD* ptr, uintV wert)
{
#if (pFN_maxlength==1)
	lspref(ptr,0) = (uintD)wert;
#elif (pFN_maxlength==2)
	lspref(ptr,1) = (uintD)(wert>>intDsize);
	lspref(ptr,0) = (uintD)(wert);
#elif (pFN_maxlength==3)
	lspref(ptr,2) = (uintD)(wert>>(2*intDsize));
	lspref(ptr,1) = (uintD)(wert>>intDsize);
	lspref(ptr,0) = (uintD)(wert);
#elif (pFN_maxlength==4)
	lspref(ptr,3) = (uintD)(wert>>(3*intDsize));
	lspref(ptr,2) = (uintD)(wert>>(2*intDsize));
	lspref(ptr,1) = (uintD)(wert>>intDsize);
	lspref(ptr,0) = (uintD)(wert);
#elif (pFN_maxlength==5)
	lspref(ptr,4) = (uintD)(wert>>(4*intDsize));
	lspref(ptr,3) = (uintD)(wert>>(3*intDsize));
	lspref(ptr,2) = (uintD)(wert>>(2*intDsize));
	lspref(ptr,1) = (uintD)(wert>>intDsize);
	lspref(ptr,0) = (uintD)(wert);
#elif (pFN_maxlength==6)
	lspref(ptr,5) = (uintD)(wert>>(5*intDsize));
	lspref(ptr,4) = (uintD)(wert>>(4*intDsize));
	lspref(ptr,3) = (uintD)(wert>>(3*intDsize));
	lspref(ptr,2) = (uintD)(wert>>(2*intDsize));
	lspref(ptr,1) = (uintD)(wert>>intDsize);
	lspref(ptr,0) = (uintD)(wert);
#elif (pFN_maxlength==7)
	lspref(ptr,6) = (uintD)(wert>>(6*intDsize));
	lspref(ptr,5) = (uintD)(wert>>(5*intDsize));
	lspref(ptr,4) = (uintD)(wert>>(4*intDsize));
	lspref(ptr,3) = (uintD)(wert>>(3*intDsize));
	lspref(ptr,2) = (uintD)(wert>>(2*intDsize));
	lspref(ptr,1) = (uintD)(wert>>intDsize);
	lspref(ptr,0) = (uintD)(wert);
#elif (pFN_maxlength==8)
	lspref(ptr,7) = (uintD)(wert>>(7*intDsize));
	lspref(ptr,6) = (uintD)(wert>>(6*intDsize));
	lspref(ptr,5) = (uintD)(wert>>(5*intDsize));
	lspref(ptr,4) = (uintD)(wert>>(4*intDsize));
	lspref(ptr,3) = (uintD)(wert>>(3*intDsize));
	lspref(ptr,2) = (uintD)(wert>>(2*intDsize));
	lspref(ptr,1) = (uintD)(wert>>intDsize);
	lspref(ptr,0) = (uintD)(wert);
#endif
}

#elif (intDsize==64)

// Holt die nächsten pFN_maxlength Digits in ein uint64.
inline uint64 pFN_maxlength_digits_at (const uintD* ptr)
{
	return (uint64)lspref(ptr,0);
}

#endif


// Umwandlungsroutinen Digit sequence --> Integer:

// Normalized Digit sequence to Integer
// NDS_to_I(MSDptr,len)
// Digit Sequence MSDptr/len/.. in Integer umwandeln.
  extern const cl_I NDS_to_I (const uintD* MSDptr, uintC len);

// Normalized Unsigned Digit Sequence to Integer
// NUDS_to_I(MSDptr,len)
// Normalized UDS MSDptr/len/.. in Integer >=0 umwandeln.
// Unterhalb von MSDptr muß 1 Digit Platz sein.
  extern const cl_I NUDS_to_I (uintD* MSDptr, uintC len);

// Unsigned Digit Sequence to Integer
// UDS_to_I(MSDptr,len)
// UDS MSDptr/len/.. in Integer >=0 umwandeln.
// Unterhalb von MSDptr muß 1 Digit Platz sein.
  extern const cl_I UDS_to_I (uintD* MSDptr, uintC len);

// Digit Sequence to Integer
// DS_to_I(MSDptr,len)
// DS MSDptr/len/.. in Integer umwandeln.
  extern const cl_I DS_to_I (const uintD* MSDptr, uintC len);


// Umwandlungsroutinen Integer --> Digit sequence:

// Unterteilung eines Fixnums in Digits:
// intDsize=8 -> MSD=LSD3,LSD2,LSD1,LSD0, sollte FN_maxlength=4 sein.
// intDsize=16 -> MSD=LSD1,LSD0, sollte FN_maxlength=2 sein.
// intDsize=32 -> MSD=LSD0, sollte FN_maxlength=1 sein.

inline sintD FN_MSD (cl_uint word)
{
	// This assumes cl_value_shift + cl_value_len == cl_pointer_size.
	return (cl_sint)word >> (cl_value_shift + (FN_maxlength-1)*intDsize);
}

#if (FN_maxlength==1)
  #define FN_LSD0(word)  FN_MSD(word)
  #define FN_LSD1(word)  (throw runtime_exception(), (uintD)0)  // never used
  #define FN_LSD2(word)  (throw runtime_exception(), (uintD)0)  // never used
  #define FN_LSD3(word)  (throw runtime_exception(), (uintD)0)  // never used
#endif
#if (FN_maxlength==2)
  inline uintD FN_LSD0 (cl_uint word)
  {
	return (uintD)(word >> cl_value_shift);
  }
  #define FN_LSD1(word)  FN_MSD(word)
  #define FN_LSD2(word)  (throw runtime_exception(), (uintD)0)  // never used
  #define FN_LSD3(word)  (throw runtime_exception(), (uintD)0)  // never used
#endif
#if (FN_maxlength==4)
  inline uintD FN_LSD0 (cl_uint word)
  {
	return (uintD)(word >> cl_value_shift);
  }
  inline uintD FN_LSD1 (cl_uint word)
  {
	return (uintD)(word >> (cl_value_shift + intDsize));
  }
  inline uintD FN_LSD2 (cl_uint word)
  {
	return (uintD)(word >> (cl_value_shift + 2*intDsize));
  }
  #define FN_LSD3(word)  FN_MSD(word)
#endif

// wird nur bei FN_maxlength >= 2 gebraucht, d.h. intDsize < cl_value_len
#define FN_MSD1_mask  ((~(cl_uint)(bitc(intDsize-1)-1)) << cl_value_shift)
// wird nur bei FN_maxlength >= 3 gebraucht, d.h. 2*intDsize < cl_value_len
#define FN_MSD2_mask  ((~(cl_uint)(bitc(2*intDsize-1)-1)) << cl_value_shift)
// wird nur bei FN_maxlength >= 4 gebraucht, d.h. 3*intDsize < cl_value_len
#define FN_MSD3_mask  ((~(cl_uint)(bitc(3*intDsize-1)-1)) << cl_value_shift)

// Store a Fixnum at destLSDptr, <= FN_maxlength digits below destLSDptr needed.
#define FN_to_NDS(destLSDptr, word, MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung, check_for_0,zero_statement)  \
{ var uintC len_from_FN_to_NDS;						\
  if (check_for_0 && (word == cl_combine(cl_FN_tag,0)))			\
    { len_from_FN_to_NDS = 0; zero_statement }				\
    else								\
    { var cl_uint testMSD;						\
      if ((FN_maxlength<=1) ||						\
          ((testMSD = (word) & FN_MSD1_mask) == 0) || (testMSD == FN_MSD1_mask) \
         )								\
        { len_from_FN_to_NDS = 1;					\
          lspref(destLSDptr,0) = FN_LSD0(word);				\
        }								\
      elif ((FN_maxlength<=2) ||					\
            ((testMSD = (word) & FN_MSD2_mask) == 0) || (testMSD == FN_MSD2_mask) \
           )								\
        { len_from_FN_to_NDS = 2;					\
          lspref(destLSDptr,0) = FN_LSD0(word);				\
          lspref(destLSDptr,1) = FN_LSD1(word);				\
        }								\
      elif ((FN_maxlength<=3) ||					\
            ((testMSD = (word) & FN_MSD3_mask) == 0) || (testMSD == FN_MSD3_mask) \
           )								\
        { len_from_FN_to_NDS = 3;					\
          lspref(destLSDptr,0) = FN_LSD0(word);				\
          lspref(destLSDptr,1) = FN_LSD1(word);				\
          lspref(destLSDptr,2) = FN_LSD2(word);				\
        }								\
      else /* (FN_maxlength<=4) */					\
        { len_from_FN_to_NDS = 4;					\
          lspref(destLSDptr,0) = FN_LSD0(word);				\
          lspref(destLSDptr,1) = FN_LSD1(word);				\
          lspref(destLSDptr,2) = FN_LSD2(word);				\
          lspref(destLSDptr,3) = FN_LSD3(word);				\
        }								\
    }									\
  unused (MSDptr_zuweisung (destLSDptr) lspop len_from_FN_to_NDS);	\
  unused (len_zuweisung len_from_FN_to_NDS);				\
  unused (LSDptr_zuweisung (destLSDptr));				\
}

// Bignum to Normalized Digit sequence, Kopieren unnötig
// BN_to_NDS_nocopy(obj, MSDptr=,len=,LSDptr=);
// > obj: ein Bignum
// < MSDptr/len/LSDptr: Normalized Digit sequence
#if CL_DS_BIG_ENDIAN_P
  #define BN_to_NDS_nocopy(obj, MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung)  \
    { var Bignum bn_from_BN_to_NDS_nocopy = TheBignum(obj);		\
      unused (MSDptr_zuweisung (const uintD*) &bn_from_BN_to_NDS_nocopy->data[0]);	\
      unused (LSDptr_zuweisung (const uintD*) &bn_from_BN_to_NDS_nocopy->data[(uintP)( len_zuweisung bn_from_BN_to_NDS_nocopy->length )]);		\
    }
#else
  #define BN_to_NDS_nocopy(obj, MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung)  \
    { var Bignum bn_from_BN_to_NDS_nocopy = TheBignum(obj);		\
      unused (LSDptr_zuweisung (const uintD*) &bn_from_BN_to_NDS_nocopy->data[0]);	\
      unused (MSDptr_zuweisung (const uintD*) &bn_from_BN_to_NDS_nocopy->data[(uintP)( len_zuweisung bn_from_BN_to_NDS_nocopy->length )]);		\
    }
#endif
  inline const uintD* BN_MSDptr (const cl_I& obj)
    { var Bignum bn = TheBignum(obj); return (const uintD*) arrayMSDptr(bn->data,bn->length); }
  inline const uintD* BN_LSDptr (const cl_I& obj)
    { var Bignum bn = TheBignum(obj); return (const uintD*) arrayLSDptr(bn->data,bn->length); }

// Bignum to Normalized Digit sequence
// BN_to_NDS(obj, MSDptr=,len=,LSDptr=);
// > obj: ein Bignum
// < MSDptr/len/LSDptr: Normalized Digit sequence, darf modifiziert werden.
// Dabei wird num_stack erniedrigt.
  #define BN_to_NDS(obj, MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung)  \
    { var const cl_I& obj_from_BN_to_NDS = (obj);			\
      var uintD* MSDptr_from_BN_to_NDS;					\
      var uintC len_from_BN_to_NDS;					\
      len_zuweisung len_from_BN_to_NDS = TheBignum(obj_from_BN_to_NDS)->length; \
      num_stack_alloc(len_from_BN_to_NDS, MSDptr_zuweisung MSDptr_from_BN_to_NDS = , LSDptr_zuweisung); \
      copy_loop_msp(BN_MSDptr(obj_from_BN_to_NDS),MSDptr_from_BN_to_NDS,len_from_BN_to_NDS); \
    }

// Bignum to Normalized Digit sequence
// BN_to_NDS_1(obj, MSDptr=,len=,LSDptr=);
// > obj: ein Bignum
// < MSDptr/len/LSDptr: Normalized Digit sequence, darf modifiziert werden.
// Unterhalb von MSDptr ist noch 1 Digit Platz.
// Dabei wird num_stack erniedrigt.
  #define BN_to_NDS_1(obj, MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung)  \
    { var const cl_I& obj_from_BN_to_NDS = (obj);			\
      var uintD* MSDptr_from_BN_to_NDS;					\
      var uintC len_from_BN_to_NDS;					\
      len_zuweisung len_from_BN_to_NDS = TheBignum(obj_from_BN_to_NDS)->length; \
      num_stack_alloc_1(len_from_BN_to_NDS, MSDptr_zuweisung MSDptr_from_BN_to_NDS = , LSDptr_zuweisung); \
      copy_loop_msp(BN_MSDptr(obj_from_BN_to_NDS),MSDptr_from_BN_to_NDS,len_from_BN_to_NDS); \
    }

// Integer to Normalized Digit sequence, Kopieren unnötig.
// I_to_NDS_nocopy(obj, MSDptr=,len=,LSDptr=,check_for_0,zero_statement);
// > obj: ein Integer
// > check_for_0: ob obj möglicherweise =0 sein kann
// > zero_statement: wird bei obj=0 ausgeführt
// < MSDptr/len/LSDptr: Normalized Digit sequence
  #define I_to_NDS_nocopy(obj, MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung, check_for_0,zero_statement)  \
    var uintD CONCAT(FN_store_,__LINE__) [FN_maxlength];		\
    { var const cl_I& obj_from_I_to_NDS_nocopy = (obj);			\
      if (fixnump(obj_from_I_to_NDS_nocopy))				\
        { FN_to_NDS(arrayLSDptr(CONCAT(FN_store_,__LINE__),FN_maxlength), cl_FN_word(obj_from_I_to_NDS_nocopy), MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung, check_for_0,zero_statement); } \
        else								\
        { BN_to_NDS_nocopy(obj_from_I_to_NDS_nocopy,MSDptr_zuweisung,len_zuweisung, LSDptr_zuweisung); } \
    }

// Integer to Normalized Digit sequence
// I_to_NDS(obj, MSDptr=,len=,LSDptr=);
// > obj: ein Integer
// < MSDptr/len/LSDptr: Normalized Digit sequence, darf modifiziert werden.
// Dabei wird num_stack erniedrigt.
  #define I_to_NDS(obj, MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung)  \
    var uintD CONCAT(FN_store_,__LINE__) [FN_maxlength];		\
    { var const cl_I& obj_from_I_to_NDS = (obj);			\
      if (fixnump(obj_from_I_to_NDS))					\
        { FN_to_NDS(arrayLSDptr(CONCAT(FN_store_,__LINE__),FN_maxlength), cl_FN_word(obj_from_I_to_NDS), MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung, true,); } \
        else                                                              \
        { BN_to_NDS(obj_from_I_to_NDS,MSDptr_zuweisung,len_zuweisung, LSDptr_zuweisung); } \
    }

// Integer to Normalized Digit sequence
// I_to_NDS_1(obj, MSDptr=,len=,LSDptr=);
// > obj: ein Integer
// < MSDptr/len/LSDptr: Normalized Digit sequence, darf modifiziert werden.
// Unterhalb von MSDptr ist noch 1 Digit Platz.
// Dabei wird num_stack erniedrigt.
  #define I_to_NDS_1(obj, MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung)  \
    var uintD CONCAT(FN_store_,__LINE__) [1+FN_maxlength];		\
    { var const cl_I& obj_from_I_to_NDS = (obj);			\
      if (fixnump(obj_from_I_to_NDS))					\
        { FN_to_NDS(arrayLSDptr(CONCAT(FN_store_,__LINE__),1+FN_maxlength), cl_FN_word(obj_from_I_to_NDS), MSDptr_zuweisung,len_zuweisung,LSDptr_zuweisung, true,); } \
        else                                                              \
        { BN_to_NDS_1(obj_from_I_to_NDS,MSDptr_zuweisung,len_zuweisung, LSDptr_zuweisung); } \
    }


// Division ganzer Zahlen

// Dividiert zwei Integers x,y >=0 und liefert Quotient und Rest
// der Division x/y. Bei y=0 Error.
// cl_divide(x,y)
// > x,y: Integers >=0
// < q,r: Quotient q, Rest r
  extern const cl_I_div_t cl_divide (const cl_I& x, const cl_I& y);


// ggT und kgV von Integers

  // Teilfunktion für die Durchführung des Euklid-Algorithmus auf
  // den führenden Ziffern a' und b':
  // partial_gcd(a',b',&erg); mit a'>b'
  // liefert in erg: x1,y1,x2,y2 mit den in cl_I_gcd.cc angegebenen Invarianten.
  typedef struct { uintD x1,y1,x2,y2; } partial_gcd_result;
  extern void partial_gcd (uintD z1, uintD z2, partial_gcd_result* erg);
  #if HAVE_DD
    extern void partial_gcd (uintDD z1, uintDD z2, partial_gcd_result* erg);
  #else
    extern void partial_gcd (uintD z1hi, uintD z1lo, uintD z2hi, uintD z2lo, partial_gcd_result* erg);
  #endif


// Wurzel aus ganzen Zahlen

// Stellt fest, ob ein Integer >=0 eine n-te Potenz ist.
// cl_rootp_aux(x,n,&w)
// > x: ein Integer >=0
// > n: ein Integer >0
// > Annahme: x > 1 und n < (integer-length x).
// < w: Integer (expt x (/ n)) falls x eine n-te Potenz
// < ergebnis: true            ........................, false sonst
  extern bool cl_rootp_aux (cl_I x, uintL n, cl_I* w);


// Hilfsfunktion zur Eingabe von Integers

// Wandelt eine Ziffernfolge in ein Integer >=0 um.
// digits_to_I(MSBptr,len,base)
// > base: Stellenwertsystem-Basis, >=2, <=36
// > MSBptr/len/..: Ziffernfolge, bestehend aus Punkten (werden überlesen)
//     und Ziffern/Buchstaben mit Wert < base.
// < ergebnis: der dargestellte Integer >=0
  extern const cl_I digits_to_I (const char * MSBptr, uintC len, uintD base);


// Hilfsfunktion zur Ausgabe von Integers

// cl_digits_need(len,base) liefert eine obere Abschätzung für die Anzahl der
// Ziffern im Stellenwertsystem der Basis base, die x >= 0 braucht.
  extern uintC cl_digits_need (const cl_I& x, uintL base);

// Wandelt ein Integer in ein Stellensystem um.
// I_to_digits(x,base, &ergebnis);
// > x: ein Integer >=0
// > base: Stellensystem-Basis, 2 <= base <= 36.
// > ergebnis.LSBptr: darunter ist mindestens digits_need(len) Bytes Platz
// < ergebnis: fertige Folge MSBptr/len/LSBptr von Ziffern
  typedef struct { uintB* MSBptr; uintC len; uintB* LSBptr; } cl_digits;
  extern void I_to_digits (const cl_I& x, uintD base, cl_digits* erg);


// Hash code.
  extern unsigned long hashcode (const cl_I& x);


// A fixnum (cl_FN) is an immediate integer.

// typedef
class cl_FN : public cl_I {
public:
// Optimization of method pointer_p().
	bool pointer_p() const
		{ return false; }
};

inline bool fixnump (const cl_FN& x)
	{ unused x; return true; }
inline bool bignump (const cl_FN& x)
	{ unused x; return false; }

inline bool minusp (const cl_FN& x)
{
	// This assumes cl_value_shift + cl_value_len == cl_pointer_size.
	return (cl_sint) x.word < 0;
}


// A bignum (cl_BN) is a heap-allocated integer.

// typedef
class cl_BN : public cl_I {
public:
// Optimization of method pointer_p().
	bool pointer_p() const
		{ return true; }
};

inline bool fixnump (const cl_BN& x)
	{ unused x; return false; }
inline bool bignump (const cl_BN& x)
	{ unused x; return true; }

inline bool minusp (const cl_BN& x)
{
	return (sintD)mspref(arrayMSDptr(TheBignum(x)->data,TheBignum(x)->length),0) < 0;
}
inline bool zerop (const cl_BN& x)
	{ unused x; return false; }

}  // namespace cln

#endif /* _CL_I_H */
