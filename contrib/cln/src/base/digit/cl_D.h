// Digit level arithmetic

#ifndef _CL_D_H
#define _CL_D_H

#include "cln/types.h"
#include "base/cl_low.h"

// Aus cln/types.h importiere:
// intDsize        Anzahl Bits in einem Digit
// uintD, sintD    Integer-Typen für ein Digit
// log2_intDsize   log2(intDsize)
// HAVE_DD         Flag, das anzeigt, ob ein Integertyp für Doppel-Digits da ist
// intDDsize       Anzahl Bits in einem Doppel-Digit
// uintDD,sintDD   Integer-Typen für ein Doppel-Digit

#ifdef HAVE_FAST_LONGLONG
  #if !((64%intDsize)==0)
    #error "intDsize should be a divisor of 64!"
  #endif
#else
  #if !((32%intDsize)==0)
    #error "intDsize should be a divisor of 32!"
  #endif
#endif

namespace cln {

// Vorzeichen eines Digit bestimmen
// sign_of_sintD(wert)
// > wert: ein Digit
// < sintD ergebnis: 0 falls wert>=0, -1 falls wert<0.
inline sint32 sign_of_sintD (sintD wert)
{
	return sign_of(wert);
}

#if HAVE_DD

// High-Digit eines Doppel-Digit bestimmen
// highD(wert)
  #if (!(intDsize==16))
    #define highD(x)  ((uintD)((uintDD)(x)>>intDsize))
  #else
    #define highD  high16
  #endif

// Low-Digit eines Doppel-Digit bestimmen
// lowD(wert)
  #define lowD(x)  ((uintD)(uintDD)(x))

// Ein Doppel-Digit aus ihrem High-Digit und ihrem Low-Digit bestimmen:
// highlowDD(uintD high, uintD low)
  #if (!(intDsize==16))
    #define highlowDD(x,y)  (((uintDD)(uintD)(x)<<intDsize)|(uintDD)(uintD)(y))
  #else
    #define highlowDD  highlow32
  #endif

// Ein Doppel-Digit aus ihrem High-Digit und ihrem Low-Digit 0 bestimmen:
// highlowDD_0(uintD high)
  #if (!(intDsize==16))
    #define highlowDD_0(x)  ((uintDD)(uintD)(x)<<intDsize)
  #else
    #define highlowDD_0  highlow32_0
  #endif

#endif

// Zwei Digits multiplizieren:
// (uintDD)hilo = muluD(uintD arg1, uintD arg2)
// bzw.
// muluD(uintD arg1, uintD arg2, uintD hi =, uintD lo =);
#if HAVE_DD
  #if (intDsize==8)
    #ifdef __GNUC__
      #define muluD(arg1,arg2)  ((uintDD)((uintD)(arg1)*(uintD)(arg2)))
    #else
      #define muluD(arg1,arg2)  ((uintDD)(uintD)(arg1)*(uintDD)(uintD)(arg2))
    #endif
  #endif
  #if (intDsize==16)
    #define muluD  mulu16
  #endif
  #if (intDsize==32) && defined(HAVE_LONGLONG)
    #define muluD(arg1,arg2)  ((uintDD)(uintD)(arg1)*(uintDD)(uintD)(arg2))
  #endif
#else
  #if (intDsize==32)
    #define muluD  mulu32
  #endif
  #if (intDsize==64)
    #define muluD  mulu64
  #endif
#endif

// Zwei Digits multiplizieren, mit einem Digit als Ergebnis.
// (uintD)lo = muluD_unchecked(uintD arg1, uintD arg2)
// Es wird vorausgesetzt, daß arg1*arg2 < 2^intDsize.
  #if (intDsize==8) || (intDsize==16) || (intDsize==64)
    #define muluD_unchecked(arg1,arg2)  ((uintD)((uintD)(arg1)*(uintD)(arg2)))
  #endif
  #if (intDsize==32)
    #define muluD_unchecked(arg1,arg2)  mulu32_unchecked(arg1,arg2)
  #endif

// Durch ein Digit dividieren:
// divuD(uintDD x, uintD y, uintD q =, uintD r =);
// bzw.
// divuD(uintD xhi, uintD xlo, uintD y, uintD q =, uintD r =);
// dividiert x/y und liefert q = floor(x/y) und r = (x mod y). x = q*y+r.
// Es wird vorausgesetzt, daß 0 <= x < 2^intDsize*y.
#if HAVE_DD
  #if (intDsize==8)
    #define divuD  divu_1616_1616
  #endif
  #if (intDsize==16)
    #define divuD  divu_3216_1616
  #endif
  #if (intDsize==32) && defined(HAVE_LONGLONG)
    #define divuD(x,y,q_zuweisung,r_zuweisung) \
      { var uint64 __x = (x);                                 \
        var uint32 __y = (y);                                 \
        var uint32 __q = floor(__x,(uint64)__y);              \
        q_zuweisung __q; r_zuweisung (uint32)__x - __q * __y; \
      }
  #endif
#else
  #if (intDsize==32)
    #define divuD  divu_6432_3232
  #endif
  #if (intDsize==64)
    #define divuD  divu_12864_6464
  #endif
#endif

// Durch ein Digit dividieren:
// floorD(uintD x, uintD y)
// dividiert x/y und liefert q = floor(x/y).
// Es wird vorausgesetzt, daß y > 0.
  #if (intDsize==8) || (intDsize==16) || (intDsize==64)
    #define floorD(arg1,arg2)  (floor((uintD)(arg1),(uintD)(arg2)))
  #endif
  #if (intDsize==32)
    #define floorD  divu_3232_3232_
  #endif

// Ganzzahl-Wurzel eines Doppel-Digits berechnen.
// isqrtD(xhi,xlo,y=,sqrtp=);
// > uintD xhi,xlo: Radikand x = 2^intDsize*xhi+xlo,
//                  >= 2^(2*intDsize-2), < 2^(2*intDsize)
// < uintD y: floor(sqrt(x)), >= 2^(intDsize-1), < 2^intDsize
// < boolean sqrtp: /=0, falls x=y^2
#if (intDsize==8)
  #define isqrtD(xhi,xlo,y_zuweisung,sqrtp_zuweisung)  \
    { var uint32 _z;								\
      isqrt_32_16((((uint32)xhi<<8) | (uint32)xlo) << 16, _z=,sqrtp_zuweisung);	\
      y_zuweisung (_z >> 8);							\
    }
#endif
#if (intDsize==16)
  #define isqrtD(xhi,xlo,y_zuweisung,sqrtp_zuweisung)  \
    isqrt_32_16(highlow32(xhi,xlo),y_zuweisung,sqrtp_zuweisung)
#endif
#if (intDsize==32)
  #define isqrtD  isqrt_64_32
#endif
#if (intDsize==64)
  #define isqrtD  isqrt_128_64
#endif

// Bits eines Digit zählen:
// integerlengthD(digit,size=);
// setzt size auf die höchste in digit vorkommende Bitnummer.
// > digit: ein uintD >0
// < size: >0, <=intDsize, mit 2^(size-1) <= digit < 2^size
#if (intDsize==8)
  #define integerlengthD  integerlength8
#endif
#if (intDsize==16)
  #define integerlengthD  integerlength16
#endif
#if (intDsize==32)
  #define integerlengthD  integerlength32
#endif
#if (intDsize==64)
  #define integerlengthD  integerlength64
#endif

// Hintere Nullbits eines Digits zählen:
// ord2_D(digit,count=);
// setzt size auf die kleinste in digit vorkommende Bitnummer.
// > digit: ein uintD >0
// < count: >=0, <intDsize, mit 2^count | digit, digit/2^count ungerade
  #if defined(FAST_ORD2)
    #define ord2_D(digit,count_zuweisung)  \
      ord2_32((uint32)(digit),count_zuweisung)
  #else
    // Sei n = ord2(x). Dann ist logxor(x,x-1) = 2^n + (2^n-1) = 2^(n+1)-1.
    // Also  (ord2 x) = (1- (integer-length (logxor x (1- x)))) .
    #define ord2_D(digit,count_zuweisung)  \
      { var uintD _digit = digit ^ (digit - 1);		\
        integerlengthD(_digit,count_zuweisung -1 + )	\
      }
  #endif

// Bits eines Wortes zählen.
// logcountD(x)
// > x: ein uintD
// < ergebnis: Anzahl der darin gesetzten Bits
#if (intDsize==8)
  inline uint8 logcountD (uint8 x8) { logcount_8(); return x8; }
#endif
#if (intDsize==16)
  inline uint16 logcountD (uint16 x16) { logcount_16(); return x16; }
#endif
#if (intDsize==32)
  inline uint32 logcountD (uint32 x32) { logcount_32(); return x32; }
#endif
#if (intDsize==64)
  inline uint64 logcountD (uint64 x64) { logcount_64(); return x64; }
#endif

}  // namespace cln

#endif /* _CL_D_H */
