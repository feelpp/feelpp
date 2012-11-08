// Low-level arithmetic: operations on 16-bit and 32-bit words

#ifndef _CL_LOW_H
#define _CL_LOW_H

namespace cln {

// Determines the sign of a 16-bit number.
// sign_of(wert)
// > wert: eine 16-Bit-Zahl
// < sint16 ergebnis: 0 falls wert>=0, -1 falls wert<0.
inline sint16 sign_of (sint16 wert)
{
#if defined(__sparc64__)
	return (sint64)wert >> 63;
#elif defined(__sparc__) || defined(__arm__)
	return (sint32)wert >> 31;
#else
	return (wert >= 0 ? 0 : -1);
#endif
}

// Determines the sign of a 32-bit number.
// sign_of(wert)
// > wert: eine 32-Bit-Zahl
// < sint32 ergebnis: 0 falls wert>=0, -1 falls wert<0.
inline sint32 sign_of (sint32 wert)
{
#if defined(__sparc64__)
	return (sint64)wert >> 63;
#elif defined(__sparc__) || defined(__arm__)
	return wert >> 31;
#else
	return (wert >= 0 ? 0 : -1);
#endif
}

#ifdef HAVE_FAST_LONGLONG

// Determines the sign of a 64-bit number.
// sign_of(wert)
// > wert: eine 64-Bit-Zahl
// < sint64 ergebnis: 0 falls wert>=0, -1 falls wert<0.
inline sint64 sign_of (sint64 wert)
{
	return wert >> 63;
}

#endif /* HAVE_FAST_LONGLONG */


// High-Word einer 32-Bit-Zahl bestimmen
// high16(wert)
inline uint16 high16 (uint32 wert)
{
	return wert >> 16;
}

// Low-Word einer 32-Bit-Zahl bestimmen
// low16(wert)
inline uint16 low16 (uint32 wert)
{
	return (uint16)wert;
}

// Eine 32-Bit-Zahl aus ihrem High-Word und ihrem Low-Word bestimmen:
// highlow32(uint16 high, uint16 low)
inline uint32 highlow32 (uint16 high, uint16 low)
{
	return ((uint32)high << 16) | (uint32)low;
}

// Eine 32-Bit-Zahl aus ihrem High-Word und ihrem Low-Word 0 bestimmen:
// highlow32_0(uint16 high)
inline uint32 highlow32_0 (uint16 high)
{
	return (uint32)high << 16;
}

#ifdef HAVE_LONGLONG

// High-Word einer 64-Bit-Zahl bestimmen
// high32(wert)
inline uint32 high32 (uint64 wert)
{
	return wert >> 32;
}

// Low-Word einer 64-Bit-Zahl bestimmen
// low32(wert)
inline uint32 low32 (uint64 wert)
{
	return (uint32)wert;
}

// Eine 64-Bit-Zahl aus ihrem High-Word und ihrem Low-Word bestimmen:
// highlow64(uint32 high, uint32 low)
inline uint64 highlow64 (uint32 high, uint32 low)
{
	return ((uint64)high << 32) | (uint64)low;
}

// Eine 64-Bit-Zahl aus ihrem High-Word und ihrem Low-Word 0 bestimmen:
// highlow64_0(uint32 high)
inline uint64 highlow64_0 (uint32 high)
{
	return (uint64)high << 32;
}

#endif /* HAVE_LONGLONG */


// Multipliziert zwei 16-Bit-Zahlen miteinander und liefert eine 32-Bit-Zahl:
// mulu16(arg1,arg2)
// > arg1, arg2 : zwei 16-Bit-Zahlen
// < ergebnis: eine 32-Bit-Zahl
#if defined(__GNUC__) && defined(__sparc__) && !defined(__sparc64__) && defined(FAST_DOUBLE)
// Ist das schneller als mulu16_ ??
inline uint32 mulu16 (uint16 arg1, uint16 arg2)
{
	union { double f; uint32 i[2]; } __fi;
	__fi.f = (double)(sint32)arg1 * (double)(sint32)arg2
		 + (double)(4503599627370496.0L); // + 2^52, zum Normalisieren
	return __fi.i[1]; // untere 32 Bit herausholen (benutzt CL_CPU_BIG_ENDIAN_P !)
}
#elif defined(__GNUC__) && defined(__sparc64__) && !defined(NO_ASM)
inline uint32 mulu16 (uint16 arg1, uint16 arg2)
{
	register uint64 _prod;
	__asm__("umul %1,%2,%0"
		: "=r" (_prod)
		: "r" (arg1), "r" (arg2)
	       );
	return _prod;
}
#elif defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__)) && !defined(NO_ASM)
inline uint32 mulu16 (uint16 arg1, uint16 arg2)
{
	register uint16 _hi;
	register uint16 _lo;
	__asm__("mulw %2"
		: "=d" /* %dx */ (_hi), "=a" /* %ax */ (_lo)
		: "rm" (arg1), "1" /* %eax */ (arg2)
	       );
	return highlow32(_hi,_lo);
}
#elif (defined(__sparc__) || defined(__sparc64__)) && !defined(NO_ASM)
  extern "C" uint32 mulu16_ (uint16 arg1, uint16 arg2);
  #define mulu16  mulu16_  // extern in Assembler
#else
inline uint32 mulu16 (uint16 arg1, uint16 arg2)
{
	return arg1 * arg2;
}
#endif

// Multipliziert zwei 24-Bit-Zahlen zusammen und liefert eine 48-Bit-Zahl.
// mulu24(arg1,arg2,hi=,lo=);
// > arg1, arg2 : zwei 24-Bit-Zahlen
// < 2^32*hi+lo : eine 48-Bit-Zahl
#if defined(__sparc__) && !defined(__sparc64__) && defined(FAST_DOUBLE)
  #define mulu24(x,y,hi_zuweisung,lo_zuweisung)  \
    { var uint32 _x = (x);					\
      var uint32 _y = (y);					\
      var union { double f; uint32 i[2]; uint16 s[4]; } __fi;	\
      __fi.f = (double)(sint32)(_x)*(double)(sint32)(_y)	\
               + (double)(4503599627370496.0L); /* + 2^52, zum Normalisieren */\
      unused (hi_zuweisung __fi.s[1]); /* mittlere 16 Bit herausholen, (benutzt CL_CPU_BIG_ENDIAN_P !) */\
      lo_zuweisung __fi.i[1]; /* untere 32 Bit herausholen (benutzt CL_CPU_BIG_ENDIAN_P !)    */\
    }
#else
  #define mulu24  mulu32
#endif

// Multipliziert zwei 32-Bit-Zahlen miteinander und liefert eine 32-Bit-Zahl:
// mulu32_unchecked(arg1,arg2)
// > arg1, arg2 : zwei 32-Bit-Zahlen
// < ergebnis : eine 32-Bit-Zahl
// Es wird vorausgesetzt, daß arg1*arg2 < 2^32.
#if defined(__GNUC__) && defined(__sparc64__) && !defined(NO_ASM)
inline uint32 mulu32_unchecked (uint32 arg1, uint32 arg2)
{
	register uint64 _prod;
	__asm__("umul %1,%2,%0"
		: "=r" (_prod)
		: "r" (arg1), "r" (arg2)
	       );
	return _prod;
}
#elif defined(__sparc__) && !defined(NO_ASM)
  extern "C" uint32 mulu32_unchecked (uint32 x, uint32 y); // extern in Assembler
#else
  // Wir können dafür auch die Bibliotheksroutine des C-Compilers nehmen:
  inline uint32 mulu32_unchecked (uint32 arg1, uint32 arg2)
  {
	return arg1 * arg2;
  }
#endif

// Multipliziert zwei 32-Bit-Zahlen miteinander und liefert eine 64-Bit-Zahl:
// mulu32(arg1,arg2,hi=,lo=);
// > arg1, arg2 : zwei 32-Bit-Zahlen
// < 2^32*hi+lo : eine 64-Bit-Zahl
  extern "C" uint32 mulu32_ (uint32 arg1, uint32 arg2); // -> Low-Teil
#ifdef _MSC_VER
  // Workaround MSVC compiler bug: extern "C" results in wrong symbols, when
  // declared inside a namespace!
} extern "C" uint32 mulu32_high; namespace cln {        // -> High-Teil
#else
  extern "C" uint32 mulu32_high;                        // -> High-Teil
#endif
#if defined(__GNUC__) && defined(__m68k__) && !defined(NO_ASM)
  #define mulu32(x,y,hi_zuweisung,lo_zuweisung)  \
    ({ var uint32 _x = (x);       \
       var uint32 _y = (y);       \
       var uint32 _hi;            \
       var uint32 _lo;            \
       __asm__("mulul %3,%0:%1" : "=d" (_hi), "=d"(_lo) : "1" (_x), "dm" (_y) ); \
       unused (hi_zuweisung _hi); \
       lo_zuweisung _lo;          \
     })
#elif defined(__GNUC__) && defined(__m68k__)
  #define mulu32(x,y,hi_zuweisung,lo_zuweisung)  \
    ({ var uint32 _x = (x);						\
       var uint32 _y = (y);						\
       var uint16 _x1 = high16(_x);					\
       var uint16 _x0 = low16(_x);					\
       var uint16 _y1 = high16(_y);					\
       var uint16 _y0 = low16(_y);					\
       var uint32 _hi = mulu16(_x1,_y1); /* obere Portion */		\
       var uint32 _lo = mulu16(_x0,_y0); /* untere Portion */		\
       {var uint32 _mid = mulu16(_x0,_y1); /* 1. mittlere Portion */	\
        _hi += high16(_mid); _mid = highlow32_0(low16(_mid));		\
        _lo += _mid; if (_lo < _mid) { _hi += 1; } /* 64-Bit-Addition */\
       }								\
       {var uint32 _mid = mulu16(_x1,_y0); /* 2. mittlere Portion */	\
        _hi += high16(_mid); _mid = highlow32_0(low16(_mid));		\
        _lo += _mid; if (_lo < _mid) { _hi += 1; } /* 64-Bit-Addition */\
       }								\
       unused (hi_zuweisung _hi);					\
       lo_zuweisung _lo;						\
     })
#elif defined(__GNUC__) && defined(__sparc64__) && !defined(NO_ASM)
  #define mulu32(x,y,hi_zuweisung,lo_zuweisung)  \
    ({ var register uint64 _prod;				\
       __asm__("umul %1,%2,%0"					\
	       : "=r" (_prod)					\
	       : "r" ((uint32)(x)), "r" ((uint32)(y))		\
	      );						\
       unused (hi_zuweisung (uint32)(_prod>>32));		\
       lo_zuweisung (uint32)(_prod);				\
     })
#elif defined(__GNUC__) && defined(__sparc__) && !defined(NO_ASM)
  #define mulu32(x,y,hi_zuweisung,lo_zuweisung)  \
    ({ lo_zuweisung mulu32_(x,y); /* extern in Assembler */	\
      {var register uint32 _hi __asm__("%g1");			\
       unused (hi_zuweisung _hi);				\
     }})
#elif defined(__GNUC__) && defined(__arm__) && 0 // see comment cl_asm_arm.cc
  #define mulu32(x,y,hi_zuweisung,lo_zuweisung)  \
    ({ lo_zuweisung mulu32_(x,y); /* extern in Assembler */	\
      {var register uint32 _hi __asm__("%r1"/*"%a2"*/);		\
       unused (hi_zuweisung _hi);				\
     }})
#elif defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__)) && !defined(NO_ASM)
  #define mulu32(x,y,hi_zuweisung,lo_zuweisung)  \
    ({ var register uint32 _hi;                                  \
       var register uint32 _lo;                                  \
       __asm__("mull %2"                                         \
               : "=d" /* %edx */ (_hi), "=a" /* %eax */ (_lo)    \
               : "g" ((uint32)(x)), "1" /* %eax */ ((uint32)(y)) \
              );                                                 \
       unused (hi_zuweisung _hi); lo_zuweisung _lo;              \
     })
#elif defined(__GNUC__) && defined(__mips__) && !defined(NO_ASM)
  #define mulu32(x,y,hi_zuweisung,lo_zuweisung)  \
    ({ var register uint32 _hi;                       \
       var register uint32 _lo;                       \
       __asm__("multu %3,%2 ; mfhi %0 ; mflo %1"      \
               : "=r" (_hi), "=r" (_lo)               \
               : "r" ((uint32)(x)), "r" ((uint32)(y)) \
              );                                      \
       unused (hi_zuweisung _hi); lo_zuweisung _lo;   \
     })
#elif defined(__GNUC__) && defined(HAVE_LONGLONG) && !defined(__arm__)
  #define mulu32(x,y,hi_zuweisung,lo_zuweisung)  \
    ({ var register uint64 _prod = (uint64)(uint32)(x) * (uint64)(uint32)(y); \
       unused (hi_zuweisung (uint32)(_prod>>32));                             \
       lo_zuweisung (uint32)(_prod);                                          \
     })
#elif defined(WATCOM) && defined(__i386__) && !defined(NO_ASM)
  #define mulu32(x,y,hi_zuweisung,lo_zuweisung)  \
    { var register uint32 _hi;                     \
      var register uint32 _lo;                     \
      _lo = mulu32_(x,y), _hi = mulu32_high_();    \
      unused (hi_zuweisung _hi); lo_zuweisung _lo; \
    }
  extern "C" uint32 mulu32_high_ (void);
  #pragma aux mulu32_ = 0xF7 0xE2 /* mull %edx */ parm [eax] [edx] value [eax] modify [eax edx];
  #pragma aux mulu32_high_ = /* */ value [edx] modify [];
#else
  #define mulu32(x,y,hi_zuweisung,lo_zuweisung)  \
    { lo_zuweisung mulu32_(x,y); unused (hi_zuweisung mulu32_high); }
  #if (defined(__m68k__) || defined(__sparc__) || defined(__sparc64__) || defined(__arm__) || (defined(__i386__) && !defined(WATCOM) && !defined(MICROSOFT)) || defined(__x86_64__) || defined(__mips__) || defined(__hppa__)) && !defined(NO_ASM)
    // mulu32_ extern in Assembler
    #if defined(__sparc__) || defined(__sparc64__)
      extern "C" uint32 _get_g1 (void);
      #define mulu32_high  (_get_g1()) // Rückgabe im Register %g1
    #elif !defined(__hppa__)
      #define NEED_VAR_mulu32_high
    #endif
  #else
    #define NEED_FUNCTION_mulu32_
  #endif
#endif

#ifdef HAVE_FAST_LONGLONG

// Multipliziert zwei 32-Bit-Zahlen miteinander und liefert eine 64-Bit-Zahl:
// mulu32_w(arg1,arg2)
// > arg1, arg2 : zwei 32-Bit-Zahlen
// < result : eine 64-Bit-Zahl
#if defined(__GNUC__) && defined(__sparc64__) && !defined(NO_ASM)
  // Prefer the umul instruction over the mulx instruction (overkill).
  #define mulu32_w(x,y)  \
    ({ var register uint64 _prod;				\
       __asm__("umul %1,%2,%0"					\
	       : "=r" (_prod)					\
	       : "r" ((uint32)(x)), "r" ((uint32)(y))		\
	      );						\
       _prod;							\
     })
#elif defined(__GNUC__)
  #define mulu32_w(x,y)  ((uint64)(uint32)(x) * (uint64)(uint32)(y))
#else
  extern "C" uint64 mulu32_w (uint32 arg1, uint32 arg2);
  #define NEED_FUNCTION_mulu32_w
#endif

// Multipliziert zwei 64-Bit-Zahlen miteinander und liefert eine 128-Bit-Zahl:
// mulu64(arg1,arg2,hi=,lo=);
// > arg1, arg2 : zwei 64-Bit-Zahlen
// < 2^64*hi+lo : eine 128-Bit-Zahl
  extern "C" uint64 mulu64_ (uint64 arg1, uint64 arg2); // -> Low-Teil
#ifdef _MSC_VER
  // Workaround MSVC compiler bug.
} extern "C" uint64 mulu64_high; namespace cln {        // -> High-Teil
#else
  extern "C" uint64 mulu64_high;                        // -> High-Teil
#endif
#if defined(__GNUC__) && defined(__alpha__) && !defined(NO_ASM)
  #define mulu64(x,y,hi_zuweisung,lo_zuweisung)  \
    ({ var register uint64 _x = (x);	\
       var register uint64 _y = (y);	\
       var register uint64 _hi;		\
       var register uint64 _lo;		\
       __asm__("mulq %1,%2,%0"		\
               : "=r" (_lo)		\
               : "r" (_x), "r" (_y)	\
              );			\
       __asm__("umulh %1,%2,%0"		\
               : "=r" (_hi)		\
               : "r" (_x), "r" (_y)	\
              );			\
       hi_zuweisung _hi;		\
       lo_zuweisung _lo;		\
     })
#elif defined(__GNUC__) && defined(__sparc64__) && !defined(NO_ASM)
  #define mulu64(x,y,hi_zuweisung,lo_zuweisung)  \
    ({ lo_zuweisung mulu64_(x,y); /* extern in Assembler */	\
      {var register uint64 _hi __asm__("%g2");			\
       hi_zuweisung _hi;					\
     }})
#elif defined(__GNUC__) && defined(__x86_64__) && !defined(NO_ASM)
  #define mulu64(x,y,hi_zuweisung,lo_zuweisung)	 \
    ({ var register uint64 _hi;                                  \
       var register uint64 _lo;                                  \
       __asm__("mulq %2"                                         \
               : "=d" /* %rdx */ (_hi), "=a" /* %rax */ (_lo)    \
               : "rm" ((uint64)(x)), "1" /* %rax */ ((uint64)(y)) \
              );                                                 \
       hi_zuweisung _hi; lo_zuweisung _lo;                       \
     })
#elif defined(__GNUC__) && defined(__ia64__) && !defined(NO_ASM)
  #define mulu64(x,y,hi_zuweisung,lo_zuweisung)	 \
    ({ var register uint64 _x = (x);				  \
       var register uint64 _y = (y);				  \
       var register uint64 _hi;					  \
       __asm__("xma.hu %0 = %1, %2, f0"				  \
               : "=f" (_hi)					  \
               : "f" ((uint64)(_x)), "f" ((uint64)(_y))		  \
              );						  \
       hi_zuweisung _hi; lo_zuweisung ((uint64)(_x)*(uint64)(_y));\
     })
#else
  #define mulu64(x,y,hi_zuweisung,lo_zuweisung)  \
    { lo_zuweisung mulu64_(x,y); hi_zuweisung mulu64_high; }
  #if defined(__sparc64__) && !defined(NO_ASM)
    // mulu64_ extern in Assembler
    extern "C" uint64 _get_g2 (void);
    #define mulu64_high  (_get_g2()) // Rückgabe im Register %g2
  #else
    #define NEED_FUNCTION_mulu64_
  #endif
#endif

#endif /* HAVE_FAST_LONGLONG */


// Dividiert eine 16-Bit-Zahl durch eine 16-Bit-Zahl und
// liefert einen 16-Bit-Quotienten und einen 16-Bit-Rest.
// divu_1616_1616(x,y,q=,r=);
// > uint16 x: Zähler
// > uint16 y: Nenner
// < uint16 q: floor(x/y)
// < uint16 r: x mod y
// < x = q*y+r
  #define divu_1616_1616(x,y,q_zuweisung,r_zuweisung)  \
    { var uint16 __x = (x);					\
      var uint16 __y = (y);					\
      q_zuweisung floor(__x,__y);				\
      r_zuweisung (__x % __y);					\
    }

// Dividiert eine 32-Bit-Zahl durch eine 16-Bit-Zahl und
// liefert einen 16-Bit-Quotienten und einen 16-Bit-Rest.
// divu_3216_1616(x,y,q=,r=);
// > uint32 x: Zähler
// > uint16 y: Nenner
// > Es sei bekannt, daß 0 <= x < 2^16*y .
// < uint16 q: floor(x/y)
// < uint16 r: x mod y
// < x = q*y+r
#if defined(__sparc__)
  extern "C" uint32 divu_3216_1616_ (uint32 x, uint16 y); // -> Quotient q, Rest r
#else
  extern "C" uint16 divu_3216_1616_ (uint32 x, uint16 y); // -> Quotient q
#ifdef _MSC_VER
  // Workaround MSVC compiler bug.
} extern "C" uint16 divu_16_rest; namespace cln {         // -> Rest r
#else
  extern "C" uint16 divu_16_rest;                         // -> Rest r
#endif
#endif
#if defined(__GNUC__) && defined(__sparc64__) && !defined(NO_ASM)
  #define divu_3216_1616(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint32 __x = (x);        \
      var uint16 __y = (y);        \
      var uint64 __q;              \
      var uint64 __r;              \
      __asm__ __volatile__ (       \
        "wr %%g0,%%g0,%%y\n\t"     \
        "udiv %2,%3,%0\n\t"        \
        "umul %0,%3,%1\n\t"        \
        "sub %2,%1,%1"             \
        : "=&r" (__q), "=&r" (__r) \
        : "r" (__x), "r" (__y));   \
      q_zuweisung (uint16)__q;     \
      r_zuweisung (uint16)__r;     \
     })
#elif defined(__GNUC__) && (defined(__sparc__) || defined(__sparc64__)) && !defined(NO_ASM)
  #define divu_3216_1616(x,y,q_zuweisung,r_zuweisung)  \
    ({ var uint32 __qr = divu_3216_1616_(x,y); /* extern in Assembler */\
       q_zuweisung low16(__qr);						\
       r_zuweisung high16(__qr);					\
     })
#elif defined(__GNUC__) && defined(__m68k__) && !defined(NO_ASM)
  #define divu_3216_1616(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint32 __x = (x);						\
      var uint16 __y = (y);						\
      var uint32 __qr;							\
      __asm__ __volatile__ ("						\
        divu %2,%0							\
        " : "=d" (__qr) : "0" (__x), "dm" (__y));			\
      q_zuweisung low16(__qr);						\
      r_zuweisung high16(__qr);						\
     })
#elif defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__)) && !defined(NO_ASM)
  #define divu_3216_1616(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint32 __x = (x);						\
      var uint16 __y = (y);						\
      var uint16 __q;							\
      var uint16 __r;							\
      __asm__("divw %4"							\
              : "=a" /* %ax */ (__q), "=d" /* %dx */ (__r)		\
              : "1" /* %dx */ ((uint16)(high16(__x))), "0" /* %ax */ ((uint16)(low16(__x))), "rm" (__y) \
             );								\
      q_zuweisung __q;							\
      r_zuweisung __r;							\
     })
#elif defined(__GNUC__) && defined(__arm__) && 0 // see comment cl_asm_arm.cc
  #define divu_3216_1616(x,y,q_zuweisung,r_zuweisung)  \
    { var uint32 __q = divu_3216_1616_(x,y); /* extern in Assembler */	\
      var register uint32 __r __asm__("%r1"/*"%a2"*/);			\
      q_zuweisung __q; r_zuweisung __r;					\
    }
#elif defined(__GNUC__) && !defined(__arm__)
  #define divu_3216_1616(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint32 __x = (x);						\
      var uint16 __y = (y);						\
      var uint16 __q = floor(__x,__y);					\
      q_zuweisung __q;							\
      r_zuweisung (__x - __q * __y);					\
     })
#elif (defined(__sparc__) || defined(__sparc64__)) && !defined(NO_ASM)
  #define divu_3216_1616(x,y,q_zuweisung,r_zuweisung)  \
    { var uint32 __qr = divu_3216_1616_(x,y); /* extern in Assembler */	\
      q_zuweisung low16(__qr);						\
      r_zuweisung high16(__qr);						\
    }
#elif defined(__arm__) && !defined(NO_ASM)
  #define divu_3216_1616(x,y,q_zuweisung,r_zuweisung)  \
    { q_zuweisung divu_3216_1616_(x,y); /* extern in Assembler */	\
      r_zuweisung divu_16_rest;						\
    }
  #define NEED_VAR_divu_16_rest
#else
  #define divu_3216_1616(x,y,q_zuweisung,r_zuweisung)  \
    { q_zuweisung divu_3216_1616_(x,y); r_zuweisung divu_16_rest; }
  #define NEED_FUNCTION_divu_3216_1616_
#endif

// Dividiert eine 32-Bit-Zahl durch eine 16-Bit-Zahl und
// liefert einen 32-Bit-Quotienten und einen 16-Bit-Rest.
// divu_3216_3216(x,y,q=,r=);
// > uint32 x: Zähler
// > uint16 y: Nenner
// Es sei bekannt, daß y>0.
// < uint32 q: floor(x/y)
// < uint16 r: x mod y
// < x = q*y+r
  extern "C" uint32 divu_3216_3216_ (uint32 x, uint16 y); // -> Quotient q
#ifdef _MSC_VER
  // Workaround MSVC compiler bug.
} extern "C" uint16 divu_16_rest; namespace cln {         // -> Rest r
#else
  extern "C" uint16 divu_16_rest;                         // -> Rest r
#endif
#if defined(__GNUC__) && defined(__sparc64__) && !defined(NO_ASM)
  #define divu_3216_3216(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint32 __x = (x);        \
      var uint16 __y = (y);        \
      var uint64 __q;              \
      var uint64 __r;              \
      __asm__ __volatile__ (       \
        "wr %%g0,%%g0,%%y\n\t"     \
        "udiv %2,%3,%0\n\t"        \
        "umul %0,%3,%1\n\t"        \
        "sub %2,%1,%1"             \
        : "=&r" (__q), "=&r" (__r) \
        : "r" (__x), "r" (__y));   \
      q_zuweisung (uint32)__q;     \
      r_zuweisung (uint16)__r;     \
     })
#elif defined(__sparc__) || defined(__sparc64__) || defined(__i386__) || defined(__x86_64__)
  #define divu_3216_3216  divu_3232_3232
#else
  // Methode: (beta = 2^16)
  // x = x1*beta+x0 schreiben.
  // Division mit Rest: x1 = q1*y + r1, wobei 0 <= x1 < beta <= beta*y.
  // Also 0 <= q1 < beta, 0 <= r1 < y.
  // Division mit Rest: (r1*beta+x0) = q0*y + r0, wobei 0 <= r1*beta+x0 < beta*y.
  // Also 0 <= q0 < beta, 0 <= r0 < y
  // und x = x1*beta+x0 = (q1*beta+q0)*y + r0.
  // Setze q := q1*beta+q0 und r := r0.
  #define divu_3216_3216(x,y,q_zuweisung,r_zuweisung)  \
    { var uint32 _x = (x);						\
      var uint16 _y = (y);						\
      var uint16 _q1;							\
      var uint16 _q0;							\
      var uint16 _r1;							\
      divu_3216_1616(high16(_x),_y, _q1 = , _r1 = );			\
      divu_3216_1616(highlow32(_r1,low16(_x)),_y, _q0 = , r_zuweisung);	\
      q_zuweisung highlow32(_q1,_q0);					\
    }
#endif

// Dividiert eine 32-Bit-Zahl durch eine 32-Bit-Zahl und
// liefert einen 32-Bit-Quotienten und einen 32-Bit-Rest.
// divu_3232_3232(x,y,q=,r=);
// > uint32 x: Zähler
// > uint32 y: Nenner
// Es sei bekannt, daß y>0.
// < uint32 q: floor(x/y)
// < uint32 r: x mod y
// < x = q*y+r
  extern "C" uint32 divu_3232_3232_ (uint32 x, uint32 y); // -> Quotient q
#ifdef _MSC_VER
  // Workaround MSVC compiler bug.
} extern "C" uint32 divu_32_rest; namespace cln {         // -> Rest r
#else
  extern "C" uint32 divu_32_rest;                         // -> Rest r
#endif
#if defined(__GNUC__) && defined(__sparc64__) && !defined(NO_ASM)
  #define divu_3232_3232(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint32 __x = (x);        \
      var uint32 __y = (y);        \
      var uint64 __q;              \
      var uint64 __r;              \
      __asm__ __volatile__ (       \
        "wr %%g0,%%g0,%%y\n\t"     \
        "udiv %2,%3,%0\n\t"        \
        "umul %0,%3,%1\n\t"        \
        "sub %2,%1,%1"             \
        : "=&r" (__q), "=&r" (__r) \
        : "r" (__x), "r" (__y));   \
      q_zuweisung (uint32)__q;     \
      r_zuweisung (uint32)__r;     \
     })
  #define divu_3232_3232_(x,y) divu_6432_3232_(0,x,y)
#elif defined(__sparc__) || defined(__sparc64__) || defined(__i386__) || defined(__x86_64__)
  #define divu_3232_3232(x,y,q_zuweisung,r_zuweisung)  \
    divu_6432_3232(0,x,y,q_zuweisung,r_zuweisung)
  #define divu_3232_3232_(x,y) divu_6432_3232_(0,x,y)
#else
  // Methode: (beta = 2^n = 2^16, n = 16)
  // Falls y < beta, handelt es sich um eine 32-durch-16-Bit-Division.
  // Falls y >= beta:
  // Quotient  q = floor(x/y) < beta  (da 0 <= x < beta^2, y >= beta).
  // y habe genau n+k Bits (1 <= k <= n), d.h. 2^(n+k-1) <= y < 2^(n+k).
  // Schreibe  x = 2^k*x1 + x0  mit  x1 := floor(x/2^k)
  // und       y = 2^k*y1 + y0  mit  y1 := floor(y/2^k)
  // und bilde den Näherungs-Quotienten floor(x1/y1)
  // oder (noch besser) floor(x1/(y1+1)).
  // Wegen 0 <= x1 < 2^(2n) und 0 < 2^(n-1) <= y1 < 2^n
  // und  x1/(y1+1) <= x/y < x1/(y1+1) + 2
  // (denn x1/(y1+1) = (x1*2^k)/((y1+1)*2^k) <= (x1*2^k)/y <= x/y
  // und x/y - x1/(y1+1) = (x+x*y1-x1*y)/(y*(y1+1))
  // = (x+x0*y1-x1*y0)/(y*(y1+1)) <= (x+x0*y1)/(y*(y1+1))
  // <= x/(y*(y1+1)) + x0/y
  // <= 2^(2n)/(2^(n+k-1)*(2^(n-1)+1)) + 2^k/2^(n+k-1)
  // = 2^(n-k+1)/(2^(n-1)+1) + 2^(1-n) <= 2^n/(2^(n-1)+1) + 2^(1-n) < 2 )
  // gilt  floor(x1/(y1+1)) <= floor(x/y) <= floor(x1/(y1+1)) + 2  .
  // Man bildet also  q:=floor(x1/(y1+1))  (ein Shift um n Bit oder
  // eine (2n)-durch-n-Bit-Division, mit Ergebnis q <= floor(x/y) < beta)
  // und x-q*y und muß hiervon noch höchstens 2 mal y abziehen und q
  // incrementieren, um den Quotienten  q = floor(x/y)  und den Rest
  // x-floor(x/y)*y  der Division zu bekommen.
  #define divu_3232_3232(x,y,q_zuweisung,r_zuweisung)  \
    { var uint32 _x = (x);						\
      var uint32 _y = (y);						\
      if (_y <= (uint32)(bit(16)-1))					\
        { var uint16 _q1;						\
          var uint16 _q0;						\
          var uint16 _r1;						\
          divu_3216_1616(high16(_x),_y, _q1 = , _r1 = );		\
          divu_3216_1616(highlow32(_r1,low16(_x)),_y, _q0 = , r_zuweisung); \
          q_zuweisung highlow32(_q1,_q0);				\
        }								\
        else								\
        { var uint32 _x1 = _x; /* x1 := x */				\
          var uint32 _y1 = _y; /* y1 := y */				\
          var uint16 _q;						\
          do { _x1 = floor(_x1,2); _y1 = floor(_y1,2); } /* k erhöhen */\
             until (_y1 <= (uint32)(bit(16)-1)); /* bis y1 < beta */	\
          { var uint16 _y2 = low16(_y1)+1; /* y1+1 bilden */		\
            if (_y2==0)							\
              { _q = high16(_x1); } /* y1+1=beta -> ein Shift */	\
              else							\
              { divu_3216_1616(_x1,_y2,_q=,); } /* Division von x1 durch y1+1 */\
          }								\
          /* _q = q = floor(x1/(y1+1)) */				\
          /* x-q*y bilden (eine 16-mal-32-Bit-Multiplikation ohne Überlauf): */\
          _x -= highlow32_0(mulu16(_q,high16(_y))); /* q * high16(y) * beta */\
          /* gefahrlos, da q*high16(y) <= q*y/beta <= x/beta < beta */	\
          _x -= mulu16(_q,low16(_y)); /* q * low16(y) */		\
          /* gefahrlos, da q*high16(y)*beta + q*low16(y) = q*y <= x */	\
          /* Noch höchstens 2 mal y abziehen: */			\
          if (_x >= _y)							\
            { _q += 1; _x -= _y;					\
              if (_x >= _y)						\
                { _q += 1; _x -= _y; }					\
            }								\
          r_zuweisung _x;						\
          q_zuweisung (uint32)(_q);					\
    }   }
  #define NEED_FUNCTION_divu_3232_3232_
#endif

// Dividiert eine 64-Bit-Zahl durch eine 32-Bit-Zahl und
// liefert einen 32-Bit-Quotienten und einen 32-Bit-Rest.
// divu_6432_3232(xhi,xlo,y,q=,r=);
// > uint32 xhi,xlo: x = 2^32*xhi+xlo = Zähler
// > uint32 y: Nenner
// > Es sei bekannt, daß 0 <= x < 2^32*y .
// < uint32 q: floor(x/y)
// < uint32 r: x mod y
// < x = q*y+r
  extern "C" uint32 divu_6432_3232_ (uint32 xhi, uint32 xlo, uint32 y); // -> Quotient q
#ifdef _MSC_VER
  // Workaround MSVC compiler bug.
} extern "C" uint32 divu_32_rest; namespace cln {                       // -> Rest r
#else
  extern "C" uint32 divu_32_rest;                                       // -> Rest r
#endif
#if defined(__GNUC__) && defined(__m68k__) && !defined(NO_ASM)
  #define divu_6432_3232(xhi,xlo,y,q_zuweisung,r_zuweisung)  \
    ({var uint32 __xhi = (xhi);						\
      var uint32 __xlo = (xlo);						\
      var uint32 __y = (y);						\
      var uint32 __q;							\
      var uint32 __r;							\
      __asm__ __volatile__ ("						\
        divul %4,%1:%0							\
        " : "=d" (__q), "=d" (__r) : "1" (__xhi), "0" (__xlo), "dm" (__y)); \
      q_zuweisung __q;							\
      r_zuweisung __r;							\
     })
  #define divu_6432_3232_(xhi,xlo,y) \
    ({var uint32 ___q; divu_6432_3232(xhi,xlo,y,___q=,); ___q; })
#elif defined(__GNUC__) && defined(__sparc64__) && !defined(NO_ASM)
  #define divu_6432_3232(xhi,xlo,y,q_zuweisung,r_zuweisung)  \
    ({var uint32 __xhi = (xhi);    \
      var uint32 __xlo = (xlo);    \
      var uint32 __y = (y);        \
      var uint64 __q;              \
      var uint64 __r;              \
      __asm__ __volatile__ (       \
        "wr %2,%%g0,%%y\n\t"       \
        "udiv %3,%4,%0\n\t"        \
        "umul %0,%4,%1\n\t"        \
        "sub %3,%1,%1"             \
        : "=&r" (__q), "=&r" (__r) \
        : "r" (__xhi), "r" (__xlo), "r" (__y)); \
      q_zuweisung (uint32)__q;     \
      r_zuweisung (uint32)__r;     \
     })
#elif defined(__GNUC__) && (defined(__sparc__) || defined(__sparc64__)) && !defined(NO_ASM)
  #define divu_6432_3232(xhi,xlo,y,q_zuweisung,r_zuweisung)  \
    ({ var uint32 _q = divu_6432_3232_(xhi,xlo,y); /* extern in Assembler */\
       var register uint32 _r __asm__("%g1");				    \
       q_zuweisung _q; r_zuweisung _r;					    \
     })
#elif defined(__GNUC__) && defined(__arm__) && 0 // see comment cl_asm_arm.cc
  #define divu_6432_3232(xhi,xlo,y,q_zuweisung,r_zuweisung)  \
    ({ var uint32 _q = divu_6432_3232_(xhi,xlo,y); /* extern in Assembler */\
       var register uint32 _r __asm__("%r1"/*"%a2"*/);			    \
       q_zuweisung _q; r_zuweisung _r;					    \
     })
#elif defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__)) && !defined(NO_ASM)
  #define divu_6432_3232(xhi,xlo,y,q_zuweisung,r_zuweisung)  \
    ({var uint32 __xhi = (xhi);						\
      var uint32 __xlo = (xlo);						\
      var uint32 __y = (y);						\
      var uint32 __q;							\
      var uint32 __r;							\
      __asm__ __volatile__ (						\
         "divl %4"							\
         : "=a" /* %eax */ (__q), "=d" /* %edx */ (__r)			\
         : "1" /* %edx */ (__xhi), "0" /* %eax */ (__xlo), "rm" (__y)	\
         );								\
      q_zuweisung __q;							\
      r_zuweisung __r;							\
     })
  #define divu_6432_3232_(xhi,xlo,y) \
    ({var uint32 ___q; divu_6432_3232(xhi,xlo,y,___q=,); ___q; })
#elif defined(__GNUC__) && defined(HAVE_LONGLONG) && !defined(__arm__)
  #define divu_6432_3232(xhi,xlo,y,q_zuweisung,r_zuweisung) \
    ({var uint32 __xhi = (xhi);						\
      var uint32 __xlo = (xlo);						\
      var uint64 __x = ((uint64)__xhi << 32) | (uint64)__xlo;		\
      var uint32 __y = (y);						\
      var uint32 __q = floor(__x,(uint64)__y);				\
      q_zuweisung __q; r_zuweisung __xlo - __q * __y;			\
     })
  #define divu_6432_3232_(xhi,xlo,y) \
    ({var uint32 ___q; divu_6432_3232(xhi,xlo,y,___q=,); ___q; })
#elif defined(WATCOM) && defined(__i386__) && !defined(NO_ASM)
  #define divu_6432_3232(xhi,xlo,y,q_zuweisung,r_zuweisung)  \
    { var uint32 __xhi = (xhi);						\
      var uint32 __xlo = (xlo);						\
      var uint32 __y = (y);						\
      var uint32 __q;							\
      var uint32 __r;							\
      __q = divu_6432_3232_(__xhi,__xlo,__y); __r = divu_6432_3232_rest(); \
      q_zuweisung __q;							\
      r_zuweisung __r;							\
    }
  extern "C" uint32 divu_6432_3232_rest (void);
  #pragma aux divu_6432_3232_ = 0xF7 0xF1 /* divl %ecx */ parm [edx] [eax] [ecx] value [eax] modify [eax edx];
  #pragma aux divu_6432_3232_rest = /* */ value [edx] modify [];
#else
  #define divu_6432_3232(xhi,xlo,y,q_zuweisung,r_zuweisung)  \
    { q_zuweisung divu_6432_3232_(xhi,xlo,y); r_zuweisung divu_32_rest; }
  #if (defined(__m68k__) || defined(__sparc__) || defined(__sparc64__) || defined(__arm__) || (defined(__i386__) && !defined(WATCOM) && !defined(MICROSOFT)) || defined(__x86_64__) || defined(__hppa__)) && !defined(NO_ASM)
    // divu_6432_3232_ extern in Assembler
    #if defined(__sparc__) || defined(__sparc64__)
      extern "C" uint32 _get_g1 (void);
      #define divu_32_rest  (_get_g1()) // Rückgabe im Register %g1
    #else
      #define NEED_VAR_divu_32_rest
    #endif
  #else
    #define NEED_FUNCTION_divu_6432_3232_
  #endif
#endif

#ifdef HAVE_FAST_LONGLONG

// Dividiert eine 64-Bit-Zahl durch eine 32-Bit-Zahl und
// liefert einen 32-Bit-Quotienten und einen 32-Bit-Rest.
// divu_6432_3232_w(x,y,q=,r=);
// > uint64 x: Zähler
// > uint32 y: Nenner
// > Es sei bekannt, daß 0 <= x < 2^32*y .
// < uint32 q: floor(x/y)
// < uint32 r: x mod y
// < x = q*y+r
#if defined(__GNUC__) && defined(__sparc64__) && !defined(NO_ASM)
  // Prefer the udiv and umul instructions over the udivx and mulx instructions
  // (overkill).
  #define divu_6432_3232_w(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint64 __x = (x);           \
      var uint32 __xhi = high32(__x); \
      var uint32 __xlo = low32(__x);  \
      var uint32 __y = (y);           \
      var uint64 __q;                 \
      var uint64 __r;                 \
      __asm__ __volatile__ (          \
        "wr %2,%%g0,%%y\n\t"          \
        "udiv %3,%4,%0\n\t"           \
        "umul %0,%4,%1\n\t"           \
        "sub %3,%1,%1"                \
        : "=&r" (__q), "=&r" (__r)    \
        : "r" (__xhi), "r" (__xlo), "r" (__y)); \
      q_zuweisung (uint32)__q;        \
      r_zuweisung (uint32)__r;        \
     })
#elif defined(__GNUC__) && (defined(__alpha__) || defined(__ia64__) || defined(__mips64__) || defined(__sparc64__))
  // On __alpha__, computing the remainder by multiplication is just two
  // instructions, compared to the __remqu (libc) function call for the %
  // operator.
  // On __ia64__, computing the remainder by multiplication is just four
  // instructions, compared to the __umoddi3 (libgcc) function call for the %
  // operator.
  // On __mips64__, computing the remainder by multiplication is just two
  // instructions, compared to the __umoddi3 (libgcc) function call for the %
  // operator.
  // On __sparc64__, computing the remainder by multiplication uses a 32-bit
  // multiplication instruction, compared to a 64-bit multiplication when the %
  // operator is used.
  #define divu_6432_3232_w(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint64 __x = (x);						\
      var uint32 __y = (y);						\
      var uint32 __q = floor(__x,(uint64)__y);				\
      q_zuweisung __q; r_zuweisung (uint32)__x - __q * __y;		\
     })
#elif defined(__GNUC__) && defined(__x86_64__)
  // On __x86_64__, gcc 4.0 performs both quotient and remainder computation
  // in a single instruction.
  #define divu_6432_3232_w(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint64 __x = (x);						\
      var uint32 __y = (y);						\
      var uint32 __q = floor(__x,(uint64)__y);				\
      q_zuweisung __q; r_zuweisung __x % (uint64)__y;			\
     })
#else
  #define divu_6432_3232_w(x,y,q_zuweisung,r_zuweisung)  \
    { var uint64 __x = (x);						  \
      divu_6432_3232(high32(__x),low32(__x),(y),q_zuweisung,r_zuweisung); \
    }
#endif

// Dividiert eine 64-Bit-Zahl durch eine 32-Bit-Zahl und
// liefert einen 64-Bit-Quotienten und einen 32-Bit-Rest.
// divu_6432_6432(x,y,q=,r=);
// > uint64 x: Zähler
// > uint32 y: Nenner
// > Es sei bekannt, daß y>0.
// < uint64 q: floor(x/y)
// < uint32 r: x mod y
// < x = q*y+r
#if defined(__GNUC__) && (defined(__alpha__) || defined(__ia64__) || defined(__mips64__) || defined(__sparc64__))
  // On __alpha__, computing the remainder by multiplication is just two
  // instructions, compared to the __remqu (libc) function call for the %
  // operator.
  // On __ia64__, computing the remainder by multiplication is just four
  // instructions, compared to the __umoddi3 (libgcc) function call for the %
  // operator.
  // On __mips64__, computing the remainder by multiplication is just two
  // instructions, compared to the __umoddi3 (libgcc) function call for the %
  // operator.
  // On __sparc64__, computing the remainder by multiplication uses a 32-bit
  // multiplication instruction, compared to a 64-bit multiplication when the %
  // operator is used.
  #define divu_6432_6432(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint64 _x = (x);                    \
      var uint32 _y = (y);                    \
      var uint64 _q;                          \
      q_zuweisung _q = floor(_x,(uint64)_y);  \
      r_zuweisung low32(_x) - low32(_q) * _y; \
     })
#elif defined(__GNUC__) && defined(__x86_64__)
  // On __x86_64__, gcc 4.0 performs both quotient and remainder computation
  // in a single instruction.
  #define divu_6432_6432(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint64 _x = (x);               \
      var uint32 _y = (y);               \
      q_zuweisung floor(_x,(uint64)_y);  \
      r_zuweisung _x % (uint64)_y;       \
     })
#else
  // Methode: (beta = 2^32)
  // x = x1*beta+x0 schreiben.
  // Division mit Rest: x1 = q1*y + r1, wobei 0 <= x1 < beta <= beta*y.
  // Also 0 <= q1 < beta, 0 <= r1 < y.
  // Division mit Rest: (r1*beta+x0) = q0*y + r0, wobei 0 <= r1*beta+x0 < beta*y.
  // Also 0 <= q0 < beta, 0 <= r0 < y
  // und x = x1*beta+x0 = (q1*beta+q0)*y + r0.
  // Setze q := q1*beta+q0 und r := r0.
  #if defined(__GNUC__)
    #define divu_6432_6432(x,y,q_zuweisung,r_zuweisung)  \
      ({var uint64 _x = (x);            \
        var uint32 _y = (y);            \
        var uint32 _q1;                 \
        var uint32 _q0;                 \
        var uint32 _r1;                 \
        divu_6432_3232(0,high32(_x),_y, _q1 = , _r1 = ); \
        divu_6432_3232(_r1,low32(_x),_y, _q0 = , r_zuweisung); \
        q_zuweisung highlow64(_q1,_q0); \
       })
  #else
    #define divu_6432_6432(x,y,q_zuweisung,r_zuweisung)  \
      {var uint64 _x = (x);            \
       var uint32 _y = (y);            \
       var uint32 _q1;                 \
       var uint32 _q0;                 \
       var uint32 _r1;                 \
       divu_6432_3232(0,high32(_x),_y, _q1 = , _r1 = ); \
       divu_6432_3232(_r1,low32(_x),_y, _q0 = , r_zuweisung); \
       q_zuweisung highlow64(_q1,_q0); \
      }
  #endif
#endif

// Dividiert eine 64-Bit-Zahl durch eine 64-Bit-Zahl und
// liefert einen 64-Bit-Quotienten und einen 64-Bit-Rest.
// divu_6464_6464(x,y,q=,r=);
// > uint64 x: Zähler
// > uint64 y: Nenner
// > Es sei bekannt, daß y>0.
// < uint64 q: floor(x/y)
// < uint64 r: x mod y
// < x = q*y+r
#if defined(__GNUC__) && (defined(__alpha__) || defined(__ia64__) || defined(__mips64__) || defined(__sparc64__))
  // On __alpha__, computing the remainder by multiplication is just two
  // instructions, compared to the __remqu (libc) function call for the %
  // operator.
  // On __ia64__, computing the remainder by multiplication is just four
  // instructions, compared to the __umoddi3 (libgcc) function call for the %
  // operator.
  // On __mips64__, computing the remainder by multiplication is just two
  // instructions, compared to the __umoddi3 (libgcc) function call for the %
  // operator.
  // On __sparc64__, it doesn't matter.
  #define divu_6464_6464(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint64 _x = (x);           \
      var uint64 _y = (y);           \
      var uint64 _q;                 \
      q_zuweisung _q = floor(_x,_y); \
      r_zuweisung _x - _q * _y;      \
     })
#elif defined(__GNUC__) && (defined(__sparc64__) || defined(__x86_64__))
  // On __sparc64__, it doesn't matter.
  // On __x86_64__, gcc 4.0 performs both quotient and remainder computation
  // in a single instruction.
  #define divu_6464_6464(x,y,q_zuweisung,r_zuweisung)  \
    ({var uint64 _x = (x);      \
      var uint64 _y = (y);      \
      q_zuweisung floor(_x,_y); \
      r_zuweisung _x % _y;      \
     })
#else
  // For unknown CPUs, we don't know whether gcc's __udivdi3 function plus a
  // multiplication is slower or faster than our own divu_6464_6464_ routine.
  // Anyway, call our own routine.
  extern "C" uint64 divu_6464_6464_ (uint64 x, uint64 y); // -> Quotient q
#ifdef _MSC_VER
  // Workaround MSVC compiler bug.
} extern "C" uint64 divu_64_rest; namespace cln {         // -> Rest r
#else
  extern "C" uint64 divu_64_rest;                         // -> Rest r
#endif
  #define divu_6464_6464(x,y,q_zuweisung,r_zuweisung)  \
    { q_zuweisung divu_6464_6464_(x,y); r_zuweisung divu_64_rest; }
  #define NEED_VAR_divu_64_rest
  #define NEED_FUNCTION_divu_6464_6464_
#endif

// Dividiert eine 128-Bit-Zahl durch eine 64-Bit-Zahl und
// liefert einen 64-Bit-Quotienten und einen 64-Bit-Rest.
// divu_12864_6464(xhi,xlo,y,q=,r=);
// > uint64 xhi,xlo: x = 2^64*xhi+xlo = Zähler
// > uint64 y: Nenner
// > Es sei bekannt, daß 0 <= x < 2^64*y .
// < uint64 q: floor(x/y)
// < uint64 r: x mod y
// < x = q*y+r
  extern "C" uint64 divu_12864_6464_ (uint64 xhi, uint64 xlo, uint64 y); // -> Quotient q
#ifdef _MSC_VER
  // Workaround MSVC compiler bug.
} extern "C" uint64 divu_64_rest; namespace cln {                        // -> Rest r
#else
  extern "C" uint64 divu_64_rest;                                        // -> Rest r
#endif
#if defined(__GNUC__) && defined(__x86_64__) && !defined(NO_ASM)
  #define divu_12864_6464(xhi,xlo,y,q_zuweisung,r_zuweisung)  \
    ({var uint64 __xhi = (xhi);						\
      var uint64 __xlo = (xlo);						\
      var uint64 __y = (y);						\
      var uint64 __q;							\
      var uint64 __r;							\
      __asm__ __volatile__ (						\
         "divq %4"							\
         : "=a" /* %rax */ (__q), "=d" /* %rdx */ (__r)			\
         : "1" /* %rdx */ (__xhi), "0" /* %rax */ (__xlo), "rm" (__y)	\
         );								\
      q_zuweisung __q;							\
      r_zuweisung __r;							\
     })
  #define divu_12864_64364_(xhi,xlo,y) \
    ({var uint64 ___q; divu_12864_6464(xhi,xlo,y,___q=,); ___q; })
#else
  #define divu_12864_6464(xhi,xlo,y,q_zuweisung,r_zuweisung)  \
    { q_zuweisung divu_12864_6464_(xhi,xlo,y); r_zuweisung divu_64_rest; }
  #define NEED_VAR_divu_64_rest
  #define NEED_FUNCTION_divu_12864_6464_
#endif

#endif /* HAVE_FAST_LONGLONG */


// Zieht die Ganzzahl-Wurzel aus einer 32-Bit-Zahl und
// liefert eine 16-Bit-Wurzel und einen Rest.
// isqrt_32_16(x,y=,sqrtp=);
// > uint32 x: Radikand, >= 2^30, < 2^32
// < uint16 y: floor(sqrt(x)), >= 2^15, < 2^16
// < boolean sqrtp: /=0, falls x=y^2
  // Methode:
  // y := 2^16 als Anfangswert,
  // y := floor((y + floor(x/y))/2) als nächster Wert,
  // solange z := floor(x/y) < y, setze y := floor((y+z)/2).
  // y ist fertig; x=y^2 genau dann, wenn z=y und die letzte Division aufging.
  // (Beweis:
  //  1. Die Folge der y ist streng monoton fallend.
  //  2. Stets gilt y >= floor(sqrt(x)) (denn für alle y>0 ist
  //     y + x/y >= 2*sqrt(x) und daher  floor((y + floor(x/y))/2) =
  //     floor(y/2 + x/(2*y)) >= floor(sqrt(x)) ).
  //  3. Am Schluß gilt x >= y^2.
  // )
  #define isqrt_32_16(x,y_zuweisung,sqrtp_zuweisung)  \
    { var uint32 _x = (x);						\
      var uint16 _x1 = high16(_x);					\
      var uint16 _y = floor(_x1,2) | bit(16-1);				\
      loop								\
        { var uint16 _z;						\
          var uint16 _r;						\
          if (_x1 >= _y) /* Division _x/_y ergäbe Überlauf -> _z > _y */\
            { unused (sqrtp_zuweisung FALSE); break; } 			\
          divu_3216_1616(_x,_y, _z=,_r=); /* Dividiere _x/_y */		\
          if (_z >= _y)							\
            { unused (sqrtp_zuweisung (_z == _y) && (_r == 0)); break; } \
          _y = floor((uint16)(_z+_y),2) | bit(16-1); /* _y muß >= 2^15 bleiben */\
        }								\
      y_zuweisung _y;							\
    }

// Zieht die Ganzzahl-Wurzel aus einer 64-Bit-Zahl und
// liefert eine 32-Bit-Wurzel und einen Rest.
// isqrt_64_32(xhi,xlo,y=,sqrtp=);
// > uint32 xhi,xlo: Radikand x = 2^32*xhi+xlo, >= 2^62, < 2^64
// < uint32 y: floor(sqrt(x)), >= 2^31, < 2^32
// < boolean sqrtp: /=0, falls x=y^2
#if defined(__sparc__) || defined(__sparc64__) || defined(__m68k__) || defined(__hppa__)
  // Methode:
  // y := 2^32 als Anfangswert,
  // y := floor((y + floor(x/y))/2) als nächster Wert,
  // solange z := floor(x/y) < y, setze y := floor((y+z)/2).
  // y ist fertig; x=y^2 genau dann, wenn z=y und die letzte Division aufging.
  // (Beweis:
  //  1. Die Folge der y ist streng monoton fallend.
  //  2. Stets gilt y >= floor(sqrt(x)) (denn für alle y>0 ist
  //     y + x/y >= 2*sqrt(x) und daher  floor((y + floor(x/y))/2) =
  //     floor(y/2 + x/(2*y)) >= floor(sqrt(x)) ).
  //  3. Am Schluß gilt x >= y^2.
  // )
  #define isqrt_64_32(xhi,xlo,y_zuweisung,sqrtp_zuweisung)  \
    { var uint32 _xhi = (xhi);						\
      var uint32 _xlo = (xlo);						\
      var uint32 _y = floor(_xhi,2) | bit(32-1);			\
      loop								\
        { var uint32 _z;						\
          var uint32 _rest;						\
          if (_xhi >= _y) /* Division _x/_y ergäbe Überlauf -> _z > _y */\
            { sqrtp_zuweisung FALSE; break; }				\
          divu_6432_3232(_xhi,_xlo,_y, _z=,_rest=); /* Dividiere _x/_y */\
          if (_z >= _y)							\
            { sqrtp_zuweisung (_z == _y) && (_rest == 0); break; }	\
          _y = floor(_z+_y,2) | bit(32-1); /* _y muß >= 2^31 bleiben */	\
        }								\
      y_zuweisung _y;							\
    }
#else
  // Methode:
  // Wie bei UDS_sqrt mit n=2.
  // y = 2^16*yhi + ylo ansetzen.
  // Dann muß
  //   yhi = floor(y/2^16) = floor(floor(sqrt(x))/2^16)
  //       = floor(sqrt(x)/2^16) = floor(sqrt(x/2^32)) = isqrt(xhi)
  // sein. Es folgt yhi >= 2^15.
  // Danach sucht man das größte ylo >=0 mit
  // x - 2^32*yhi^2 >= 2*2^16*yhi*ylo + ylo^2.
  // Dazu setzen wir  xhi*2^32+xlo := x - 2^32*yhi^2
  // (also xhi := xhi - yhi^2, das ist >=0, <=2*yhi).
  // Die Schätzung für die zweite Ziffer
  //     ylo' := min(2^16-1,floor((xhi*2^32+xlo)/(2*2^16*yhi)))
  // erfüllt ylo'-1 <= ylo <= ylo', ist also um höchstens 1 zu groß.
  // (Beweis: Rechte Ungleichung klar, da  ylo < 2^16  und
  //   xhi*2^32+xlo >= 2*2^16*yhi*ylo + ylo^2 >= 2*2^16*yhi*ylo
  //   ==> (xhi*2^32+xlo)/(2*2^16*yhi) >= ylo  gelten muß.
  //   Linke Ungleichung: Falls floor(...)>=2^16, ist
  //   xhi*2^32+xlo >= 2*2^16*2^16*yhi >= 2*2^16*yhi*(2^16-1) + 2^32
  //                >= 2*2^16*yhi*(2^16-1) + (2^16-1)^2
  //   und xhi*2^32+xlo < 2*2^16*2^16*yhi + (2^16)^2, also
  //   ylo = 2^16-1 = ylo'.
  //   Sonst ist ylo' = floor((xhi*2^32+xlo)/(2*2^16*yhi)), also
  //   xhi*2^32+xlo >= 2*2^16*yhi*ylo' >= 2*2^16*yhi*(ylo'-1) + 2^32
  //                >= 2*2^16*yhi*(ylo'-1) + (ylo'-1)^2,
  //   also ylo >= ylo'-1 nach Definition von ylo.)
  #define isqrt_64_32(xhi,xlo,y_zuweisung,sqrtp_zuweisung)  \
    { var uint32 _xhi = (xhi);						\
      var uint32 _xlo = (xlo);						\
      var uint16 _yhi;							\
      var uint16 _ylo;							\
      /* erste Ziffer berechnen: */					\
      isqrt_32_16(_xhi,_yhi=,); /* yhi := isqrt(xhi) */			\
      _xhi -= mulu16(_yhi,_yhi); /* jetzt 0 <= xhi <= 2*yhi */		\
      /* x = 2^32*yhi^2 + 2^32*xhi + xlo */				\
      /* Schätzung für die zweite Ziffer berechnen: */			\
      /* ylo := min(2^16-1,floor((xhi*2^32+xlo)/(2*2^16*yhi))) bilden: */\
     {var uint32 _z = (_xhi << 15) | (_xlo >> 17); /* < 2^15*(2*yhi+1) */\
      var uint32 _r = highlow32_0(_yhi);				\
      if (_z >= _r)							\
        { _ylo = bit(16)-1; _r = _z - _r + (uint32)_yhi; }		\
        else								\
        { divu_3216_1616(_z,_yhi, _ylo=,_r=); }				\
      /* x = 2^32*yhi^2 + 2*2^16*yhi*ylo + 2^17*r + (xlo mod 2^17), */	\
      /* 0 <= r < yhi + 2^15 */						\
      _xlo = (_r << 17) | (_xlo & (bit(17)-1));				\
      /* x = 2^32*yhi^2 + 2*2^16*yhi*ylo + 2^32*floor(r/2^15) + xlo */	\
      _z = mulu16(_ylo,_ylo); /* z = ylo^2 */				\
      /* Versuche vom Rest 2^32*floor(r/2^15) + xlo  z zu subtrahieren. */\
      /* Falls Rest >= z (d.h. r>=2^15 oder xlo>=z), ist ylo fertig, */	\
      /* und es gilt x=y^2 genau dann, wenn r<2^15 und xlo=z. */	\
      /* Sonst (d.h. r<2^15 und xlo<z), muß man ylo erniedrigen. Dazu */\
      /* setzt man  ylo := ylo-1, z := z-(2*ylo+1), */			\
      /* Rest := Rest + 2^17*yhi = xlo + 2^17*yhi >= 2^32 > z, also x>y^2. */\
      if (_r < bit(15))							\
        { if (_xlo < _z)						\
            { _ylo -= 1; sqrtp_zuweisung FALSE; }			\
            else							\
            { sqrtp_zuweisung (_xlo == _z); }				\
        }								\
        else								\
        { sqrtp_zuweisung FALSE; }					\
      y_zuweisung highlow32(_yhi,_ylo);					\
    }}
#endif

#ifdef HAVE_FAST_LONGLONG

// Zieht die Ganzzahl-Wurzel aus einer 128-Bit-Zahl und
// liefert eine 64-Bit-Wurzel und einen Rest.
// isqrt_128_64(xhi,xlo,y=,sqrtp=);
// > uint64 xhi,xlo: Radikand x = 2^64*xhi+xlo, >= 2^126, < 2^128
// < uint64 y: floor(sqrt(x)), >= 2^63, < 2^64
// < boolean sqrtp: /=0, falls x=y^2
  // Methode:
  // Wie bei UDS_sqrt mit n=2.
  // y = 2^32*yhi + ylo ansetzen.
  // Dann muß
  //   yhi = floor(y/2^32) = floor(floor(sqrt(x))/2^32)
  //       = floor(sqrt(x)/2^32) = floor(sqrt(x/2^64)) = isqrt(xhi)
  // sein. Es folgt yhi >= 2^31.
  // Danach sucht man das größte ylo >=0 mit
  // x - 2^64*yhi^2 >= 2*2^32*yhi*ylo + ylo^2.
  // Dazu setzen wir  xhi*2^64+xlo := x - 2^64*yhi^2
  // (also xhi := xhi - yhi^2, das ist >=0, <=2*yhi).
  // Die Schätzung für die zweite Ziffer
  //     ylo' := min(2^32-1,floor((xhi*2^64+xlo)/(2*2^32*yhi)))
  // erfüllt ylo'-1 <= ylo <= ylo', ist also um höchstens 1 zu groß.
  // (Beweis: Rechte Ungleichung klar, da  ylo < 2^32  und
  //   xhi*2^64+xlo >= 2*2^32*yhi*ylo + ylo^2 >= 2*2^32*yhi*ylo
  //   ==> (xhi*2^64+xlo)/(2*2^32*yhi) >= ylo  gelten muß.
  //   Linke Ungleichung: Falls floor(...)>=2^32, ist
  //   xhi*2^64+xlo >= 2*2^32*2^32*yhi >= 2*2^32*yhi*(2^32-1) + 2^64
  //                >= 2*2^32*yhi*(2^32-1) + (2^32-1)^2
  //   und xhi*2^64+xlo < 2*2^32*2^32*yhi + (2^32)^2, also
  //   ylo = 2^32-1 = ylo'.
  //   Sonst ist ylo' = floor((xhi*2^64+xlo)/(2*2^32*yhi)), also
  //   xhi*2^64+xlo >= 2*2^32*yhi*ylo' >= 2*2^32*yhi*(ylo'-1) + 2^64
  //                >= 2*2^32*yhi*(ylo'-1) + (ylo'-1)^2,
  //   also ylo >= ylo'-1 nach Definition von ylo.)
  #define isqrt_128_64(x_hi,x_lo,y_zuweisung,sqrtp_zuweisung)  \
    { var uint64 xhi = (x_hi);						\
      var uint64 xlo = (x_lo);						\
      var uint32 yhi;							\
      var uint32 ylo;							\
      /* erste Ziffer berechnen: */					\
      isqrt_64_32(high32(xhi),low32(xhi),yhi=,); /* yhi := isqrt(xhi) */\
      xhi -= mulu32_w(yhi,yhi); /* jetzt 0 <= xhi <= 2*yhi */		\
      /* x = 2^64*yhi^2 + 2^64*xhi + xlo */				\
      /* Schätzung für die zweite Ziffer berechnen: */			\
      /* ylo := min(2^32-1,floor((xhi*2^64+xlo)/(2*2^32*yhi))) bilden: */\
     {var uint64 z = (xhi << 31) | (xlo >> 33); /* < 2^31*(2*yhi+1) */	\
      var uint64 r = highlow64_0(yhi);					\
      if (z >= r)							\
        { ylo = bit(32)-1; r = z - r + (uint64)yhi; }			\
        else								\
        { divu_6432_3232_w(z,yhi, ylo=,r=); }				\
      /* x = 2^64*yhi^2 + 2*2^32*yhi*ylo + 2^33*r + (xlo mod 2^33), */	\
      /* 0 <= r < yhi + 2^31 */						\
      xlo = (r << 33) | (xlo & (bit(33)-1));				\
      /* x = 2^64*yhi^2 + 2*2^32*yhi*ylo + 2^64*floor(r/2^31) + xlo */	\
      z = mulu32_w(ylo,ylo); /* z = ylo^2 */				\
      /* Versuche vom Rest 2^64*floor(r/2^31) + xlo  z zu subtrahieren. */\
      /* Falls Rest >= z (d.h. r>=2^31 oder xlo>=z), ist ylo fertig, */	\
      /* und es gilt x=y^2 genau dann, wenn r<2^31 und xlo=z. */	\
      /* Sonst (d.h. r<2^31 und xlo<z), muß man ylo erniedrigen. Dazu */\
      /* setzt man  ylo := ylo-1, z := z-(2*ylo+1), */			\
      /* Rest := Rest + 2^33*yhi = xlo + 2^33*yhi >= 2^64 > z, also x>y^2. */\
      if (r < bit(31))							\
        { if (xlo < z)							\
            { ylo -= 1; sqrtp_zuweisung FALSE; }			\
            else							\
            { sqrtp_zuweisung (xlo == z); }				\
        }								\
        else								\
        { sqrtp_zuweisung FALSE; }					\
      y_zuweisung highlow64(yhi,ylo);					\
    }}

#endif /* HAVE_FAST_LONGLONG */

// Zieht die Ganzzahl-Wurzel aus einer 32-Bit-Zahl und
// liefert eine 16-Bit-Wurzel.
// isqrt(x)
// > uintL x : Radikand, >=0, <2^32
// < uintL ergebnis : Wurzel, >=0, <2^16
  extern uintL isqrt (uintL x);

#ifdef HAVE_LONGLONG
// Extracts integer root of a 64-bit number and returns a 32-bit number.
// isqrt(x)
// > uintQ x : radicand, >=0, <2^64
// < uintL result : square root, >=0, <2^32
  extern uintL isqrt (uintQ x);
#endif

// Sorry for this. We need an isqrt function taking uintC arguments but we
// cannot use overloading since this would lead to ambiguities with any of the
// two signatures above.
  inline uintL isqrtC (uintC x)
  {
#if (intCsize==32)
      return isqrt((uintL)x);
#else
      return isqrt((uintQ)x);
#endif
  }


// Zieht die Ganzzahl-Wurzel aus einer 64-Bit-Zahl und
// liefert eine 32-Bit-Wurzel.
// isqrt(x1,x0)
// > uintL2 x = x1*2^32+x0 : Radikand, >=0, <2^64
// < uintL ergebnis : Wurzel, >=0, <2^32
  extern uintL isqrt (uintL x1, uintL x0);


// Bits einer 8-Bit-Zahl zählen:
// integerlength8(digit,size=);
// setzt size auf die höchste in digit vorkommende Bitnummer.
// > digit: ein uint8 >0
// < size: >0, <=8, mit 2^(size-1) <= digit < 2^size
#if defined(__GNUC__) && defined(__m68k__) && !defined(NO_ASM)
  #define integerlength8(digit,size_zuweisung)  \
    { var uintL _zero_counter; /* zählt die führenden Nullbits in digit            */\
      __asm__("bfffo %1{#0:#8},%0" : "=d" (_zero_counter) : "dm" ((uint8)(digit)) ); \
      size_zuweisung (8-_zero_counter);                                              \
    }
#elif defined(__sparc__) && !defined(__sparc64__)
  #define integerlength8(digit,size_zuweisung)  \
    integerlength32((uint32)(digit),size_zuweisung) // siehe unten
#elif defined(__GNUC__) && defined(__i386__) && !defined(NO_ASM)
  #define integerlength8(digit,size_zuweisung)  \
    integerlength16((uint16)(digit),size_zuweisung)
#else
  #define integerlength8(digit,size_zuweisung)  \
    { var uintC _bitsize = 1;					\
      var uintL _x8 = (uint8)(digit);				\
      /* _x8 hat höchstens 8 Bits.                             */\
      if (_x8 >= bit(4)) { _x8 = _x8>>4; _bitsize += 4; }		\
      /* _x8 hat höchstens 4 Bits.                             */\
      if (_x8 >= bit(2)) { _x8 = _x8>>2; _bitsize += 2; }		\
      /* _x8 hat höchstens 2 Bits.                             */\
      if (_x8 >= bit(1)) { /* _x8 = _x8>>1; */ _bitsize += 1; }	\
      /* _x8 hat höchstens 1 Bit. Dieses Bit muß gesetzt sein. */\
      size_zuweisung _bitsize;					\
    }
#endif

// Bits einer 16-Bit-Zahl zählen:
// integerlength16(digit,size=);
// setzt size auf die höchste in digit vorkommende Bitnummer.
// > digit: ein uint16 >0
// < size: >0, <=16, mit 2^(size-1) <= digit < 2^size
#if defined(__GNUC__) && defined(__m68k__) && !defined(NO_ASM)
  #define integerlength16(digit,size_zuweisung)  \
    { var uintL _zero_counter; /* zählt die führenden Nullbits in digit              */\
      __asm__("bfffo %1{#0:#16},%0" : "=d" (_zero_counter) : "dm" ((uint16)(digit)) ); \
      size_zuweisung (16-_zero_counter);                                               \
    }
#elif defined(__sparc__) && !defined(__sparc64__)
  #define integerlength16(digit,size_zuweisung)  \
    integerlength32((uint32)(digit),size_zuweisung) // siehe unten
#elif defined(__GNUC__) && defined(__i386__) && !defined(NO_ASM)
  #define integerlength16(digit,size_zuweisung)  \
    { var uintW _one_position; /* Position der führenden 1                 */\
      __asm__("bsrw %1,%0" : "=r" (_one_position) : "r" ((uint16)(digit)) ); \
      size_zuweisung (1+_one_position);                                      \
    }
// Die weiteren kommen von gcc/longlong.h :
#elif defined(__GNUC__) && defined(__ibm032__) && !defined(NO_ASM) // RT/ROMP
  #define integerlength16(digit,size_zuweisung)  \
    { var uintL _zero_counter; /* zählt die führenden Nullbits in digit   */\
      __asm__("clz %0,%1" : "=r" (_zero_counter) : "r" ((uint32)(digit)) ); \
      size_zuweisung (16-_zero_counter);                                    \
    }
#else
  #define integerlength16(digit,size_zuweisung)  \
    { var uintC _bitsize = 1;						\
      var uintWL _x16 = (uint16)(digit);					\
      /* _x16 hat höchstens 16 Bits.                                   */\
      if (_x16 >= bit(8)) { _x16 = _x16>>8; _bitsize += 8; }		\
      /* _x16 hat höchstens 8 Bits.                                    */\
      if (_x16 >= bit(4)) { _x16 = _x16>>4; _bitsize += 4; }		\
      /* _x16 hat höchstens 4 Bits.                                    */\
      if (_x16 >= bit(2)) { _x16 = _x16>>2; _bitsize += 2; }		\
      /* _x16 hat höchstens 2 Bits.                                    */\
      if (_x16 >= bit(1)) { /* _x16 = _x16>>1; */ _bitsize += 1; }		\
      /* _x16 hat höchstens 1 Bit. Dieses Bit muß gesetzt sein.        */\
      size_zuweisung _bitsize;						\
    }
#endif

// Bits einer 32-Bit-Zahl zählen:
// integerlength32(digit,size=);
// setzt size auf die höchste in digit vorkommende Bitnummer.
// > digit: ein uint32 >0
// < size: >0, <=32, mit 2^(size-1) <= digit < 2^size
#if defined(__GNUC__) && defined(__m68k__) && !defined(NO_ASM)
  #define integerlength32(digit,size_zuweisung)  \
    { var uintL _zero_counter; /* zählt die führenden Nullbits in digit              */\
      __asm__("bfffo %1{#0:#32},%0" : "=d" (_zero_counter) : "dm" ((uint32)(digit)) ); \
      size_zuweisung (32-_zero_counter);                                               \
    }
#elif defined(__sparc__) && !defined(__sparc64__) && defined(FAST_DOUBLE)
  #define integerlength32(digit,size_zuweisung)  \
    {var union { double f; uint32 i[2]; } __fi;				\
     const int df_mant_len = 52;  /* mantissa bits (excl. hidden bit) */\
     const int df_exp_mid = 1022; /* exponent bias */                   \
     /* Bilde 2^52 + digit:                                           */\
     __fi.i[0] = (uint32)(df_mant_len+1+df_exp_mid) << (df_mant_len-32); /* Vorzeichen 0, Exponent 53 */\
     __fi.i[1] = (digit); /* untere 32 Bits setzen (benutzt CL_CPU_BIG_ENDIAN_P !) */\
     /* subtrahiere 2^52:                                             */\
     __fi.f = __fi.f - (double)(4503599627370496.0L);			\
     /* Hole davon den Exponenten:                                    */\
     size_zuweisung ((__fi.i[0] >> (df_mant_len-32)) - df_exp_mid);	\
    }
#elif defined(__GNUC__) && defined(__i386__) && !defined(NO_ASM)
  #define integerlength32(digit,size_zuweisung)  \
    { var uintL _one_position; /* Position der führenden 1                  */\
      __asm__("bsrl %1,%0" : "=r" (_one_position) : "rm" ((uint32)(digit)) ); \
      size_zuweisung (1+_one_position);                                       \
    }
#elif defined(__hppa__) && !defined(NO_ASM)
  #define integerlength32(digit,size_zuweisung)  \
    size_zuweisung length32(digit);
  extern "C" uintL length32 (uintL digit); // extern in Assembler
// Die weiteren kommen von gcc/longlong.h :
#elif defined(__GNUC__) && (defined(__a29k__) || defined(___AM29K__)) && !defined(NO_ASM)
  #define integerlength32(digit,size_zuweisung)  \
    { var uintL _zero_counter; /* zählt die führenden Nullbits in digit   */\
      __asm__("clz %0,%1" : "=r" (_zero_counter) : "r" ((uint32)(digit)) ); \
      size_zuweisung (32-_zero_counter);                                    \
    }
#elif defined(__GNUC__) && defined(__gmicro__) && !defined(NO_ASM)
  #define integerlength32(digit,size_zuweisung)  \
    { var uintL _zero_counter; /* zählt die führenden Nullbits in digit      */\
      __asm__("bsch/1 %1,%0" : "=g" (_zero_counter) : "g" ((uint32)(digit)) ); \
      size_zuweisung (32-_zero_counter);                                       \
    }
#elif defined(__GNUC__) && defined(__rs6000__) && !defined(NO_ASM)
 #ifdef _AIX
  // old assembler syntax
  #define integerlength32(digit,size_zuweisung)  \
    { var uintL _zero_counter; /* zählt die führenden Nullbits in digit     */\
      __asm__("cntlz %0,%1" : "=r" (_zero_counter) : "r" ((uint32)(digit)) ); \
      size_zuweisung (32-_zero_counter);                                      \
    }
 #else
  // new assembler syntax
  #define integerlength32(digit,size_zuweisung)  \
    { var uintL _zero_counter; /* zählt die führenden Nullbits in digit      */\
      __asm__("cntlzw %0,%1" : "=r" (_zero_counter) : "r" ((uint32)(digit)) ); \
      size_zuweisung (32-_zero_counter);                                       \
    }
 #endif
#elif defined(__GNUC__) && defined(__m88k__) && !defined(NO_ASM)
  #define integerlength32(digit,size_zuweisung)  \
    { var uintL _one_position; /* Position der führenden 1                */\
      __asm__("ff1 %0,%1" : "=r" (_one_position) : "r" ((uint32)(digit)) ); \
      size_zuweisung (1+_one_position);                                     \
    }
#elif defined(__GNUC__) && defined(__ibm032__) && !defined(NO_ASM) // RT/ROMP
  #define integerlength32(digit,size_zuweisung)  \
    { var uintL _x32 = (uint32)(digit);				\
      if (_x32 >= bit(16))					\
        { integerlength16(_x32>>16,size_zuweisung 16 + ); }	\
        else							\
        { integerlength16(_x32,size_zuweisung); }		\
    }
#else
  #define integerlength32(digit,size_zuweisung)  \
    { var uintC _bitsize = 1;						\
      var uintL _x32 = (uint32)(digit);					\
      /* _x32 hat höchstens 32 Bits.                                   */\
      if (_x32 >= bit(16)) { _x32 = _x32>>16; _bitsize += 16; }		\
      /* _x32 hat höchstens 16 Bits.                                   */\
      if (_x32 >= bit(8)) { _x32 = _x32>>8; _bitsize += 8; }		\
      /* _x32 hat höchstens 8 Bits.                                    */\
      if (_x32 >= bit(4)) { _x32 = _x32>>4; _bitsize += 4; }		\
      /* _x32 hat höchstens 4 Bits.                                    */\
      if (_x32 >= bit(2)) { _x32 = _x32>>2; _bitsize += 2; }		\
      /* _x32 hat höchstens 2 Bits.                                    */\
      if (_x32 >= bit(1)) { /* _x32 = _x32>>1; */ _bitsize += 1; }	\
      /* _x32 hat höchstens 1 Bit. Dieses Bit muß gesetzt sein.        */\
      size_zuweisung _bitsize;						\
    }
  #define GENERIC_INTEGERLENGTH32
#endif

// Bits einer 64-Bit-Zahl zählen:
// integerlength64(digit,size=);
// setzt size auf die höchste in digit vorkommende Bitnummer.
// > digit: ein uint64 >0
// < size: >0, <=64, mit 2^(size-1) <= digit < 2^size
#ifdef GENERIC_INTEGERLENGTH32
  #define integerlength64(digit,size_zuweisung)  \
    { var uintC _bitsize = 1;						\
      var uint64 _x64 = (uint64)(digit);				\
      /* _x64 hat höchstens 64 Bits.                                   */\
      if (_x64 >= bit(32)) { _x64 = _x64>>32; _bitsize += 32; }		\
      /* _x64 hat höchstens 32 Bits.                                   */\
      if (_x64 >= bit(16)) { _x64 = _x64>>16; _bitsize += 16; }		\
      /* _x64 hat höchstens 16 Bits.                                   */\
      if (_x64 >= bit(8)) { _x64 = _x64>>8; _bitsize += 8; }		\
      /* _x64 hat höchstens 8 Bits.                                    */\
      if (_x64 >= bit(4)) { _x64 = _x64>>4; _bitsize += 4; }		\
      /* _x64 hat höchstens 4 Bits.                                    */\
      if (_x64 >= bit(2)) { _x64 = _x64>>2; _bitsize += 2; }		\
      /* _x64 hat höchstens 2 Bits.                                    */\
      if (_x64 >= bit(1)) { /* _x64 = _x64>>1; */ _bitsize += 1; }	\
      /* _x64 hat höchstens 1 Bit. Dieses Bit muß gesetzt sein.        */\
      size_zuweisung _bitsize;						\
    }
#else
  #define integerlength64(digit,size_zuweisung)  \
    { var uint64 _x64 = (digit);                                        \
      var uintC _bitsize64 = 0;                                         \
      var uint32 _x32_from_integerlength64;                             \
      if (_x64 >= (1ULL << 32)) {                                       \
        _x32_from_integerlength64 = _x64>>32; _bitsize64 += 32;         \
      } else {                                                          \
        _x32_from_integerlength64 = _x64;                               \
      }                                                                 \
      integerlength32(_x32_from_integerlength64, size_zuweisung _bitsize64 + ); \
    }
#endif

// Bits einer uintC-Zahl zählen:
// integerlengthC(digit,size=);
// setzt size auf die höchste in digit vorkommende Bitnummer.
// > digit: ein uintC >0
// < size: >0, <=intCsize, mit 2^(size-1) <= digit < 2^size
  #if (intCsize==32)
    #define integerlengthC  integerlength32
  #endif
  #if (intCsize==64)
    #define integerlengthC  integerlength64
  #endif

// Hintere Nullbits eines 32-Bit-Wortes zählen:
// ord2_32(digit,count=);
// setzt size auf die kleinste in digit vorkommende Bitnummer.
// > digit: ein uint32 >0
// < count: >=0, <32, mit 2^count | digit, digit/2^count ungerade
  #if defined(__GNUC__) && defined(__i386__) && !defined(NO_ASM)
    #define ord2_32(digit,count_zuweisung)  \
      { var uintL _one_position; /* Position der letzten 1                    */\
        __asm__("bsfl %1,%0" : "=r" (_one_position) : "rm" ((uint32)(digit)) ); \
        count_zuweisung _one_position;                                          \
      }
    #define FAST_ORD2
  #elif defined(__sparc__) && !defined(__sparc64__)
    #define ord2_32(digit,count_zuweisung)  \
    { var uint32 n = (digit);                                             \
      n = n | -n;                                                         \
      n = (n<<4) + n;                                                     \
      n = (n<<6) + n;                                                     \
      n = n - (n<<16); /* or  n = n ^ (n<<16);  or  n = n &~ (n<<16);  */ \
      /* static const char ord2_tab [64] = {-1,0,1,12,2,6,-1,13,3,-1,7,-1,-1,-1,-1,14,10,4,-1,-1,8,-1,-1,25,-1,-1,-1,-1,-1,21,27,15,31,11,5,-1,-1,-1,-1,-1,9,-1,-1,24,-1,-1,20,26,30,-1,-1,-1,-1,23,-1,19,29,-1,22,18,28,17,16,-1}; */ \
      /* count_zuweisung ord2_tab[n>>26];                              */ \
      count_zuweisung "\377\000\001\014\002\006\377\015\003\377\007\377\377\377\377\016\012\004\377\377\010\377\377\031\377\377\377\377\377\025\033\017\037\013\005\377\377\377\377\377\011\377\377\030\377\377\024\032\036\377\377\377\377\027\377\023\035\377\026\022\034\021\020"[n>>26]; \
    }
    #define FAST_ORD2
  #else
    // Sei n = ord2(x). Dann ist logxor(x,x-1) = 2^n + (2^n-1) = 2^(n+1)-1.
    // Also  (ord2 x) = (1- (integer-length (logxor x (1- x)))) .
    #define ord2_32(digit,count_zuweisung)  \
      { var uint32 _digit = (digit) ^ ((digit) - 1);	\
        integerlength32(_digit,count_zuweisung -1 + )	\
      }
  #endif

// Hintere Nullbits eines 64-Bit-Wortes zählen:
// ord2_64(digit,count=);
// setzt size auf die kleinste in digit vorkommende Bitnummer.
// > digit: ein uint64 >0
// < count: >=0, <64, mit 2^count | digit, digit/2^count ungerade
  // Sei n = ord2(x). Dann ist logxor(x,x-1) = 2^n + (2^n-1) = 2^(n+1)-1.
  // Also  (ord2 x) = (1- (integer-length (logxor x (1- x)))) .
  #define ord2_64(digit,count_zuweisung)  \
    { var uint64 _digit = (digit) ^ ((digit) - 1);	\
      integerlength64(_digit,count_zuweisung -1 + )	\
    }


// Bits eines Wortes zählen.
// logcount_NN();
// > xNN: ein uintNN
// < xNN: Anzahl der darin gesetzten Bits
  // Bits von x8 zählen: (Input x8, Output x8)
  #define logcount_8()  \
    ( /* x8 besteht aus 8 1-Bit-Zählern (0,1).        */\
      x8 = (x8 & 0x55U) + ((x8 & 0xAAU) >> 1),		\
      /* x8 besteht aus 4 2-Bit-Zählern (0,1,2).      */\
      x8 = (x8 & 0x33U) + ((x8 & 0xCCU) >> 2),		\
      /* x8 besteht aus 2 4-Bit-Zählern (0,1,2,3,4).  */\
      x8 = (x8 & 0x0FU) + (x8 >> 4)			\
      /* x8 besteht aus 1 8-Bit-Zähler (0,...,8).     */\
    )
  // Bits von x16 zählen: (Input x16, Output x16)
  #define logcount_16()  \
    ( /* x16 besteht aus 16 1-Bit-Zählern (0,1).      */\
      x16 = (x16 & 0x5555U) + ((x16 & 0xAAAAU) >> 1),	\
      /* x16 besteht aus 8 2-Bit-Zählern (0,1,2).     */\
      x16 = (x16 & 0x3333U) + ((x16 & 0xCCCCU) >> 2),	\
      /* x16 besteht aus 4 4-Bit-Zählern (0,1,2,3,4). */\
      x16 = (x16 & 0x0F0FU) + ((x16 & 0xF0F0U) >> 4),	\
      /* x16 besteht aus 2 8-Bit-Zählern (0,...,8).   */\
      x16 = (x16 & 0x00FFU) + (x16 >> 8)		\
      /* x16 besteht aus 1 16-Bit-Zähler (0,...,16).  */\
    )
  // Bits von x32 zählen: (Input x32, Output x32)
  #define logcount_32()  \
    ( /* x32 besteht aus 32 1-Bit-Zählern (0,1).              */\
      x32 = (x32 & 0x55555555UL) + ((x32 & 0xAAAAAAAAUL) >> 1),	\
      /* x32 besteht aus 16 2-Bit-Zählern (0,1,2).            */\
      x32 = (x32 & 0x33333333UL) + ((x32 & 0xCCCCCCCCUL) >> 2),	\
      /* x32 besteht aus 8 4-Bit-Zählern (0,1,2,3,4).         */\
      x32 = high16(x32)+low16(x32),				\
      /* x32 besteht aus 4 4-Bit-Zählern (0,...,8).           */\
      x32 = (x32 & 0x0F0FU) + ((x32 & 0xF0F0U) >> 4),		\
      /* x32 besteht aus 2 8-Bit-Zählern (0,...,16).          */\
      x32 = (x32 & 0x00FFU) + (x32 >> 8)			\
      /* x32 besteht aus 1 16-Bit-Zähler (0,...,32).          */\
    )
  // Bits von x64 zählen: (Input x64, Output x64)
  #define logcount_64()  \
    ( /* x64 besteht aus 64 1-Bit-Zählern (0,1).                             */\
      x64 = (x64 & 0x5555555555555555ULL) + ((x64 & 0xAAAAAAAAAAAAAAAAULL) >> 1),\
      /* x64 besteht aus 32 2-Bit-Zählern (0,1,2).                           */\
      x64 = (x64 & 0x3333333333333333ULL) + ((x64 & 0xCCCCCCCCCCCCCCCCULL) >> 2),\
      /* x64 besteht aus 16 4-Bit-Zählern (0,1,2,3,4).                       */\
      x64 = (uint32)(x64 + (x64 >> 32)),				       \
      /* x64 besteht aus 8 4-Bit-Zählern (0,...,8).                          */\
      x64 = (x64 & 0x0F0F0F0FUL) + ((x64 & 0xF0F0F0F0UL) >> 4),		       \
      /* x64 besteht aus 4 8-Bit-Zählern (0,...,16).                         */\
      x64 = (x64 & 0x00FF00FFU) + ((x64 & 0xFF00FF00U) >> 8),		       \
      /* x64 besteht aus 2 16-Bit-Zählern (0,...,32).                        */\
      x64 = (x64 & 0x0000FFFFU) + (x64 >> 16)				       \
      /* x64 besteht aus 1 16-Bit-Zähler (0,...,64).                         */\
    )

}  // namespace cln

#endif /* _CL_LOW_H */
