// cl_FF internals

#ifndef _CL_FF_H
#define _CL_FF_H

#include "cln/number.h"
#include "cln/malloc.h"
#include "base/cl_low.h"
#include "float/cl_F.h"

#ifdef FAST_FLOAT
#include "base/cl_N.h"
#include "float/cl_F.h"
#endif

namespace cln {

typedef uint32 ffloat; // 32-bit float in IEEE format

union ffloatjanus {
	ffloat eksplicit;	// explicit value
	#ifdef FAST_FLOAT
	float machine_float;	// value as a C `float'
	#endif
};

#if defined(CL_WIDE_POINTERS)
#define FF_value_shift  32
inline ffloat cl_ffloat_value (const cl_FF& x)
{
	return x.word >> FF_value_shift;
}
#else
struct cl_heap_ffloat : cl_heap {
	ffloatjanus representation;
};
inline cl_heap_ffloat* TheFfloat (const cl_number& obj)
	{ return (cl_heap_ffloat*)(obj.pointer); }
inline ffloat cl_ffloat_value (const cl_FF& x)
{
	return TheFfloat(x)->representation.eksplicit;
}
#endif

// The word contains:
//   |..|.......|..........................|
//  sign exponent             mantissa

  #define FF_exp_len    8	// number of bits in the exponent
  #define FF_mant_len  23	// number of bits in the mantissa
				// (excluding the hidden bit)
  #define FF_exp_low   1		// minimum exponent
  #define FF_exp_mid   126		// exponent bias
  #define FF_exp_high  254		// maximum exponent, 255 is NaN/Inf
  #define FF_exp_shift  (FF_mant_len+FF_mant_shift) // lowest exponent bit
  #define FF_mant_shift  0			    // lowest mantissa bit
  #define FF_sign_shift  (32 - 1)	// = (FF_exp_len+FF_mant_len)

// Private constructor.
#if !defined(CL_WIDE_POINTERS)
inline cl_FF::cl_FF (cl_heap_ffloat* ptr) : cl_F ((cl_private_thing) ptr) {}
#endif

extern cl_class cl_class_ffloat;

// Builds a float from the explicit word.
#if defined(CL_WIDE_POINTERS)
inline cl_FF::cl_FF (struct cl_heap_ffloat * null, cl_uint w)
	: cl_F ((cl_private_thing) w) { unused null; }
inline const cl_FF allocate_ffloat (ffloat eksplicit)
{
	return cl_FF((struct cl_heap_ffloat *) 0, ((cl_uint)eksplicit << FF_value_shift) | (cl_FF_tag << cl_tag_shift));
}
#else
inline cl_heap_ffloat* allocate_ffloat (ffloat eksplicit)
{
	cl_heap_ffloat* p = (cl_heap_ffloat*) malloc_hook(sizeof(cl_heap_ffloat));
	p->refcount = 1;
	p->type = &cl_class_ffloat;
	p->representation.eksplicit = eksplicit;
	return p;
}
#endif

// Builds a float word from sign (0 or -1), exponent and mantissa.
inline uint32 make_FF_word (cl_sint sign, unsigned int exp, cl_uint mant)
{
	return (sign << FF_sign_shift)
	       | (exp << FF_exp_shift)
	       | ((mant & (bit(FF_mant_len)-1)) << FF_mant_shift);
}

// Builds a float from sign (0 or -1), exponent and mantissa.
inline const cl_FF make_FF (cl_sint sign, unsigned int exp, cl_uint mant)
{
	return allocate_ffloat(make_FF_word(sign,exp,mant));
}

#if defined(CL_WIDE_POINTERS)
// Single Float 0.0
  #define cl_FF_0  make_FF(0,0,0)
// Single Float 1.0
  #define cl_FF_1  make_FF(0,FF_exp_mid+1,bit(FF_mant_len))
// Single Float -1.0
  #define cl_FF_minus1  make_FF(-1,FF_exp_mid+1,bit(FF_mant_len))
#else
// Single Float 0.0
  extern const cl_FF cl_FF_0;
// Single Float 1.0
  extern const cl_FF cl_FF_1;
// Single Float -1.0
  extern const cl_FF cl_FF_minus1;
#endif


// Entpacken eines Single-Float:
// FF_decode(obj, zero_statement, sign=,exp=,mant=);
// zerlegt ein Single-Float obj.
// Ist obj=0.0, wird zero_statement ausgeführt.
// Sonst: cl_signean sign = Vorzeichen (0 = +, -1 = -),
//        sintL exp = Exponent (vorzeichenbehaftet),
//        uintL mant = Mantisse (>= 2^FF_mant_len, < 2^(FF_mant_len+1))
#define FF_uexp(x)  (((x) >> FF_mant_len) & (bit(FF_exp_len)-1))
#define FF_decode(obj, zero_statement, sign_zuweisung,exp_zuweisung,mant_zuweisung)  \
  { var ffloat _x = cl_ffloat_value(obj);				\
    var uintL uexp = FF_uexp(_x);					\
    if (uexp==0)							\
      { zero_statement } /* e=0 -> Zahl 0.0 */				\
      else								\
      { exp_zuweisung (sintL)(uexp - FF_exp_mid); /* Exponent */	\
        unused (sign_zuweisung sign_of((sint32)(_x))); /* Vorzeichen */\
        mant_zuweisung (bit(FF_mant_len) | (_x & (bit(FF_mant_len)-1))); \
  }   }

// Einpacken eines Single-Float:
// encode_FF(sign,exp,mant);
// liefert ein Single-Float.
// > cl_signean sign: Vorzeichen, 0 für +, -1 für negativ.
// > sintE exp: Exponent
// > uintL mant: Mantisse, sollte >= 2^FF_mant_len und < 2^(FF_mant_len+1) sein.
// < object ergebnis: ein Single-Float
// Der Exponent wird auf Überlauf/Unterlauf getestet.
inline const cl_FF encode_FF (cl_signean sign, sintE exp, uintL mant)
{
	if (exp < (sintE)(FF_exp_low-FF_exp_mid))
		{ if (underflow_allowed())
			{ throw floating_point_underflow_exception(); }
			else
			{ return cl_FF_0; }
		}
	else
	if (exp > (sintE)(FF_exp_high-FF_exp_mid))
		{ throw floating_point_overflow_exception(); }
	else
	return make_FF(sign, exp+FF_exp_mid, mant & (bit(FF_mant_len)-1));
}

#ifdef FAST_FLOAT
// Auspacken eines Floats:
inline float FF_to_float (const cl_FF& obj)
{
  #if defined(CL_WIDE_POINTERS) // eines der beiden 32-Bit-Wörter
    #if defined(__GNUC__)
      return ((ffloatjanus) { eksplicit: cl_ffloat_value(obj) }).machine_float;
    #else
      return *(float*)(&((uint32*)&(obj))[BIG_ENDIAN_P+(1-2*BIG_ENDIAN_P)*(FF_value_shift/32)]);
    #endif
  #else
    return TheFfloat(obj)->representation.machine_float;
  #endif
}
// Überprüfen und Einpacken eines von den 'float'-Routinen gelieferten
// IEEE-Floats.
// Klassifikation:
//   1 <= e <= 254 : normalisierte Zahl
//   e=0, m/=0: subnormale Zahl
//   e=0, m=0: vorzeichenbehaftete 0.0
//   e=255, m=0: vorzeichenbehaftete Infinity
//   e=255, m/=0: NaN
// Angabe der möglicherweise auftretenden Sonderfälle:
//   maybe_overflow: Operation läuft über, liefert IEEE-Infinity
//   maybe_subnormal: Ergebnis sehr klein, liefert IEEE-subnormale Zahl
//   maybe_underflow: Ergebnis sehr klein und /=0, liefert IEEE-Null
//   maybe_divide_0: Ergebnis unbestimmt, liefert IEEE-Infinity
//   maybe_nan: Ergebnis unbestimmt, liefert IEEE-NaN
  #define float_to_FF(expr,ergebnis_zuweisung,maybe_overflow,maybe_subnormal,maybe_underflow,maybe_divide_0,maybe_nan)  \
    { var ffloatjanus _erg; _erg.machine_float = (expr);		\
      if ((_erg.eksplicit & ((uint32)bit(FF_exp_len+FF_mant_len)-bit(FF_mant_len))) == 0) /* e=0 ? */\
        { if ((maybe_underflow						\
               || (maybe_subnormal && !((_erg.eksplicit << 1) == 0))	\
              )								\
              && underflow_allowed()					\
             )								\
            { throw floating_point_underflow_exception(); } /* subnormal oder noch kleiner -> Underflow */\
            else							\
            { ergebnis_zuweisung cl_FF_0; } /* +/- 0.0 -> 0.0 */	\
        }								\
      elif ((maybe_overflow || maybe_divide_0)				\
            && (((~_erg.eksplicit) & ((uint32)bit(FF_exp_len+FF_mant_len)-bit(FF_mant_len))) == 0) /* e=255 ? */\
           )								\
        { if (maybe_nan && !((_erg.eksplicit << (32-FF_mant_len)) == 0)) \
            { throw division_by_0_exception(); } /* NaN, also Singularität -> "Division durch 0" */\
          else /* Infinity */						\
          if (!maybe_overflow || maybe_divide_0)			\
            { throw division_by_0_exception(); } /* Infinity, Division durch 0 */\
            else							\
            { throw floating_point_overflow_exception(); } /* Infinity, Overflow */\
        }								\
      else								\
        { ergebnis_zuweisung allocate_ffloat(_erg.eksplicit); }		\
    }
#endif

// Liefert zu einem Single-Float x : (futruncate x), ein FF.
// x wird von der 0 weg zur nächsten ganzen Zahl gerundet.
extern const cl_FF futruncate (const cl_FF& x);

// FF_to_I(x) wandelt ein Single-Float x, das eine ganze Zahl darstellt,
// in ein Integer um.
extern const cl_I cl_FF_to_I (const cl_FF& x);

// cl_I_to_FF(x) wandelt ein Integer x in ein Single-Float um und rundet dabei.
extern const cl_FF cl_I_to_FF (const cl_I& x);

// cl_RA_to_FF(x) wandelt eine rationale Zahl x in ein Single-Float um
// und rundet dabei.
extern const cl_FF cl_RA_to_FF (const cl_RA& x);


// IEEE-Single-Float:
// Bit 31 = s, Bits 30..23 = e, Bits 22..0 = m.
//   e=0, m=0: vorzeichenbehaftete 0.0
//   e=0, m/=0: subnormale Zahl,
//     Wert = (-1)^s * 2^(1-126) * [ 0 . 0 m22 ... m0 ]
//   1 <= e <= 254 : normalisierte Zahl,
//     Wert = (-1)^s * 2^(e-126) * [ 0 . 1 m22 ... m0 ]
//   e=255, m=0: vorzeichenbehaftete Infinity
//   e=255, m/=0: NaN

// cl_float_to_FF_pointer(val) wandelt ein IEEE-Single-Float val in ein Single-Float um.
extern cl_private_thing cl_float_to_FF_pointer (const float val);

// cl_FF_to_float(obj,&val);
// wandelt ein Single-Float obj in ein IEEE-Single-Float val um.
extern void cl_FF_to_float (const cl_FF& obj, ffloatjanus* val_);

}  // namespace cln

#endif /* _CL_FF_H */
