// cl_DF internals

#ifndef _CL_DF_H
#define _CL_DF_H

#include "cln/number.h"
#include "cln/malloc.h"
#include "base/cl_low.h"
#include "float/cl_F.h"

#ifdef FAST_DOUBLE
#include "base/cl_N.h"
#include "float/cl_F.h"
#endif

namespace cln {

typedef // 64-bit float in IEEE format
	#if (cl_word_size==64)
	  // Sign/Exponent/Mantissa
	  uint64
	#else
	  // Sign/Exponent/MantissaHigh and MantissaLow
	  #if defined(double_wordorder_bigendian_p)
	    #if double_wordorder_bigendian_p
	      struct { uint32 semhi, mlo; }
	    #else
	      struct { uint32 mlo, semhi; }
	    #endif
	  #else
	    #if CL_CPU_BIG_ENDIAN_P
	      struct { uint32 semhi, mlo; }
	    #else
	      struct { uint32 mlo, semhi; }
	    #endif
	  #endif
	#endif
	dfloat;

union dfloatjanus {
	dfloat eksplicit;	// explicit value
	#ifdef FAST_DOUBLE
	double machine_double;	// value as a C `double'
	#endif
};
struct cl_heap_dfloat : cl_heap {
	dfloatjanus representation;
};
inline cl_heap_dfloat* TheDfloat (const cl_number& obj)
	{ return (cl_heap_dfloat*)(obj.pointer); }
inline dfloat& cl_dfloat_value (/* const ?? */ cl_DF& x)
{
	return TheDfloat(x)->representation.eksplicit;
}

// The double-word contains:
//   |..|.......|........................................|
//  sign exponent                  mantissa

  #define DF_exp_len   11	// number of bits in the exponent
  #define DF_mant_len  52	// number of bits in the mantissa
				// (excluding the hidden bit)
  #define DF_exp_low   1		// minimum exponent
  #define DF_exp_mid   1022		// exponent bias
  #define DF_exp_high  2046		// maximum exponent, 2047 is NaN/Inf
  #define DF_exp_shift  (DF_mant_len+DF_mant_shift) // lowest exponent bit
  #define DF_mant_shift  0			    // lowest mantissa bit
  #define DF_sign_shift  (64 - 1)	// = (DF_exp_len+DF_mant_len)

// Private constructor.
inline cl_DF::cl_DF (cl_heap_dfloat* ptr) : cl_F ((cl_private_thing) ptr) {}

extern cl_class cl_class_dfloat;

// Builds a float from the explicit words.
#if (cl_word_size==64)
inline cl_heap_dfloat* allocate_dfloat (dfloat eksplicit)
{
	cl_heap_dfloat* p = (cl_heap_dfloat*) malloc_hook(sizeof(cl_heap_dfloat));
	p->refcount = 1;
	p->type = &cl_class_dfloat;
	p->representation.eksplicit = eksplicit;
	return p;
}
#else
inline cl_heap_dfloat* allocate_dfloat (uint32 semhi, uint32 mlo)
{
	cl_heap_dfloat* p = (cl_heap_dfloat*) malloc_hook(sizeof(cl_heap_dfloat));
	p->refcount = 1;
	p->type = &cl_class_dfloat;
	p->representation.eksplicit.semhi = semhi;
	p->representation.eksplicit.mlo   = mlo;
	return p;
}
#endif

// Double Float 0.0
  extern const cl_DF cl_DF_0;
// Double Float 1.0
  extern const cl_DF cl_DF_1;
// Double Float -1.0
  extern const cl_DF cl_DF_minus1;


#define dfloat_value  representation.eksplicit

// Entpacken eines Double-Float:
#if (cl_word_size==64)
// DF_decode(obj, zero_statement, sign=,exp=,mant=);
// zerlegt ein Double-Float obj.
// Ist obj=0.0, wird zero_statement ausgeführt.
// Sonst: cl_signean sign = Vorzeichen (0 = +, -1 = -),
//        sintL exp = Exponent (vorzeichenbehaftet),
//        uintQ mant = Mantisse (>= 2^DF_mant_len, < 2^(DF_mant_len+1))
  #define dfloat_value_semhi  dfloat_value
  #define DF_uexp(x)  (((x) >> DF_mant_len) & (bit(DF_exp_len)-1))
  #define DF_decode(obj, zero_statement, sign_zuweisung,exp_zuweisung,mant_zuweisung)  \
    { var dfloat _x = TheDfloat(obj)->dfloat_value;			\
      var uintL uexp = DF_uexp(_x);					\
      if (uexp==0)							\
        { zero_statement } /* e=0 -> Zahl 0.0 */			\
        else								\
        { exp_zuweisung (sintL)(uexp - DF_exp_mid); /* Exponent */	\
          unused (sign_zuweisung ((sint64)_x >> 63)); /* Vorzeichen */	\
          mant_zuweisung (bit(DF_mant_len) | (_x & (bit(DF_mant_len)-1))); \
    }   }
#else
// DF_decode2(obj, zero_statement, sign=,exp=,manthi=,mantlo=);
// zerlegt ein Double-Float obj.
// Ist obj=0.0, wird zero_statement ausgeführt.
// Sonst: cl_signean sign = Vorzeichen (0 = +, -1 = -),
//        sintL exp = Exponent (vorzeichenbehaftet),
//        uintL manthi,mantlo = Mantisse 2^32*manthi+mantlo
//                              (>= 2^DF_mant_len, < 2^(DF_mant_len+1))
  #define dfloat_value_semhi  dfloat_value.semhi
  #define DF_uexp(semhi)  (((semhi) >> (DF_mant_len-32)) & (bit(DF_exp_len)-1))
  #define DF_decode2(obj, zero_statement, sign_zuweisung,exp_zuweisung,manthi_zuweisung,mantlo_zuweisung)  \
    { var uint32 semhi = TheDfloat(obj)->dfloat_value.semhi;		\
      var uint32 mlo = TheDfloat(obj)->dfloat_value.mlo;		\
      var uintL uexp = DF_uexp(semhi);					\
      if (uexp==0)							\
        { zero_statement } /* e=0 -> Zahl 0.0 */			\
        else								\
        { exp_zuweisung (sintL)(uexp - DF_exp_mid); /* Exponent */	\
          unused (sign_zuweisung sign_of((sint32)(semhi))); /* Vorzeichen */\
          manthi_zuweisung (bit(DF_mant_len-32) | (semhi & (bit(DF_mant_len-32)-1))); \
          mantlo_zuweisung mlo;						\
    }   }
#endif

// Einpacken eines Double-Float:
#if (cl_word_size==64)
// encode_DF(sign,exp,mant)
// liefert ein Double-Float.
// > cl_signean sign: Vorzeichen, 0 für +, -1 für negativ.
// > sintE exp: Exponent
// > uintQ mant: Mantisse, sollte >= 2^DF_mant_len und < 2^(DF_mant_len+1) sein.
// < cl_DF ergebnis: ein Double-Float
// Der Exponent wird auf Überlauf/Unterlauf getestet.
inline const cl_DF encode_DF (cl_signean sign, sintE exp, uintQ mant)
{
      if (exp < (sintE)(DF_exp_low-DF_exp_mid))
        { if (underflow_allowed())
            { throw floating_point_underflow_exception(); }
            else
            { return cl_DF_0; }
        }
      else
      if (exp > (sintE)(DF_exp_high-DF_exp_mid))
        { throw floating_point_overflow_exception(); }
      else
      return allocate_dfloat
        (  ((sint64)sign & bit(63))                  /* Vorzeichen */
         | ((uint64)(exp+DF_exp_mid) << DF_mant_len) /* Exponent   */
         | ((uint64)mant & (bit(DF_mant_len)-1))     /* Mantisse   */
        );
}
#else
// encode_DF(sign,exp,manthi,mantlo)
// liefert ein Double-Float.
// > cl_signean sign: Vorzeichen, 0 für +, -1 für negativ.
// > sintE exp: Exponent
// > uintL manthi,mantlo: Mantisse 2^32*manthi+mantlo,
//                        sollte >= 2^DF_mant_len und < 2^(DF_mant_len+1) sein.
// < cl_DF ergebnis: ein Double-Float
// Der Exponent wird auf Überlauf/Unterlauf getestet.
inline const cl_DF encode_DF (cl_signean sign, sintE exp, uintL manthi, uintL mantlo)
{
      if (exp < (sintE)(DF_exp_low-DF_exp_mid))
        { if (underflow_allowed())
            { throw floating_point_underflow_exception(); }
            else
            { return cl_DF_0; }
        }
      else
      if (exp > (sintE)(DF_exp_high-DF_exp_mid))
        { throw floating_point_overflow_exception(); }
      else
      return allocate_dfloat
        (  ((sint32)sign & bit(31))                       /* Vorzeichen */
         | ((uint32)(exp+DF_exp_mid) << (DF_mant_len-32)) /* Exponent   */
         | ((uint32)manthi & (bit(DF_mant_len-32)-1))     /* Mantisse   */
         , mantlo
        );
}
#endif

#ifdef FAST_DOUBLE
// Auspacken eines Double:
inline double DF_to_double (const cl_DF& obj)
{
	return TheDfloat(obj)->representation.machine_double;
}
// Überprüfen und Einpacken eines von den 'double'-Routinen gelieferten
// IEEE-Floats.
// Klassifikation:
//   1 <= e <= 2046 : normalisierte Zahl
//   e=0, m/=0: subnormale Zahl
//   e=0, m=0: vorzeichenbehaftete 0.0
//   e=2047, m=0: vorzeichenbehaftete Infinity
//   e=2047, m/=0: NaN
// Angabe der möglicherweise auftretenden Sonderfälle:
//   maybe_overflow: Operation läuft über, liefert IEEE-Infinity
//   maybe_subnormal: Ergebnis sehr klein, liefert IEEE-subnormale Zahl
//   maybe_underflow: Ergebnis sehr klein und /=0, liefert IEEE-Null
//   maybe_divide_0: Ergebnis unbestimmt, liefert IEEE-Infinity
//   maybe_nan: Ergebnis unbestimmt, liefert IEEE-NaN
#if (cl_word_size==64)
  #define double_to_DF(expr,ergebnis_zuweisung,maybe_overflow,maybe_subnormal,maybe_underflow,maybe_divide_0,maybe_nan)  \
    { var dfloatjanus _erg; _erg.machine_double = (expr);		\
      if ((_erg.eksplicit & ((uint64)bit(DF_exp_len+DF_mant_len)-bit(DF_mant_len))) == 0) /* e=0 ? */\
        { if ((maybe_underflow						\
               || (maybe_subnormal && !((_erg.eksplicit << 1) == 0))	\
              )								\
              && underflow_allowed()					\
             )								\
            { throw floating_point_underflow_exception(); } /* subnormal oder noch kleiner -> Underflow */\
            else							\
            { ergebnis_zuweisung cl_DF_0; } /* +/- 0.0 -> 0.0 */	\
        }								\
      elif ((maybe_overflow || maybe_divide_0)				\
            && (((~_erg.eksplicit) & ((uint64)bit(DF_exp_len+DF_mant_len)-bit(DF_mant_len))) == 0) /* e=2047 ? */\
           )								\
        { if (maybe_nan && !((_erg.eksplicit<<(64-DF_mant_len)) == 0))	\
            { throw division_by_0_exception(); } /* NaN, also Singularität -> "Division durch 0" */\
          else /* Infinity */						\
          if (!maybe_overflow || maybe_divide_0)			\
            { throw division_by_0_exception(); } /* Infinity, Division durch 0 */\
            else							\
            { throw floating_point_overflow_exception(); } /* Infinity, Overflow */\
        }								\
      else								\
        { ergebnis_zuweisung allocate_dfloat(_erg.eksplicit); }		\
    }
#else
  #define double_to_DF(expr,ergebnis_zuweisung,maybe_overflow,maybe_subnormal,maybe_underflow,maybe_divide_0,maybe_nan)  \
    { var dfloatjanus _erg; _erg.machine_double = (expr);                 \
      if ((_erg.eksplicit.semhi & ((uint32)bit(DF_exp_len+DF_mant_len-32)-bit(DF_mant_len-32))) == 0) /* e=0 ? */\
        { if ((maybe_underflow                                            \
               || (maybe_subnormal                                        \
                   && !(((_erg.eksplicit.semhi << 1) == 0) && (_erg.eksplicit.mlo == 0)) \
              )   )                                                       \
              && underflow_allowed()                                      \
             )                                                            \
            { throw floating_point_underflow_exception(); } /* subnormal oder noch kleiner -> Underflow */\
            else                                                          \
            { ergebnis_zuweisung cl_DF_0; } /* +/- 0.0 -> 0.0           */\
        }                                                                 \
      elif ((maybe_overflow || maybe_divide_0)                            \
            && (((~_erg.eksplicit.semhi) & ((uint32)bit(DF_exp_len+DF_mant_len-32)-bit(DF_mant_len-32))) == 0) /* e=2047 ?  */\
           )                                                              \
        { if (maybe_nan && !(((_erg.eksplicit.semhi<<(64-DF_mant_len)) == 0) && (_erg.eksplicit.mlo==0))) \
            { throw division_by_0_exception(); } /* NaN, also Singularität -> "Division durch 0"  */\
          else /* Infinity                                              */\
          if (!maybe_overflow || maybe_divide_0)                          \
            { throw division_by_0_exception(); } /* Infinity, Division durch 0 */\
            else                                                          \
            { throw floating_point_overflow_exception(); } /* Infinity, Overflow */\
        }                                                                 \
      else                                                                \
        { ergebnis_zuweisung allocate_dfloat(_erg.eksplicit.semhi,_erg.eksplicit.mlo); }  \
    }
#endif
#endif

// Liefert zu einem Double-Float x : (futruncate x), ein DF.
// x wird von der 0 weg zur nächsten ganzen Zahl gerundet.
extern const cl_DF futruncate (const cl_DF& x);

// DF_to_I(x) wandelt ein Double-Float x, das eine ganze Zahl darstellt,
// in ein Integer um.
extern const cl_I cl_DF_to_I (const cl_DF& x);

// cl_I_to_DF(x) wandelt ein Integer x in ein Double-Float um und rundet dabei.
extern const cl_DF cl_I_to_DF (const cl_I& x);

// cl_RA_to_DF(x) wandelt eine rationale Zahl x in ein Double-Float um
// und rundet dabei.
extern const cl_DF cl_RA_to_DF (const cl_RA& x);


// IEEE-Double-Float:
// Bit 63 = s, Bits 62..52 = e, Bits 51..0 = m.
//   e=0, m=0: vorzeichenbehaftete 0.0
//   e=0, m/=0: subnormale Zahl,
//     Wert = (-1)^s * 2^(1-1022) * [ 0 . 0 m51 ... m0 ]
//   1 <= e <= 2046 : normalisierte Zahl,
//     Wert = (-1)^s * 2^(e-1022) * [ 0 . 1 m51 ... m0 ]
//   e=2047, m=0: vorzeichenbehaftete Infinity
//   e=2047, m/=0: NaN

// cl_double_to_DF_pointer(val) wandelt ein IEEE-Double-Float val in ein Double-Float um.
extern cl_heap_dfloat* cl_double_to_DF_pointer (const double val);

// cl_DF_to_double(obj,&val);
// wandelt ein Double-Float obj in ein IEEE-Double-Float val um.
extern void cl_DF_to_double (const cl_DF& obj, dfloatjanus* val_);

}  // namespace cln

#endif /* _CL_DF_H */
