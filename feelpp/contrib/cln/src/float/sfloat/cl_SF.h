// cl_SF internals

#ifndef _CL_SF_H
#define _CL_SF_H

#include "cln/number.h"
#include "float/cl_F.h"

namespace cln {

// The immediate word contains:
//   |..|.......|..........................|....|
//  sign exponent             mantissa      tag

  #define SF_value_shift 7	// could also be = cl_value_shift
  #define SF_exp_len    8	// number of bits in the exponent
  #define SF_mant_len  16	// number of bits in the mantissa
				// (excluding the hidden bit)
  #define SF_mant_hiddenbit 1	// yes, we have a hidden bit representation
				// (this is hardwired in some of the code)
  #define SF_exp_low   1			// minimum exponent
  #define SF_exp_mid   bit(SF_exp_len-1)	// exponent bias
  #define SF_exp_high  (bit(SF_exp_len)-1)	// maximum exponent
  #define SF_exp_shift  (SF_mant_len+SF_mant_shift) // lowest exponent bit
  #define SF_mant_shift  SF_value_shift		    // lowest mantissa bit
  #define SF_sign_shift  (cl_pointer_size - 1)

// Builds a float from the immediate word.
inline cl_SF::cl_SF (struct cl_sfloat * null, cl_uint w)
	: cl_F ((cl_private_thing) w) { unused null; }
inline const cl_SF cl_SF_from_word (cl_uint word)
{
	return cl_SF((struct cl_sfloat *) 0, word);
}

// Builds a float word from sign (0 or -1), exponent and mantissa.
inline cl_uint make_SF_word (cl_sint sign, unsigned int exp, cl_uint mant)
{
	return (sign & ((cl_uint)1 << SF_sign_shift))
	       | (exp << SF_exp_shift)
	       #if SF_mant_hiddenbit
	       | ((mant & (bit(SF_mant_len)-1)) << SF_mant_shift)
	       #else
	       | (mant << SF_mant_shift)
	       #endif
	       | (cl_SF_tag << cl_tag_shift);
}

// Builds a float from sign (0 or -1), exponent and mantissa.
inline const cl_SF make_SF (cl_sint sign, unsigned int exp, cl_uint mant)
{
	return cl_SF_from_word(make_SF_word(sign,exp,mant));
}

// Short Float 0.0
  #define SF_0  make_SF(0,0,0)
// Short Float 1.0
  #define SF_1  make_SF(0,SF_exp_mid+1,bit(SF_mant_len))
// Short Float -1.0
  #define SF_minus1  make_SF(-1,SF_exp_mid+1,bit(SF_mant_len))


// Entpacken eines Short-Float:
// SF_decode(obj, zero_statement, sign=,exp=,mant=);
// zerlegt ein Short-Float obj.
// Ist obj=0.0, wird zero_statement ausgeführt.
// Sonst: cl_signean sign = Vorzeichen (0 = +, -1 = -),
//        sintL exp = Exponent (vorzeichenbehaftet),
//        uintL mant = Mantisse (>= 2^SF_mant_len, < 2^(SF_mant_len+1))
inline uintL SF_uexp (const cl_SF& x)
{
	return (x.word >> SF_exp_shift) & (bit(SF_exp_len)-1);
}
inline cl_signean SF_sign (const cl_SF& x)
{
	return ((cl_sint)x.word << (cl_pointer_size-1 - SF_sign_shift)) >> (cl_pointer_size-1);
}
inline uintL SF_mant (const cl_SF& x)
{
	return
	       #if SF_mant_hiddenbit
	       bit(SF_mant_len) |
	       #endif
	       ((uintL)(x.word >> SF_mant_shift) & (bit(SF_mant_len)-1));
}
#define SF_decode(_x, zero_statement, sign_zuweisung,exp_zuweisung,mant_zuweisung)  \
  { var uintL uexp = SF_uexp(_x);					\
    if (uexp==0)							\
      { zero_statement } /* e=0 -> Zahl 0.0 */				\
      else								\
      { exp_zuweisung (sintL)(uexp - SF_exp_mid);	/* Exponent */	\
        unused (sign_zuweisung SF_sign(_x));		/* Vorzeichen */\
        mant_zuweisung SF_mant(_x);			/* Mantisse */  \
  }   }

// Einpacken eines Short-Float:
// encode_SF(sign,exp,mant)
// liefert ein Short-Float.
// > cl_signean sign: Vorzeichen, 0 für +, -1 für negativ.
// > sintE exp: Exponent
// > uintL mant: Mantisse, sollte >= 2^SF_mant_len und < 2^(SF_mant_len+1) sein.
// < object ergebnis: ein Short-Float
// Der Exponent wird auf Überlauf/Unterlauf getestet.
inline const cl_SF encode_SF (cl_signean sign, sintE exp, uintL mant)
{
	if (exp < (sintE)(SF_exp_low-SF_exp_mid))
	  { if (underflow_allowed())
	      { throw floating_point_underflow_exception(); }
	      else
	      { return SF_0; }
	  }
	else
	if (exp > (sintE)(SF_exp_high-SF_exp_mid))
	  { throw floating_point_overflow_exception(); }
	else
	return make_SF(sign, exp+SF_exp_mid, mant);
}


// Liefert zu einem Short-Float x : (futruncate x), ein SF.
// x wird von der 0 weg zur nächsten ganzen Zahl gerundet.
extern const cl_SF futruncate (const cl_SF& x);

// SF_to_I(x) wandelt ein Short-Float x, das eine ganze Zahl darstellt,
// in ein Integer um.
extern const cl_I cl_SF_to_I (const cl_SF& x);

// cl_I_to_SF(x) wandelt ein Integer x in ein Short-Float um und rundet dabei.
extern const cl_SF cl_I_to_SF (const cl_I& x);

// cl_RA_to_SF(x) wandelt eine rationale Zahl x in ein Short-Float um
// und rundet dabei.
extern const cl_SF cl_RA_to_SF (const cl_RA& x);

}  // namespace cln

#endif /* _CL_SF_H */
