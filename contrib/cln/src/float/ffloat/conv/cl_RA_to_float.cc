// float_approx().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "rational/cl_RA.h"
#include "cln/integer.h"
#include "integer/cl_I.h"
#include "float/cl_F.h"

namespace cln {

float float_approx (const cl_RA& x)
{
// Method: same as cl_RA_to_FF().
      if (integerp(x)) {
        DeclareType(cl_I,x);
        return float_approx(x);
      }
 {    // x Ratio
      DeclareType(cl_RT,x);
      union { ffloat eksplicit; float machine_float; } u;
      var cl_I a = numerator(x); // +/- a
      var const cl_I& b = denominator(x); // b
      var cl_signean sign = -(cl_signean)minusp(a); // Vorzeichen
      if (!(sign==0)) { a = -a; } // Betrag nehmen, liefert a
      var sintC lendiff = (sintC)integer_length(a) // (integer-length a)
                          - (sintC)integer_length(b); // (integer-length b)
      if (lendiff > FF_exp_high-FF_exp_mid) // Exponent >= n-m > Obergrenze ?
        { u.eksplicit = make_FF_word(sign,bit(FF_exp_len)-1,0); // Infinity
          return u.machine_float;
        }
      if (lendiff < FF_exp_low-FF_exp_mid-2) // Exponent <= n-m+2 < Untergrenze ?
        { u.eksplicit = make_FF_word(sign,0,0); // 0.0
          return u.machine_float;
        }
      var cl_I zaehler;
      var cl_I nenner;
      if (lendiff >= FF_mant_len+2)
        // n-m-25>=0
        { nenner = ash(b,lendiff - (FF_mant_len+2)); // (ash b n-m-25)
          zaehler = a; // a
        }
        else
        { zaehler = ash(a,(FF_mant_len+2) - lendiff); // (ash a -n+m+25)
          nenner = b; // b
        }
      // Division zaehler/nenner durchführen:
      var cl_I_div_t q_r = cl_divide(zaehler,nenner);
      var cl_I& q = q_r.quotient;
      var cl_I& r = q_r.remainder;
      // 2^24 <= q < 2^26, also ist q Fixnum oder Bignum mit bn_minlength Digits.
      var uint32 mant = ((FF_mant_len+3 < cl_value_len)
                          ? FN_to_UV(q)
                          : cl_I_to_UL(q)
                        );
      if (mant >= bit(FF_mant_len+2))
        // 2^25 <= q < 2^26, schiebe um 2 Bits nach rechts
        { var uintL rounding_bits = mant & (bit(2)-1);
          lendiff = lendiff+1; // Exponent := n-m+1
          mant = mant >> 2;
          if ( (rounding_bits < bit(1)) // 00,01 werden abgerundet
               || ( (rounding_bits == bit(1)) // 10
                    && (eq(r,0)) // und genau halbzahlig (r=0)
                    && ((mant & bit(0)) ==0) // -> round-to-even
             )    )
            // abrunden
            goto ab;
            else
            // aufrunden
            goto auf;
        }
        else
        { var uintL rounding_bit = mant & bit(0);
          mant = mant >> 1;
          if ( (rounding_bit == 0) // 0 wird abgerundet
               || ( (eq(r,0)) // genau halbzahlig (r=0)
                    && ((mant & bit(0)) ==0) // -> round-to-even
             )    )
            // abrunden
            goto ab;
            else
            // aufrunden
            goto auf;
        }
      auf:
      mant += 1;
      if (mant >= bit(FF_mant_len+1)) // rounding overflow?
        { mant = mant>>1; lendiff = lendiff+1; }
      ab:
      // Fertig.
      if (lendiff < (sintL)(FF_exp_low-FF_exp_mid))
        { u.eksplicit = make_FF_word(sign,0,0); }
      else if (lendiff > (sintL)(FF_exp_high-FF_exp_mid))
        { u.eksplicit = make_FF_word(sign,bit(FF_exp_len)-1,0); } // Infinity
      else
        { u.eksplicit = make_FF_word(sign,lendiff+FF_exp_mid,mant); }
      return u.machine_float;
}}

}  // namespace cln
