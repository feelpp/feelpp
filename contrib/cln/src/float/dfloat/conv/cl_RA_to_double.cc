// double_approx().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "rational/cl_RA.h"
#include "cln/integer.h"
#include "integer/cl_I.h"
#include "float/cl_F.h"

namespace cln {

double double_approx (const cl_RA& x)
{
// Method: same as cl_RA_to_DF().
      if (integerp(x)) {
        DeclareType(cl_I,x);
        return double_approx(x);
      }
 {    // x Ratio
      DeclareType(cl_RT,x);
      union { dfloat eksplicit; double machine_double; } u;
      var cl_I a = numerator(x); // +/- a
      var const cl_I& b = denominator(x); // b
      var cl_signean sign = -(cl_signean)minusp(a); // Vorzeichen
      if (!(sign==0)) { a = -a; } // Betrag nehmen, liefert a
      var sintC lendiff = (sintC)integer_length(a) // (integer-length a)
                          - (sintC)integer_length(b); // (integer-length b)
      if (lendiff > DF_exp_high-DF_exp_mid) // Exponent >= n-m > Obergrenze ?
        {
          #if (cl_word_size==64)
          u.eksplicit =
              ((sint64)sign & bit(63))
            | ((uint64)(bit(DF_exp_len)-1) << DF_mant_len); // Infinity
          #else
          u.eksplicit.semhi =
              ((sint32)sign & bit(31))
            | ((uint32)(bit(DF_exp_len)-1) << (DF_mant_len-32)); // Infinity
          u.eksplicit.mlo = 0;
          #endif
          return u.machine_double;
        }
      if (lendiff < DF_exp_low-DF_exp_mid-2) // Exponent <= n-m+2 < Untergrenze ?
        {
          #if (cl_word_size==64)
          u.eksplicit = ((sint64)sign & bit(63)); // 0.0
          #else
          u.eksplicit.semhi = ((sint32)sign & bit(31)); // 0.0
          u.eksplicit.mlo = 0;
          #endif
          return u.machine_double;
        }
      var cl_I zaehler;
      var cl_I nenner;
      if (lendiff >= DF_mant_len+2)
        // n-m-54>=0
        { nenner = ash(b,lendiff - (DF_mant_len+2)); // (ash b n-m-54)
          zaehler = a; // a
        }
        else
        { zaehler = ash(a,(DF_mant_len+2) - lendiff); // (ash a -n+m+54)
          nenner = b; // b
        }
      // Division zaehler/nenner durchführen:
      var cl_I_div_t q_r = cl_divide(zaehler,nenner);
      var cl_I& q = q_r.quotient;  // 2^53 <= q < 2^55
      var cl_I& r = q_r.remainder;
      #if (cl_word_size==64)
      # if (cl_value_len-1 > 55)
      var uint64 mant = FN_to_V(q);  // q is a fixnum!
      # elif (cl_value_len-1 <= 53)
      var uint64 mant = get_max64_Dptr(55,BN_MSDptr(q));  // q is a bignum!
      # else
      var uint64 mant = fixnump(q) ? FN_to_V(q) : get_max64_Dptr(55,BN_MSDptr(q));
      # endif
      if (mant >= bit(DF_mant_len+2))
        // 2^54 <= q < 2^55, schiebe um 2 Bits nach rechts
        { var uint64 rounding_bits = mant & (bit(2)-1);
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
        { var uint64 rounding_bit = mant & bit(0);
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
      if (mant >= bit(DF_mant_len+1)) // rounding overflow?
        { mant = mant>>1; lendiff = lendiff+1; }
      ab:
      // Fertig.
      if (lendiff < (sintL)(DF_exp_low-DF_exp_mid))
        { u.eksplicit = ((sint64)sign & bit(63)); } // 0.0
      else if (lendiff > (sintL)(DF_exp_high-DF_exp_mid))
        { u.eksplicit =
              ((sint64)sign & bit(63))
            | ((uint64)(bit(DF_exp_len)-1) << DF_mant_len); // Infinity
        }
      else
        { u.eksplicit =
              ((sint64)sign & bit(63))                      /* Vorzeichen */
            | ((uint64)(lendiff+DF_exp_mid) << DF_mant_len) /* Exponent   */
            | ((uint64)mant & (bit(DF_mant_len)-1));        /* Mantisse   */
        }
      return u.machine_double;
      #else
      // q is bignum with ceiling(55/intDsize) Digits.
      var const uintD* ptr = BN_MSDptr(q);
      var uint32 manthi = get_max32_Dptr(23,ptr);
      var uint32 mantlo = get_32_Dptr(ptr mspop ceiling(23,intDsize));
      if (manthi >= bit(DF_mant_len-32+2))
        // 2^54 <= q < 2^55, schiebe um 2 Bits nach rechts
        { var uintL rounding_bits = mantlo & (bit(2)-1);
          lendiff = lendiff+1; // Exponent := n-m+1
          mantlo = (mantlo >> 2) | (manthi << 30); manthi = manthi >> 2;
          if ( (rounding_bits < bit(1)) // 00,01 werden abgerundet
               || ( (rounding_bits == bit(1)) // 10
                    && (eq(r,0)) // und genau halbzahlig (r=0)
                    && ((mantlo & bit(0)) ==0) // -> round-to-even
             )    )
            // abrunden
            goto ab;
            else
            // aufrunden
            goto auf;
        }
        else
        { var uintL rounding_bit = mantlo & bit(0);
          mantlo = (mantlo >> 1) | (manthi << 31); manthi = manthi >> 1;
          if ( (rounding_bit == 0) // 0 wird abgerundet
               || ( (eq(r,0)) // genau halbzahlig (r=0)
                    && ((mantlo & bit(0)) ==0) // -> round-to-even
             )    )
            // abrunden
            goto ab;
            else
            // aufrunden
            goto auf;
        }
      auf:
      mantlo += 1;
      if (mantlo==0)
        { manthi += 1;
          if (manthi >= bit(DF_mant_len-32+1)) // rounding overflow?
            { manthi = manthi>>1; lendiff = lendiff+1; }
        }
      ab:
      // Fertig.
      if (lendiff < (sintL)(DF_exp_low-DF_exp_mid))
        { u.eksplicit.semhi = ((sint32)sign & bit(31)); // 0.0
          u.eksplicit.mlo = 0;
        }
      else if (lendiff > (sintL)(DF_exp_high-DF_exp_mid))
        { u.eksplicit.semhi =
              ((sint32)sign & bit(31))
            | ((uint32)(bit(DF_exp_len)-1) << (DF_mant_len-32)); // Infinity
          u.eksplicit.mlo = 0;
        }
      else
        { u.eksplicit.semhi =
              ((sint32)sign & bit(31))                           /* Vorzeichen */
            | ((uint32)(lendiff+DF_exp_mid) << (DF_mant_len-32)) /* Exponent   */
            | ((uint32)manthi & (bit(DF_mant_len-32)-1));        /* Mantisse   */
          u.eksplicit.mlo = mantlo;
        }
      return u.machine_double;
      #endif
}}

}  // namespace cln
