// binary operator *

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "base/cl_low.h"

namespace cln {

const cl_SF operator* (const cl_SF& x1, const cl_SF& x2)
{
// Methode:
// Falls x1=0.0 oder x2=0.0 -> Ergebnis 0.0
// Sonst: Ergebnis-Vorzeichen = VZ von x1 xor VZ von x2.
//        Ergebnis-Exponent = Summe der Exponenten von x1 und x2.
//        Ergebnis-Mantisse = Produkt der Mantissen von x1 und x2, gerundet:
//          2^-17 * (2^16 + m1)  *  2^-17 * (2^16 + m2)
//          = 2^-34 * (2^32 + 2^16*m1 + 2^16*m2 + m1*m2),
//          die Klammer ist >=2^32, <=(2^17-1)^2<2^34 .
//          Falls die Klammer >=2^33 ist, um 17 Bit nach rechts schieben und
//            runden: Falls Bit 16 Null, abrunden; falls Bit 16 Eins und
//            Bits 15..0 alle Null, round-to-even; sonst aufrunden.
//          Falls die Klammer <2^33 ist, um 16 Bit nach rechts schieben und
//            runden: Falls Bit 15 Null, abrunden; falls Bit 15 Eins und
//            Bits 14..0 alle Null, round-to-even; sonst aufrunden. Nach
//            Aufrunden: Falls =2^17, um 1 Bit nach rechts schieben. Sonst
//            Exponenten um 1 erniedrigen.
      // x1,x2 entpacken:
      var cl_signean sign1;
      var sintL exp1;
      var uintL mant1;
      var cl_signean sign2;
      var sintL exp2;
      var uintL mant2;
      SF_decode(x1, { return x1; }, sign1=,exp1=,mant1=);
      SF_decode(x2, { return x2; }, sign2=,exp2=,mant2=);
      exp1 = exp1 + exp2; // Summe der Exponenten
      sign1 = sign1 ^ sign2; // Ergebnis-Vorzeichen
      var uintL manthi;
      var uintL mantlo;
      // Mantissen mant1 und mant2 multiplizieren:
      #if (SF_mant_len<16)
      mantlo = mulu16(mant1,mant2);
      manthi = mantlo >> SF_mant_len;
      mantlo = mantlo & (bit(SF_mant_len)-1);
      #elif (SF_mant_len==16)
      manthi = mulu16(low16(mant1),low16(mant2));
      mantlo = low16(manthi);
      manthi = (uint32)(high16(manthi)) + (uint32)(low16(mant1)) + mant2;
      #else // (SF_mant_len>16)
      mulu24(mant1,mant2, manthi=,mantlo=);
      manthi = (manthi << (32-SF_mant_len)) | (mantlo >> SF_mant_len);
      mantlo = mantlo & (bit(SF_mant_len)-1);
      #endif
      // Nun ist 2^SF_mant_len * manthi + mantlo = mant1 * mant2.
      if (manthi >= bit(SF_mant_len+1))
        // mant1*mant2 >= 2^(2*SF_mant_len+1)
        { if ( ((manthi & bit(0)) ==0) // Bit SF_mant_len =0 -> abrunden
               || ( (mantlo ==0) // Bit SF_mant_len =1 und Bits SF_mant_len-1..0 >0 -> aufrunden
                    // round-to-even, je nach Bit SF_mant_len+1 :
                    && ((manthi & bit(1)) ==0)
             )    )
            // abrunden
            { manthi = manthi >> 1; goto ab; }
            else
            // aufrunden
            { manthi = manthi >> 1; goto auf; }
        }
        else
        // mant1*mant2 < 2^(2*SF_mant_len+1)
        { exp1 = exp1-1; // Exponenten decrementieren
          if ( ((mantlo & bit(SF_mant_len-1)) ==0) // Bit SF_mant_len-1 =0 -> abrunden
               || ( ((mantlo & (bit(SF_mant_len-1)-1)) ==0) // Bit SF_mant_len-1 =1 und Bits SF_mant_len-2..0 >0 -> aufrunden
                    // round-to-even, je nach Bit SF_mant_len :
                    && ((manthi & bit(0)) ==0)
             )    )
            // abrunden
            goto ab;
            else
            // aufrunden
            goto auf;
        }
      auf:
      manthi = manthi+1;
      // Hier ist 2^SF_mant_len <= manthi <= 2^(SF_mant_len+1)
      if (manthi >= bit(SF_mant_len+1)) // rounding overflow?
        { manthi = manthi>>1; exp1 = exp1+1; } // Shift nach rechts
      ab:
      // Runden fertig, 2^SF_mant_len <= manthi < 2^(SF_mant_len+1)
      return encode_SF(sign1,exp1,manthi);
}

}  // namespace cln
