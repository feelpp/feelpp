// binary operator *

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "float/cl_F.h"
#include "base/cl_low.h"

#include "base/cl_inline.h"
#include "float/ffloat/elem/cl_FF_zerop.cc"

namespace cln {


const cl_FF operator* (const cl_FF& x1, const cl_FF& x2)
{
// Methode:
// Falls x1=0.0 oder x2=0.0 -> Ergebnis 0.0
// Sonst: Ergebnis-Vorzeichen = VZ von x1 xor VZ von x2.
//        Ergebnis-Exponent = Summe der Exponenten von x1 und x2.
//        Ergebnis-Mantisse = Produkt der Mantissen von x1 und x2, gerundet:
//          2^-24 * mant1  *  2^-24 * mant2  =  2^-48 * (mant1*mant2),
//          die Klammer ist >=2^46, <=(2^24-1)^2<2^48 .
//          Falls die Klammer >=2^47 ist, um 24 Bit nach rechts schieben und
//            runden: Falls Bit 23 Null, abrunden; falls Bit 23 Eins und
//            Bits 22..0 alle Null, round-to-even; sonst aufrunden.
//          Falls die Klammer <2^47 ist, um 23 Bit nach rechts schieben und
//            runden: Falls Bit 22 Null, abrunden; falls Bit 22 Eins und
//            Bits 21..0 alle Null, round-to-even; sonst aufrunden. Nach
//            Aufrunden: Falls =2^24, um 1 Bit nach rechts schieben. Sonst
//            Exponenten um 1 erniedrigen.
  #ifdef FAST_FLOAT
      float_to_FF(FF_to_float(x1) * FF_to_float(x2), return ,
                  TRUE, TRUE, // Overflow und subnormale Zahl abfangen
                  !(zerop_inline(x1) || zerop_inline(x2)), // ein Ergebnis +/- 0.0
                              // ist genau dann in Wirklichkeit ein Underflow
                  FALSE, FALSE // keine Singularität, kein NaN als Ergebnis möglich
                 );
  #else
      // x1,x2 entpacken:
      var cl_signean sign1;
      var sintL exp1;
      var uintL mant1;
      var cl_signean sign2;
      var sintL exp2;
      var uintL mant2;
      FF_decode(x1, { return x1; }, sign1=,exp1=,mant1=);
      FF_decode(x2, { return x2; }, sign2=,exp2=,mant2=);
      exp1 = exp1 + exp2; // Summe der Exponenten
      sign1 = sign1 ^ sign2; // Ergebnis-Vorzeichen
      var uintL manthi;
      var uintL mantlo;
      // Mantissen mant1 und mant2 multiplizieren:
      mulu24(mant1,mant2, manthi=,mantlo=);
      manthi = (manthi << (32-FF_mant_len)) | (mantlo >> FF_mant_len);
      mantlo = mantlo & (bit(FF_mant_len)-1);
      // Nun ist 2^FF_mant_len * manthi + mantlo = mant1 * mant2.
      if (manthi >= bit(FF_mant_len+1))
        // mant1*mant2 >= 2^(2*FF_mant_len+1)
        { if ( ((manthi & bit(0)) ==0) // Bit FF_mant_len =0 -> abrunden
               || ( (mantlo ==0) // Bit FF_mant_len =1 und Bits FF_mant_len-1..0 >0 -> aufrunden
                    // round-to-even, je nach Bit FF_mant_len+1 :
                    && ((manthi & bit(1)) ==0)
             )    )
            // abrunden
            { manthi = manthi >> 1; goto ab; }
            else
            // aufrunden
            { manthi = manthi >> 1; goto auf; }
        }
        else
        // mant1*mant2 < 2^(2*FF_mant_len+1)
        { exp1 = exp1-1; // Exponenten decrementieren
          if ( ((mantlo & bit(FF_mant_len-1)) ==0) // Bit FF_mant_len-1 =0 -> abrunden
               || ( ((mantlo & (bit(FF_mant_len-1)-1)) ==0) // Bit FF_mant_len-1 =1 und Bits FF_mant_len-2..0 >0 -> aufrunden
                    // round-to-even, je nach Bit FF_mant_len :
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
      // Hier ist 2^FF_mant_len <= manthi <= 2^(FF_mant_len+1)
      if (manthi >= bit(FF_mant_len+1)) // rounding overflow?
        { manthi = manthi>>1; exp1 = exp1+1; } // Shift nach rechts
      ab:
      // Runden fertig, 2^FF_mant_len <= manthi < 2^(FF_mant_len+1)
      return encode_FF(sign1,exp1,manthi);
  #endif
}

}  // namespace cln
