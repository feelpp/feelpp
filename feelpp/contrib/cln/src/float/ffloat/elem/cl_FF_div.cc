// binary operator /

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "base/cl_N.h"
#include "float/cl_F.h"
#include "base/cl_low.h"

#include "base/cl_inline.h"
#include "float/ffloat/elem/cl_FF_zerop.cc"

namespace cln {


const cl_FF operator/ (const cl_FF& x1, const cl_FF& x2)
{
// Methode:
// x2 = 0.0 -> Error
// x1 = 0.0 -> Ergebnis 0.0
// Sonst:
// Ergebnis-Vorzeichen = xor der beiden Vorzeichen von x1 und x2
// Ergebnis-Exponent = Differenz der beiden Exponenten von x1 und x2
// Ergebnis-Mantisse = Mantisse mant1 / Mantisse mant2, gerundet.
//   mant1/mant2 > 1/2, mant1/mant2 < 2;
//   nach Rundung mant1/mant2 >=1/2, <=2*mant1<2.
//   Bei mant1/mant2 >=1 brauche 23 Nachkommabits,
//   bei mant1/mant2 <1 brauche 24 Nachkommabits.
//   Fürs Runden: brauche ein Rundungsbit (Rest gibt an, ob exakt).
//   Brauche daher insgesamt 25 Nachkommabits von mant1/mant2.
//   Dividiere daher (als Unsigned Integers) 2^25*(2^24*mant1) durch (2^24*mant2).
//   Falls der Quotient >=2^25 ist, runde die letzten zwei Bits weg und
//     erhöhe den Exponenten um 1.
//   Falls der Quotient <2^25 ist, runde das letzte Bit weg. Bei rounding
//     overflow schiebe um ein weiteres Bit nach rechts, incr. Exponenten.
  #if defined(FAST_FLOAT) && !defined(__i386__)
      float_to_FF(FF_to_float(x1) / FF_to_float(x2), return ,
                  TRUE, TRUE, // Overflow und subnormale Zahl abfangen
                  !zerop_inline(x1), // ein Ergebnis +/- 0.0
                              // ist genau dann in Wirklichkeit ein Underflow
                  zerop_inline(x2), // Division durch Null abfangen
                  FALSE // kein NaN als Ergebnis möglich
                 );
  #else
      // x1,x2 entpacken:
      var cl_signean sign1;
      var sintL exp1;
      var uintL mant1;
      var cl_signean sign2;
      var sintL exp2;
      var uintL mant2;
      FF_decode(x2, { throw division_by_0_exception(); }, sign2=,exp2=,mant2=);
      FF_decode(x1, { return x1; }, sign1=,exp1=,mant1=);
      exp1 = exp1 - exp2; // Differenz der Exponenten
      sign1 = sign1 ^ sign2; // Ergebnis-Vorzeichen
      // Dividiere 2^25*mant1 durch mant2 oder (äquivalent)
      // 2^i*2^25*mant1 durch 2^i*mant2 für irgendein i mit 0 <= i <= 32-24 :
      var uintL mant;
      var uintL rest;
      // wähle i = 32-(FF_mant_len+1), also i+(FF_mant_len+2) = 33.
      divu_6432_3232(mant1<<1,0, mant2<<(32-(FF_mant_len+1)), mant=,rest=);
      if (mant >= bit(FF_mant_len+2))
        // Quotient >=2^25 -> 2 Bits wegrunden
        { var uintL rounding_bits = mant & (bit(2)-1);
          exp1 += 1; // Exponenten incrementieren
          mant = mant >> 2;
          if ( (rounding_bits < bit(1)) // 00,01 werden abgerundet
               || ( (rounding_bits == bit(1)) // 10
                    && (rest == 0) // und genau halbzahlig
                    && ((mant & bit(0)) ==0) // -> round-to-even
             )    )
            // abrunden
            {}
            else
            // aufrunden
            { mant += 1; }
        }
        else
        // Quotient <2^25 -> 1 Bit wegrunden
        { var uintL rounding_bit = mant & bit(0);
          mant = mant >> 1;
          if ( (rounding_bit == 0) // 0 wird abgerundet
               || ( (rest == 0) // genau halbzahlig
                    && ((mant & bit(0)) ==0) // -> round-to-even
             )    )
            // abrunden
            {}
            else
            // aufrunden
            { mant += 1;
              if (mant >= bit(FF_mant_len+1)) // rounding overflow?
                { mant = mant>>1; exp1 = exp1+1; }
        }   }
      return encode_FF(sign1,exp1,mant);
  #endif
}

}  // namespace cln
