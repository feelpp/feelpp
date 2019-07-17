// sqrt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "float/cl_F.h"
#include "base/cl_low.h"

namespace cln {

const cl_FF sqrt (const cl_FF& x)
{
// Methode:
// x = 0.0 -> Ergebnis 0.0
// Ergebnis-Vorzeichen := positiv,
// Ergebnis-Exponent := ceiling(e/2),
// Ergebnis-Mantisse:
//   Bilde aus [1,m22,...,m0,(26 Nullbits)] bei geradem e,
//         aus [0,1,m22,...,m0,(25 Nullbits)] bei ungeradem e
//   die Ganzzahl-Wurzel, eine 25-Bit-Zahl mit einer führenden 1.
//   Runde das letzte Bit weg:
//     Bit 0 = 0 -> abrunden,
//     Bit 0 = 1 und Wurzel exakt -> round-to-even,
//     Bit 0 = 1 und Rest >0 -> aufrunden.
//   Dabei um ein Bit nach rechts schieben.
//   Bei Aufrundung auf 2^24 (rounding overflow) Mantisse um 1 Bit nach rechts
//     schieben und Exponent incrementieren.
      // x entpacken:
      var sintL exp;
      var uint32 mant;
      FF_decode(x, { return x; }, ,exp=,mant=);
      // Um die 64-Bit-Ganzzahl-Wurzel ausnutzen zu können, fügen wir beim
      // Radikanden 39 bzw. 40 statt 25 bzw. 26 Nullbits an.
      if (exp & bit(0))
        // e ungerade
        { mant = mant << (31-(FF_mant_len+1)); exp = exp+1; }
        else
        // e gerade
        { mant = mant << (32-(FF_mant_len+1)); }
      exp = exp >> 1; // exp := exp/2
      var bool exactp;
      isqrt_64_32(mant,0, mant=,exactp=); // mant := isqrt(mant*2^32), eine 32-Bit-Zahl
      // Die hinteren 31-FF_mant_len Bits wegrunden:
      if ( ((mant & bit(30-FF_mant_len)) ==0) // Bit 7 =0 -> abrunden
           || ( ((mant & (bit(30-FF_mant_len)-1)) ==0) // Bit 7 =1 und Bits 6..0 >0 -> aufrunden
                && exactp                   // Bit 7 =1 und Bits 6..0 =0, aber Rest -> aufrunden
                // round-to-even, je nach Bit 8 :
                && ((mant & bit(31-FF_mant_len)) ==0)
         )    )
        // abrunden
        { mant = mant >> (31-FF_mant_len); }
        else
        // aufrunden
        { mant = mant >> (31-FF_mant_len);
          mant += 1;
          if (mant >= bit(FF_mant_len+1)) // rounding overflow?
            { mant = mant>>1; exp = exp+1; }
        }
      return encode_FF(0,exp,mant);
}

}  // namespace cln
