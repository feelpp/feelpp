// fround().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"

namespace cln {

const cl_SF fround (const cl_SF& x)
{
// Methode:
// x = 0.0 oder e<0 -> Ergebnis 0.0
// 0<=e<=16 -> letzte (17-e) Bits der Mantisse wegrunden,
//             Exponent und Vorzeichen beibehalten.
// e>16 -> Ergebnis x
      var uintL uexp = SF_uexp(x); // e + SF_exp_mid
      if (uexp < SF_exp_mid) // x = 0.0 oder e<0 ?
        { return SF_0; }
        else
        { if (uexp > SF_exp_mid+SF_mant_len) // e > 16 ?
            { return x; }
            else
            if (uexp > SF_exp_mid+1) // e>1 ?
              { var cl_uint bitmask = // Bitmaske: Bit 16-e gesetzt, alle anderen gelöscht
                  bit(SF_mant_len+SF_mant_shift + SF_exp_mid-uexp);
                var cl_uint mask = // Bitmaske: Bits 15-e..0 gesetzt, alle anderen gelöscht
                  bitmask - bit(SF_mant_shift);
                if ( ((x.word & bitmask) ==0) // Bit 16-e =0 -> abrunden
                     || ( ((x.word & mask) ==0) // Bit 16-e =1 und Bits 15-e..0 >0 -> aufrunden
                          // round-to-even, je nach Bit 17-e :
                          && ((x.word & (bitmask<<1)) ==0)
                   )    )
                  // abrunden
                  { mask |= bitmask; // Bitmaske: Bits 16-e..0 gesetzt, alle anderen gelöscht
                    return cl_SF_from_word(x.word & ~mask);
                  }
                  else
                  // aufrunden
                  { return cl_SF_from_word(
                      (x.word | mask) // alle diese Bits 15-e..0 setzen (Bit 16-e schon gesetzt)
                      + bit(SF_mant_shift) // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                      );
                  }
              }
            elif (uexp == SF_exp_mid+1) // e=1 ?
              // Wie bei 1 < e <= 16, nur daß Bit 17-e stets gesetzt ist.
              { if ((x.word & bit(SF_mant_len+SF_mant_shift-1)) ==0) // Bit 16-e =0 -> abrunden
                  // abrunden
                  { return cl_SF_from_word(x.word & ~(bit(SF_mant_len+SF_mant_shift)-bit(SF_mant_shift))); }
                  else
                  // aufrunden
                  { return cl_SF_from_word(
                      (x.word | (bit(SF_mant_len+SF_mant_shift)-bit(SF_mant_shift))) // alle diese Bits 16-e..0 setzen
                      + bit(SF_mant_shift) // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                      );
                  }
              }
            else // e=0 ?
              // Wie bei 1 < e <= 16, nur daß Bit 16-e stets gesetzt
              // und Bit 17-e stets gelöscht ist.
              { if ((x.word & (bit(SF_mant_len+SF_mant_shift)-bit(SF_mant_shift))) ==0)
                  // abrunden von +-0.5 zu 0.0
                  { return SF_0; }
                  else
                  // aufrunden
                  { return cl_SF_from_word(
                      (x.word | (bit(SF_mant_len+SF_mant_shift)-bit(SF_mant_shift))) // alle Bits 15-e..0 setzen
                      + bit(SF_mant_shift) // letzte Stelle erhöhen, dabei Exponenten incrementieren
                      );
              }   }
        }
}

}  // namespace cln
