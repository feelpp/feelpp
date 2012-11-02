// fround().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

namespace cln {

const cl_FF fround (const cl_FF& x)
{
// Methode:
// x = 0.0 oder e<0 -> Ergebnis 0.0
// 0<=e<=23 -> letzte (24-e) Bits der Mantisse wegrunden,
//             Exponent und Vorzeichen beibehalten.
// e>23 -> Ergebnis x
      var ffloat x_ = cl_ffloat_value(x);
      var uintL uexp = FF_uexp(x_); // e + FF_exp_mid
      if (uexp < FF_exp_mid) // x = 0.0 oder e<0 ?
        { return cl_FF_0; }
        else
        { if (uexp > FF_exp_mid+FF_mant_len) // e > 23 ?
            { return x; }
            else
            if (uexp > FF_exp_mid+1) // e>1 ?
              { var uint32 bitmask = // Bitmaske: Bit 23-e gesetzt, alle anderen gelöscht
                  bit(FF_mant_len+FF_exp_mid-uexp);
                var uint32 mask = // Bitmaske: Bits 22-e..0 gesetzt, alle anderen gelöscht
                  bitmask-1;
                if ( ((x_ & bitmask) ==0) // Bit 23-e =0 -> abrunden
                     || ( ((x_ & mask) ==0) // Bit 23-e =1 und Bits 22-e..0 >0 -> aufrunden
                          // round-to-even, je nach Bit 24-e :
                          && ((x_ & (bitmask<<1)) ==0)
                   )    )
                  // abrunden
                  { mask |= bitmask; // Bitmaske: Bits 23-e..0 gesetzt, alle anderen gelöscht
                    return allocate_ffloat( x_ & ~mask );
                  }
                  else
                  // aufrunden
                  { return allocate_ffloat
                      ((x_ | mask) // alle diese Bits 22-e..0 setzen (Bit 23-e schon gesetzt)
                       + 1 // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                      );
                  }
              }
            elif (uexp == FF_exp_mid+1) // e=1 ?
              // Wie bei 1 < e <= 23, nur daß Bit 24-e stets gesetzt ist.
              { if ((x_ & bit(FF_mant_len-1)) ==0) // Bit 23-e =0 -> abrunden
                  // abrunden
                  { return allocate_ffloat( x_ & ~(bit(FF_mant_len)-1) ); }
                  else
                  // aufrunden
                  { return allocate_ffloat
                      ((x_ | (bit(FF_mant_len)-1)) // alle diese Bits 23-e..0 setzen
                       + 1 // letzte Stelle erhöhen, dabei evtl. Exponenten incrementieren
                      );
                  }
              }
            else // e=0 ?
              // Wie bei 1 < e <= 23, nur daß Bit 23-e stets gesetzt
              // und Bit 24-e stets gelöscht ist.
              { if ((x_ & (bit(FF_mant_len)-1)) ==0)
                  // abrunden von +-0.5 zu 0.0
                  { return cl_FF_0; }
                  else
                  // aufrunden
                  { return allocate_ffloat
                      ((x_ | (bit(FF_mant_len)-1)) // alle Bits 22-e..0 setzen
                       + 1 // letzte Stelle erhöhen, dabei Exponenten incrementieren
                      );
              }   }
        }
}

}  // namespace cln
