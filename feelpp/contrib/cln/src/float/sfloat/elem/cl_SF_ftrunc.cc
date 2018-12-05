// ftruncate().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"

namespace cln {

const cl_SF ftruncate (const cl_SF& x)
{
// Methode:
// x = 0.0 oder e<=0 -> Ergebnis 0.0
// 1<=e<=16 -> letzte (17-e) Bits der Mantisse auf 0 setzen,
//             Exponent und Vorzeichen beibehalten
// e>=17 -> Ergebnis x
      var uintL uexp = SF_uexp(x); // e + SF_exp_mid
      if (uexp <= SF_exp_mid) // 0.0 oder e<=0 ?
        { return SF_0; }
        else
        { if (uexp > SF_exp_mid+SF_mant_len) // e > 16 ?
            { return x; }
            else
            { return cl_SF_from_word(
                x.word & // Bitmaske: Bits 16-e..0 gelöscht, alle anderen gesetzt
                ~(bit(SF_mant_len+SF_mant_shift + 1+SF_exp_mid-uexp) - bit(SF_mant_shift))
                );
            }
        }
}

}  // namespace cln
