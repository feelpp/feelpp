// ftruncate().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

namespace cln {

const cl_FF ftruncate (const cl_FF& x)
{
// Methode:
// x = 0.0 oder e<=0 -> Ergebnis 0.0
// 1<=e<=23 -> letzte (24-e) Bits der Mantisse auf 0 setzen,
//             Exponent und Vorzeichen beibehalten
// e>=24 -> Ergebnis x
      var ffloat x_ = cl_ffloat_value(x);
      var uintL uexp = FF_uexp(x_); // e + FF_exp_mid
      if (uexp <= FF_exp_mid) // 0.0 oder e<=0 ?
        { return cl_FF_0; }
        else
        { if (uexp > FF_exp_mid+FF_mant_len) // e > 23 ?
            { return x; }
            else
            { return allocate_ffloat
                ( x_ & // Bitmaske: Bits 23-e..0 gelöscht, alle anderen gesetzt
                  ~(bit(FF_mant_len+1+FF_exp_mid-uexp)-1)
                );
        }   }
}

}  // namespace cln
