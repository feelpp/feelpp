// ftruncate().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

namespace cln {

const cl_DF ftruncate (const cl_DF& x)
{
// Methode:
// x = 0.0 oder e<=0 -> Ergebnis 0.0
// 1<=e<=52 -> letzte (53-e) Bits der Mantisse auf 0 setzen,
//             Exponent und Vorzeichen beibehalten
// e>=53 -> Ergebnis x
#if (cl_word_size==64)
      var dfloat x_ = TheDfloat(x)->dfloat_value;
      var uintL uexp = DF_uexp(x_); // e + DF_exp_mid
      if (uexp <= DF_exp_mid) // 0.0 oder e<=0 ?
        { return cl_DF_0; }
        else
        { if (uexp > DF_exp_mid+DF_mant_len) // e > 52 ?
            { return x; }
            else
            // 1<=e<=52
            { return allocate_dfloat
                ( x_ & // Bitmaske: Bits 52-e..0 gelöscht, alle anderen gesetzt
                  ~(bit(DF_mant_len+1+DF_exp_mid-uexp)-1)
                );
        }   }
#else
      var uint32 semhi = TheDfloat(x)->dfloat_value.semhi;
      var uint32 mlo = TheDfloat(x)->dfloat_value.mlo;
      var uintL uexp = DF_uexp(semhi); // e + DF_exp_mid
      if (uexp <= DF_exp_mid) // 0.0 oder e<=0 ?
        { return cl_DF_0; }
        else
        { if (uexp > DF_exp_mid+DF_mant_len) // e > 52 ?
            { return x; }
            else
            // 1<=e<=52
            if (uexp > DF_exp_mid+DF_mant_len+1-32) // e > 21 ?
              { return allocate_dfloat
                  ( semhi,
                    mlo & // Bitmaske: Bits 52-e..0 gelöscht, alle anderen gesetzt
                    ~(bit(DF_mant_len+1+DF_exp_mid-uexp)-1)
                  );
              }
              else
              { return allocate_dfloat
                  ( semhi & // Bitmaske: Bits 20-e..0 gelöscht, alle anderen gesetzt
                    ~(bit(DF_mant_len+1+DF_exp_mid-32-uexp)-1),
                    0
                  );
        }     }
#endif
}

}  // namespace cln
