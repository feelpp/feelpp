// compare().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "base/cl_inline.h"
#include "float/sfloat/elem/cl_SF_minusp.cc"

namespace cln {

cl_signean CL_FLATTEN compare (const cl_SF& x, const cl_SF& y)
{
// Methode:
// x und y haben verschiedenes Vorzeichen ->
//    x < 0 -> x < y
//    x >= 0 -> x > y
// x und y haben gleiches Vorzeichen ->
//    x >=0 -> vergleiche x und y (die rechten 24 Bits)
//    x <0 -> vergleiche y und x (die rechten 24 Bits)
      if (!minusp_inline(y))
        // y>=0
        { if (!minusp_inline(x))
            // y>=0, x>=0
            { if (x.word < y.word) return signean_minus; // x<y
              if (x.word > y.word) return signean_plus; // x>y
              return signean_null;
            }
            else
            // y>=0, x<0
            { return signean_minus; } // x<y
        }
        else
        { if (!minusp_inline(x))
            // y<0, x>=0
            { return signean_plus; } // x>y
            else
            // y<0, x<0
            { if (x.word > y.word) return signean_minus; // |x|>|y| -> x<y
              if (x.word < y.word) return signean_plus; // |x|<|y| -> x>y
              return signean_null;
            }
        }
}

}  // namespace cln
