// compare().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

namespace cln {

cl_signean compare (const cl_FF& x, const cl_FF& y)
{
// Methode:
// x und y haben verschiedenes Vorzeichen ->
//    x < 0 -> x < y
//    x >= 0 -> x > y
// x und y haben gleiches Vorzeichen ->
//    x >=0 -> vergleiche x und y (die rechten 24 Bits)
//    x <0 -> vergleiche y und x (die rechten 24 Bits)
      var uint32 x_ = cl_ffloat_value(x);
      var uint32 y_ = cl_ffloat_value(y);
      if ((sint32)y_ >= 0)
        // y>=0
        { if ((sint32)x_ >= 0)
            // y>=0, x>=0
            { if (x_ < y_) return signean_minus; // x<y
              if (x_ > y_) return signean_plus; // x>y
              return signean_null;
            }
            else
            // y>=0, x<0
            { return signean_minus; } // x<y
        }
        else
        { if ((sint32)x_ >= 0)
            // y<0, x>=0
            { return signean_plus; } // x>y
            else
            // y<0, x<0
            { if (x_ > y_) return signean_minus; // |x|>|y| -> x<y
              if (x_ < y_) return signean_plus; // |x|<|y| -> x>y
              return signean_null;
            }
        }
}

}  // namespace cln
