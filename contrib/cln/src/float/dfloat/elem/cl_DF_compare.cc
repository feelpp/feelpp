// compare().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

namespace cln {

cl_signean compare (const cl_DF& x, const cl_DF& y)
{
// Methode:
// x und y haben verschiedenes Vorzeichen ->
//    x < 0 -> x < y
//    x >= 0 -> x > y
// x und y haben gleiches Vorzeichen ->
//    x >=0 -> vergleiche x und y (die rechten 53 Bits)
//    x <0 -> vergleiche y und x (die rechten 53 Bits)
#if (cl_word_size==64)
      var dfloat x_ = TheDfloat(x)->dfloat_value;
      var dfloat y_ = TheDfloat(y)->dfloat_value;
      if ((sint64)y_ >= 0)
        // y>=0
        { if ((sint64)x_ >= 0)
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
        { if ((sint64)x_ >= 0)
            // y<0, x>=0
            { return signean_plus; } // x>y
            else
            // y<0, x<0
            { if (x_ > y_) return signean_minus; // |x|>|y| -> x<y
              if (x_ < y_) return signean_plus; // |x|<|y| -> x>y
              return signean_null;
            }
        }
#else
      var uint32 x_semhi = TheDfloat(x)->dfloat_value.semhi;
      var uint32 y_semhi = TheDfloat(y)->dfloat_value.semhi;
      var uint32 x_mlo = TheDfloat(x)->dfloat_value.mlo;
      var uint32 y_mlo = TheDfloat(y)->dfloat_value.mlo;
      if ((sint32)y_semhi >= 0)
        // y>=0
        { if ((sint32)x_semhi >= 0)
            // y>=0, x>=0
            { if (x_semhi < y_semhi) return signean_minus; // x<y
              if (x_semhi > y_semhi) return signean_plus; // x>y
              if (x_mlo < y_mlo) return signean_minus; // x<y
              if (x_mlo > y_mlo) return signean_plus; // x>y
              return signean_null;
            }
            else
            // y>=0, x<0
            { return signean_minus; } // x<y
        }
        else
        { if ((sint32)x_semhi >= 0)
            // y<0, x>=0
            { return signean_plus; } // x>y
            else
            // y<0, x<0
            { if (x_semhi > y_semhi) return signean_minus; // |x|>|y| -> x<y
              if (x_semhi < y_semhi) return signean_plus; // |x|<|y| -> x>y
              if (x_mlo > y_mlo) return signean_minus; // |x|>|y| -> x<y
              if (x_mlo < y_mlo) return signean_plus; // |x|<|y| -> x>y
              return signean_null;
            }
        }
#endif
}

}  // namespace cln
