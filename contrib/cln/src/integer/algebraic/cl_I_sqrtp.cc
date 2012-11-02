// sqrtp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

bool sqrtp (const cl_I& x, cl_I* w)
{
// Methode:
// [Henri Cohen: A course in computational algebraic number theory, 2nd prnt.,
//  section 1.7.2.]
// If x is a square, it has not only to be ==0,1 mod 4, but also
//   - one of the 12 squares mod 64,
//   - one of the 16 squares mod 63,
//   - one of the 21 squares mod 65,
//   - one of the 6 squares mod 11.
// Then we check whether the ISQRT gives the remainder 0.
   static char squares_mod_11 [11] =
       // 0 1 3 4 5 9
       { 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0 };
   static char squares_mod_63 [63] =
       // 0 1 4 7 9 16 18 22 25 28 36 37 43 46 49 58
       { 1, 1, 0, 0, 1, 0, 0, 1, 0, 1,  0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
         0, 0, 1, 0, 0, 1, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
         0, 0, 0, 1, 0, 0, 1, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0
       };
    static char squares_mod_64 [64] =
       // 0 1 4 9 16 17 25 33 36 41 49 57
       { 1, 1, 0, 0, 1, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0
       };
    static char squares_mod_65 [65] =
       // 0 1 4 9 10 14 16 25 26 29 30 35 36 39 40 49 51 55 56 61 64
       { 1, 1, 0, 0, 1, 0, 0, 0, 0, 1,  1, 0, 0, 0, 1, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 1, 0, 0, 1,  1, 0, 0, 0, 0, 1, 1, 0, 0, 1,
         1, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 1, 0, 0, 0, 1, 1, 0, 0, 0,
         0, 1, 0, 0, 1
       };
    { CL_ALLOCA_STACK;
      var const uintD* x_MSDptr;
      var uintC x_len;
      var const uintD* x_LSDptr;
      I_to_NDS_nocopy(x, x_MSDptr=,x_len=,x_LSDptr=, // Digit sequence >=0 zu x
                      true, { *w = 0; return true; } // 0 is a square
                     );
      // Check mod 64.
      { var uintD lsd = lspref(x_LSDptr,0);
        if (!squares_mod_64[lsd & 63])
          { return false; } // not a square mod 64 -> not a square
      }
      // Check mod 63.
      { var cl_I_div_t div63 = cl_divide(x,L_to_FN(63));
        if (!squares_mod_63[FN_to_UV(div63.remainder)])
          { return false; } // not a square mod 63 -> not a square
      }
      // Check mod 65.
      { var cl_I_div_t div65 = cl_divide(x,L_to_FN(65));
        if (!squares_mod_65[FN_to_UV(div65.remainder)])
          { return false; } // not a square mod 65 -> not a square
      }
      // Check mod 11.
      { var cl_I_div_t div11 = cl_divide(x,L_to_FN(11));
        if (!squares_mod_11[FN_to_UV(div11.remainder)])
          { return false; } // not a square mod 11 -> not a square
      }
      // Check with full precision.
      { var DS y;
        var bool squarep;
        UDS_sqrt(x_MSDptr,x_len,x_LSDptr, &y, squarep=); // Wurzel ziehen
        if (squarep)
          { *w = NUDS_to_I(y.MSDptr,y.len); } // als Integer
        return squarep;
    } }
}
// Bit complexity (x of length N): O(M(N)).

}  // namespace cln
