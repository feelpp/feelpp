// compare().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "base/digitseq/cl_DS.h"

#include "base/cl_inline.h"
#include "float/lfloat/elem/cl_LF_minusp.cc"

namespace cln {

cl_signean compare (const cl_LF& x, const cl_LF& y)
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
            { // Vergleiche Exponenten und Mantissen:
              { var uintE x_uexp = TheLfloat(x)->expo;
                var uintE y_uexp = TheLfloat(y)->expo;
                if (x_uexp < y_uexp) return signean_minus; // x<y
                if (x_uexp > y_uexp) return signean_plus; // x>y
              }
              { var uintC x_len = TheLfloat(x)->len;
                var uintC y_len = TheLfloat(y)->len;
                var uintC len = (x_len<y_len ? x_len : y_len); // min(x_len,y_len)
                // len Digits vergleichen:
                var cl_signean erg =
                  compare_loop_msp(arrayMSDptr(TheLfloat(x)->data,x_len),arrayMSDptr(TheLfloat(y)->data,y_len),len);
                if (!(erg==0)) { return erg; } // verschieden -> fertig
                // gemeinsames Teilstück war gleich
                if (x_len == y_len) { return signean_null; } // gleiche Länge -> fertig
                if (x_len > y_len)
                  // x länger als y
                  { if (DS_test_loop(arrayMSDptr(TheLfloat(x)->data,x_len) mspop y_len,x_len-y_len,arrayLSDptr(TheLfloat(x)->data,x_len)))
                      { return signean_plus; } // x>y
                      else
                      { return signean_null; }
                  }
                  else
                  // y länger als x
                  { if (DS_test_loop(arrayMSDptr(TheLfloat(y)->data,y_len) mspop x_len,y_len-x_len,arrayLSDptr(TheLfloat(y)->data,y_len)))
                      { return signean_minus; } // x<y
                      else
                      { return signean_null; }
                  }
            } }
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
            { // Vergleiche Exponenten und Mantissen:
              { var uintE x_uexp = TheLfloat(x)->expo;
                var uintE y_uexp = TheLfloat(y)->expo;
                if (x_uexp < y_uexp) return signean_plus; // |x|<|y| -> x>y
                if (x_uexp > y_uexp) return signean_minus; // |x|>|y| -> x<y
              }
              { var uintC x_len = TheLfloat(x)->len;
                var uintC y_len = TheLfloat(y)->len;
                var uintC len = (x_len<y_len ? x_len : y_len); // min(x_len,y_len)
                // len Digits vergleichen:
                var cl_signean erg =
                  compare_loop_msp(arrayMSDptr(TheLfloat(y)->data,y_len),arrayMSDptr(TheLfloat(x)->data,x_len),len);
                if (!(erg==0)) { return erg; } // verschieden -> fertig
                // gemeinsames Teilstück war gleich
                if (x_len == y_len) { return signean_null; } // gleiche Länge -> fertig
                if (x_len > y_len)
                  // x länger als y
                  { if (DS_test_loop(arrayMSDptr(TheLfloat(x)->data,x_len) mspop y_len,x_len-y_len,arrayLSDptr(TheLfloat(x)->data,x_len)))
                      { return signean_minus; } // |x|>|y| -> x<y
                      else
                      { return signean_null; }
                  }
                  else
                  // y länger als x
                  { if (DS_test_loop(arrayMSDptr(TheLfloat(y)->data,y_len) mspop x_len,y_len-x_len,arrayLSDptr(TheLfloat(y)->data,y_len)))
                      { return signean_plus; } // |x|<|y| -> x>y
                      else
                      { return signean_null; }
                  }
            } }
        }
}

}  // namespace cln
