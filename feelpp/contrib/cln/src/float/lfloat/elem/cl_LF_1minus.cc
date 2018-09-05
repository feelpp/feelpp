// LF_LF_minus_LF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_LF LF_LF_minus_LF (const cl_LF& x1, const cl_LF& x2)
{
// Methode:
// (- x1 x2) = (+ x1 (- x2))
      if (TheLfloat(x2)->expo == 0)
        { return x1; }
        else
        { var uintC len2 = TheLfloat(x2)->len;
          var Lfloat mx2 = allocate_lfloat(len2, TheLfloat(x2)->expo, ~ TheLfloat(x2)->sign);
          copy_loop_up(&TheLfloat(x2)->data[0],&TheLfloat(mx2)->data[0],len2);
          return LF_LF_plus_LF(x1,mx2);
        }
}

}  // namespace cln
