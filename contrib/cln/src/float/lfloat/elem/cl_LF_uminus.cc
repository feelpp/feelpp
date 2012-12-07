// unary operator -

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_LF operator- (const cl_LF& x)
{
// Methode:
// Falls x=0.0, fertig. Sonst Vorzeichenbit umdrehen und Pointer beibehalten.
      if (TheLfloat(x)->expo == 0)
        { return x; }
        else
        { var uintC len = TheLfloat(x)->len;
          var Lfloat mx = allocate_lfloat(len, TheLfloat(x)->expo, ~ TheLfloat(x)->sign);
          copy_loop_up(&TheLfloat(x)->data[0],&TheLfloat(mx)->data[0],len);
          return mx;
        }
}

}  // namespace cln
