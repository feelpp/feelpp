// scale_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "float/cl_F.h"
#include "integer/cl_I.h"

namespace cln {

const cl_SF scale_float (const cl_SF& x, const cl_I& delta)
{
  // Methode:
  // x=0.0 -> x als Ergebnis
  // delta muß ein Fixnum betragsmäßig <= SF_exp_high-SF_exp_low sein.
  // Neues SF mit um delta vergrößertem Exponenten bilden.
      // x entpacken:
      var cl_signean sign;
      var sintL exp;
      var uint32 mant;
      SF_decode(x, { return x; }, sign=,exp=,mant=);
      if (!minusp(delta))
        // delta>=0
        { var uintV udelta;
          if (fixnump(delta)
              && ((udelta = FN_to_V(delta)) <= (uintV)(SF_exp_high-SF_exp_low))
             )
            { exp = exp+udelta;
              return encode_SF(sign,exp,mant);
            }
            else
            { throw floating_point_overflow_exception(); }
        }
        else
        // delta<0
        { var uintV udelta;
          if (fixnump(delta)
              && ((udelta = -FN_to_V(delta)) <= (uintV)(SF_exp_high-SF_exp_low))
              && ((cl_value_len+1<intVsize) || !(udelta==0))
             )
            { exp = exp-udelta;
              return encode_SF(sign,exp,mant);
            }
            else
            if (underflow_allowed())
              { throw floating_point_underflow_exception(); }
              else
              { return SF_0; }
        }
}

}  // namespace cln
