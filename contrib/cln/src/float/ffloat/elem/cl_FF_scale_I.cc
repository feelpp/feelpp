// scale_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "float/cl_F.h"
#include "integer/cl_I.h"

namespace cln {

const cl_FF scale_float (const cl_FF& x, const cl_I& delta)
{
  // Methode:
  // x=0.0 -> x als Ergebnis
  // delta muß ein Fixnum betragsmäßig <= FF_exp_high-FF_exp_low sein.
  // Neues FF mit um delta vergrößertem Exponenten bilden.
      // x entpacken:
      var cl_signean sign;
      var sintL exp;
      var uint32 mant;
      FF_decode(x, { return x; }, sign=,exp=,mant=);
      if (!minusp(delta))
        // delta>=0
        { var uintV udelta;
          if (fixnump(delta)
              && ((udelta = FN_to_V(delta)) <= (uintV)(FF_exp_high-FF_exp_low))
             )
            { exp = exp+udelta;
              return encode_FF(sign,exp,mant);
            }
            else
            { throw floating_point_overflow_exception(); }
        }
        else
        // delta<0
        { var uintV udelta;
          if (fixnump(delta)
              && ((udelta = -FN_to_V(delta)) <= (uintV)(FF_exp_high-FF_exp_low))
              && ((cl_value_len+1<intVsize) || !(udelta==0))
             )
            { exp = exp-udelta;
              return encode_FF(sign,exp,mant);
            }
            else
            if (underflow_allowed())
              { throw floating_point_underflow_exception(); }
              else
              { return cl_FF_0; }
        }
}

}  // namespace cln
