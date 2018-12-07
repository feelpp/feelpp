// scale_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "float/cl_F.h"

namespace cln {

const cl_FF scale_float (const cl_FF& x, sintC delta)
{
  // Methode:
  // x=0.0 -> x als Ergebnis
  // delta muß betragsmäßig <= FF_exp_high-FF_exp_low sein.
  // Neues FF mit um delta vergrößertem Exponenten bilden.
      // x entpacken:
      var cl_signean sign;
      var sintL exp;
      var uint32 mant;
      FF_decode(x, { return x; }, sign=,exp=,mant=);
      if (delta >= 0)
        // delta>=0
        { var uintC udelta = delta;
          if (udelta <= (uintL)(FF_exp_high-FF_exp_low))
            { exp = exp+udelta;
              return encode_FF(sign,exp,mant);
            }
            else
            { throw floating_point_overflow_exception(); }
        }
        else
        // delta<0
        { var uintC udelta = -delta;
          if (udelta <= (uintL)(FF_exp_high-FF_exp_low))
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
