// scale_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "float/cl_F.h"
#include "integer/cl_I.h"

namespace cln {

const cl_DF scale_float (const cl_DF& x, const cl_I& delta)
{
  // Methode:
  // x=0.0 -> x als Ergebnis
  // delta muß ein Fixnum betragsmäßig <= DF_exp_high-DF_exp_low sein.
  // Neues DF mit um delta vergrößertem Exponenten bilden.
      // x entpacken:
      var cl_signean sign;
      var sintL exp;
#if (cl_word_size==64)
      var uint64 mant;
      DF_decode(x, { return x; }, sign=,exp=,mant=);
#else
      var uint32 manthi;
      var uint32 mantlo;
      DF_decode2(x, { return x; }, sign=,exp=,manthi=,mantlo=);
#endif
      if (!minusp(delta))
        // delta>=0
        { var uintV udelta;
          if (fixnump(delta)
              && ((udelta = FN_to_V(delta)) <= (uintV)(DF_exp_high-DF_exp_low))
             )
            { exp = exp+udelta;
#if (cl_word_size==64)
              return encode_DF(sign,exp,mant);
#else
              return encode_DF(sign,exp,manthi,mantlo);
#endif
            }
            else
            { throw floating_point_overflow_exception(); }
        }
        else
        // delta<0
        { var uintL udelta;
          if (fixnump(delta)
              && ((udelta = -FN_to_V(delta)) <= (uintV)(DF_exp_high-DF_exp_low))
              && ((cl_value_len+1<intVsize) || !(udelta==0))
             )
            { exp = exp-udelta;
#if (cl_word_size==64)
              return encode_DF(sign,exp,mant);
#else
              return encode_DF(sign,exp,manthi,mantlo);
#endif
            }
            else
            if (underflow_allowed())
              { throw floating_point_underflow_exception(); }
              else
              { return cl_DF_0; }
        }
}

}  // namespace cln
