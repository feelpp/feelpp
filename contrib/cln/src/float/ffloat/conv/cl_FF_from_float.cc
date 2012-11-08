// cl_float_to_FF_pointer().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/ffloat/cl_FF.h"

namespace cln {

// Implementation.

cl_private_thing cl_float_to_FF_pointer (const float x)
{
      var union { ffloat eksplicit; float machine_float; } u;
      u.machine_float = x;
      var ffloat val = u.eksplicit;
      var uintL exp = (val >> FF_mant_len) & (bit(FF_exp_len)-1); // e
      if (exp == 0) // e=0 ?
        // vorzeichenbehaftete 0.0 oder subnormale Zahl
        { if (!((val << 1) == 0) && underflow_allowed())
            { throw floating_point_underflow_exception(); }
            else
            { return as_cl_private_thing(cl_FF_0); } // +/- 0.0 -> 0.0
        }
      elif (exp == 255) // e=255 ?
        { if (!((val << (32-FF_mant_len)) == 0))
            { throw floating_point_nan_exception(); } // NaN
            else
            { throw floating_point_overflow_exception(); } // Infinity, Overflow
        }
      else
        { // Der Exponent muß um FF_exp_mid-126 erhöht werden.
          if ((FF_exp_mid>126) && (exp > FF_exp_high-FF_exp_mid+126))
            { throw floating_point_overflow_exception(); } // Overflow
          val += (FF_exp_mid - 126) << FF_mant_len;
          #if defined(CL_WIDE_POINTERS)
          return as_cl_private_thing(allocate_ffloat(val));
          #else
          return (cl_private_thing)allocate_ffloat(val);
          #endif
        }
}

}  // namespace cln
