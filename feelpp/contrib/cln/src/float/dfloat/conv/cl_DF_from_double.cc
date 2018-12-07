// cl_double_to_DF_pointer().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/dfloat/cl_DF.h"


// Implementation.

namespace cln {

cl_heap_dfloat* cl_double_to_DF_pointer (const double x)
{
      var union { dfloat eksplicit; double machine_double; } u;
      u.machine_double = x;
      var dfloat val = u.eksplicit;
      #if (cl_word_size==64)
      var uintL exp = (val >> DF_mant_len) & (bit(DF_exp_len)-1); // e
      if (exp == 0) // e=0 ?
        // vorzeichenbehaftete 0.0 oder subnormale Zahl
        { if (!((val << 1) == 0) && underflow_allowed())
            { throw floating_point_underflow_exception(); }
            else
            { return cl_DF_0; } // +/- 0.0 -> 0.0
        }
      elif (exp == 2047) // e=2047 ?
        { if (!((val << (64-DF_mant_len)) == 0))
            { throw floating_point_nan_exception(); } // NaN
            else
            { throw floating_point_overflow_exception(); } // Infinity, Overflow
        }
      else
        { // Der Exponent muß um DF_exp_mid-1022 erhöht werden.
          if ((DF_exp_mid>1022) && (exp > DF_exp_high-DF_exp_mid+1022))
            { throw floating_point_overflow_exception(); } // Overflow
          val += (sint64)(DF_exp_mid - 1022) << DF_mant_len;
          return allocate_dfloat(val);
        }
      #else
      var uintL exp = (val.semhi >> (DF_mant_len-32)) & (bit(DF_exp_len)-1); // e
      if (exp == 0) // e=0 ?
        // vorzeichenbehaftete 0.0 oder subnormale Zahl
        { if (!(((val.semhi << 1) == 0) && (val.mlo == 0)) && underflow_allowed())
            { throw floating_point_underflow_exception(); }
            else
            { return cl_DF_0; } // +/- 0.0 -> 0.0
        }
      elif (exp == 2047) // e=2047 ?
        { if (!(((val.semhi << (64-DF_mant_len)) == 0) && (val.mlo == 0)))
            { throw floating_point_nan_exception(); } // NaN
            else
            { throw floating_point_overflow_exception(); } // Infinity, Overflow
        }
      else
        { // Der Exponent muß um DF_exp_mid-1022 erhöht werden.
          if ((DF_exp_mid>1022) && (exp > DF_exp_high-DF_exp_mid+1022))
            { throw floating_point_overflow_exception(); } // Overflow
          val.semhi += (sint32)(DF_exp_mid - 1022) << (DF_mant_len-32);
          return allocate_dfloat(val.semhi,val.mlo);
        }
      #endif
}

}  // namespace cln
