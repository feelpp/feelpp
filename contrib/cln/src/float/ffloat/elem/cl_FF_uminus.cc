// unary operator -

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

namespace cln {

const cl_FF operator- (const cl_FF& x)
{
// Methode:
// Falls x=0.0, fertig. Sonst Vorzeichenbit umdrehen.
      var ffloat x_ = cl_ffloat_value(x);
      if (FF_uexp(x_) == 0)
        return x;
      else
        return allocate_ffloat( x_ ^ bit(31) );
}

}  // namespace cln
