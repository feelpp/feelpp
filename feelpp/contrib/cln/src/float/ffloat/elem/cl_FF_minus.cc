// binary operator -

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

namespace cln {


const cl_FF operator- (const cl_FF& x1, const cl_FF& x2)
{
  #ifdef FAST_FLOAT
      float_to_FF(FF_to_float(x1) - FF_to_float(x2), return ,
                  TRUE, TRUE, // Overflow und subnormale Zahl abfangen
                  FALSE, // kein Underflow mit Ergebnis +/- 0.0 möglich
                         // (nach Definition der subnormalen Zahlen)
                  FALSE, FALSE // keine Singularität, kein NaN als Ergebnis möglich
                 );
  #else
      var ffloat x2_ = cl_ffloat_value(x2);
      if (FF_uexp(x2_) == 0)
        { return x1; }
        else
        { return x1 + allocate_ffloat(x2_ ^ bit(31)); }
  #endif
}

}  // namespace cln
