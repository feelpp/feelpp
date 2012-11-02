// binary operator -

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

namespace cln {


const cl_DF operator- (const cl_DF& x1, const cl_DF& x2)
{
// Methode:
// (- x1 x2) = (+ x1 (- x2))
#ifdef FAST_DOUBLE
      double_to_DF(DF_to_double(x1) - DF_to_double(x2), return ,
                   TRUE, TRUE, // Overflow und subnormale Zahl abfangen
                   FALSE, // kein Underflow mit Ergebnis +/- 0.0 möglich
                          // (nach Definition der subnormalen Zahlen)
                   FALSE, FALSE // keine Singularität, kein NaN als Ergebnis möglich
                  );
#else
#if (cl_word_size==64)
      var dfloat x2_ = TheDfloat(x2)->dfloat_value;
      if (DF_uexp(x2_) == 0)
        { return x1; }
        else
        { return x1 + allocate_dfloat(x2_ ^ bit(63)); }
#else
      var uint32 x2_semhi = TheDfloat(x2)->dfloat_value.semhi;
      var uint32 x2_mlo = TheDfloat(x2)->dfloat_value.mlo;
      if (DF_uexp(x2_semhi) == 0)
        { return x1; }
        else
        { return x1 + allocate_dfloat(x2_semhi ^ bit(31), x2_mlo); }
#endif
#endif
}

}  // namespace cln
