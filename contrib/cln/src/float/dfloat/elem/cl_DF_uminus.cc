// unary operator -

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

namespace cln {

const cl_DF operator- (const cl_DF& x)
{
// Methode:
// Falls x=0.0, fertig. Sonst Vorzeichenbit umdrehen.
#if (cl_word_size==64)
      var dfloat x_ = TheDfloat(x)->dfloat_value;
      if (DF_uexp(x_) == 0)
        return x;
      else
        return allocate_dfloat( x_ ^ bit(63) );
#else
      var uint32 semhi = TheDfloat(x)->dfloat_value.semhi;
      var uint32 mlo = TheDfloat(x)->dfloat_value.mlo;
      if (DF_uexp(semhi) == 0)
        return x;
      else
        return allocate_dfloat( semhi ^ bit(31), mlo );
#endif
}

}  // namespace cln
