// dpb().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"

namespace cln {

const cl_I dpb (const cl_I& newbyte, const cl_I& n, const cl_byte& b)
{
      // Methode:
      // (DPB newbyte (byte s p) integer)
      // = (DEPOSIT-FIELD (ASH newbyte p) (byte s p) integer)
      return deposit_field(ash(newbyte,b.position),n,b);
}

}  // namespace cln
