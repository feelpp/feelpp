// deposit_field().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"

namespace cln {

const cl_I deposit_field (const cl_I& newbyte, const cl_I& n, const cl_byte& b)
{
      // Methode:
      // (DEPOSIT-FIELD newbyte (byte s p) integer)
      //  = (logxor integer
      //            (ash (logxor (ldb (byte s p) newbyte) (ldb (byte s p) integer))
      //                 p
      //    )       )
      return logxor(n, ash(logxor(ldb(newbyte,b),ldb(n,b)), b.position));
}

}  // namespace cln
