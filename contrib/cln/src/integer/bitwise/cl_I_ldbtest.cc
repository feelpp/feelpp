// ldb_test().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"
#include "integer/bitwise/cl_I_byte.h"

namespace cln {

bool ldb_test (const cl_I& n, const cl_byte& b)
{
      // Methode:
      // (ldb-test (byte s p) n)
      // Falls s=0: =0.
      // Falls s>0:
      //   l:=(integer-length n)
      //   Falls l <= p : Falls n>=0, =0, denn Bits p+s-1..p sind =0.
      //                  Falls n<0, /=0, denn Bits p+s-1..p sind =1.
      //   Falls p < l :
      //     Falls p+s>l, /=0, denn bei n>=0 ist Bit l-1 =1,
      //                       und bei n<0 sind Bits p+s-1..l =1.
      //     Falls p+s<=l,
      //       extrahiere die Bits p,...,p+s-1 von n und teste sie.
      var uintC s = b.size;
      var uintC p = b.position;
      if (s==0) return false;
      var uintC l = integer_length(n); // l = (integer-length n)
      if (l<=p)
        // l<=p
        if (!minusp(n))
          return false; // n>=0
          else
          return true; // n<0
        else
        // l>p
        { var uintC ps = p+s;
          if (ps>l) // p+s>l ?
            return true;
          // Bits p,...,q-1 mit q = min(p+s,l) = p+s extrahieren und testen:
          return ldb_extract_test(n,p,ps);
        }
}

}  // namespace cln
