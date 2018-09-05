// mask_field().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"
#include "integer/bitwise/cl_I_byte.h"

namespace cln {

const cl_I mask_field (const cl_I& n, const cl_byte& b)
{
      // Methode:
      // (mask-field (byte s p) n) extrahiere die Bits p,...,p+s-1 von n.
      // l:=(integer-length n)
      // Falls l <= p :
      //   Falls n>=0: 0, falls n<0: 2^(p+s) - 2^p (s Einsenbits).
      // Falls p <= l :
      //   q:=min(p+s,l).
      //   Extrahiere die Bits p,...,q-1 von n.
      //   Falls p+s>l und n<0, füge p+s-l Einsenbits an (addiere 2^(p+s)-2^l).
      var uintC s = b.size;
      var uintC p = b.position;
     {var uintC ps = p+s;
      var uintC l = integer_length(n); // l = (integer-length n)
      if (l<=p)
        // l<=p
        if (!minusp(n))
          // n>=0
          return 0; // 0 als Ergebnis
          else
          // n<0
          return cl_fullbyte(p,ps); // 2^(p+s)-2^p als Ergebnis
        else
        // l>p
        { // Bits p,...,q-1 mit q = min(p+s,l) extrahieren:
          var cl_I erg = mkf_extract(n,p,(ps<l ? ps : l));
          if ((ps>l) && minusp(n)) // p+s>l und n<0 ?
            { return logior(erg,cl_fullbyte(l,ps)); } // setze Bits l,...,p+s-1
            // (logisches Exklusiv-Oder oder Addition ginge auch)
            else
            return erg;
     }  }
}

}  // namespace cln
