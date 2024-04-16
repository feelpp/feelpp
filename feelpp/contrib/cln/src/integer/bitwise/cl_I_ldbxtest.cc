// ldb_extract_test().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/bitwise/cl_I_byte.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

bool ldb_extract_test (const cl_I& x, uintC p, uintC q)
    { var const uintD* MSDptr;
      var uintC len;
      var const uintD* LSDptr;
      I_to_NDS_nocopy(x, MSDptr=,len=,LSDptr=,true, { return false; } ); // NDS zu x bilden
      // MSDptr erhöhen und len erniedrigen, so daß len = ceiling(q/intDsize) wird:
      { var uintC qD = ceiling(q,intDsize); // ceiling(q/intDsize)
        // wegen q<=l ist qD = ceiling(q/intDsize) <= ceiling((l+1)/intDsize) = len, also
        // paßt qD ebenso wie len in ein uintC.
        MSDptr = MSDptr mspop (len - qD); // MSDptr um len-qD Digits erhöhen
        len = qD; // len um len-qD erniedrigen
      }
      // LSDptr und len um floor(p/intDsize) erniedrigen:
      { var uintC pD = p/intDsize; // floor(p/intDsize)
        LSDptr = LSDptr lspop pD;
        len -= pD;
      }
      // Jetzt enthält MSDptr/len/LSDptr genau die maßgeblichen Digits.
      if (len==0) return false; // len=0 -> keine Bits abzutesten
      q = ((q-1)%intDsize); // q := intDsize - (intDsize*ceiling(q/intDsize) - q) - 1
      p = p%intDsize; // p := p - intDsize*floor(p/intDsize)
      // Jetzt ist 0 <= q < intDsize, 0 <= p < intDsize.
      // Vom ersten Digit müssen die vorderen intDsize-1-q Bits unberücksichtigt bleiben.
      // Ein AND 2^(q+1)-1 erreicht dies.
      // Vom letzten Digit müssen die hinteren p Bits unberücksichtigt bleiben.
      // Ein AND -2^p erreicht dies.
      if (--len==0) {
        // 1 Digit maßgeblich, wird von beiden Seiten angeschnitten:
        // Ein AND 2^(q+1)-2^p erreicht dies.
        if (!(((uintD)(bitm(q+1)-bit(p)) & mspref(MSDptr,0)) == 0))
          return true;
        else
          return false;
      }
      // mindestens 2 Digits. Teste erst die Randdigits, dann die inneren:
      if (!(((msprefnext(MSDptr) & (uintD)(bitm(q+1)-1)) == 0) &&
            ((lsprefnext(LSDptr) & (uintD)(minus_bit(p))) == 0)
         ) )
        return true;
      len--; // die beiden Randdigits sind jetzt abgezogen.
      if (DS_test_loop(MSDptr,len,LSDptr)) { return true; } else { return false; }
    }

}  // namespace cln
