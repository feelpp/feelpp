// square().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "base/cl_low.h"

namespace cln {

const cl_I square (const cl_I& x)
{
  // Methode:
  // x Fixnum -> direkt multiplizieren
  // sonst: zu DS machen, multiplizieren.
      if (fixnump(x))
        { var sintV x_ = FN_to_V(x);
          #if (cl_value_len > 32)
          // nur falls x ein Integer mit höchstens 32 Bit ist:
          if ((uintV)((sintV)sign_of(x_) ^ x_) < bit(31))
          #endif
            {
              // Werte direkt multiplizieren:
              var uint32 hi;
              var uint32 lo;
              mulu32((uint32)x_,(uint32)x_,hi=,lo=); // erst unsigned multiplizieren
              if (x_ < 0) { hi -= 2*(uint32)x_; } // dann Korrektur für Vorzeichen
              return L2_to_I(hi,lo);
            }
        }
      CL_ALLOCA_STACK;
      var const uintD* xMSDptr;
      var uintC xlen;
      var const uintD* xLSDptr;
      I_to_NDS_nocopy(x, xMSDptr = , xlen = , xLSDptr = , false,);
      var uintD* ergMSDptr;
      var uintC erglen = 2*xlen;
      var uintD* ergLSDptr;
      num_stack_alloc(erglen,ergMSDptr=,ergLSDptr=);
      var uintC len = xlen;
      var uintD MSD = mspref(xMSDptr,0);
      if (MSD==0)
        { mspref(ergMSDptr,0) = 0; mspref(ergMSDptr,1) = 0; len--; }
      cl_UDS_mul_square(xLSDptr,len,ergLSDptr);
      if ((sintD)MSD < 0)
        { subfrom_loop_lsp(xLSDptr,ergLSDptr lspop xlen,xlen);
          subfrom_loop_lsp(xLSDptr,ergLSDptr lspop xlen,xlen);
        }
      return DS_to_I(ergMSDptr,erglen);
}
// Bit complexity (x of length N): O(M(N)).

}  // namespace cln
