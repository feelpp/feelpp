// lognot().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "integer/bitwise/cl_I_log.h"

namespace cln {

const cl_I lognot (const cl_I& x)
    { if (fixnump(x)) // Fixnum -> ganz einfach:
        { // bitweise als Fixnum zurück
          return cl_I_from_word(x.word ^ cl_combine(0,~(cl_uint)0));
        }
        else
        // Bignum:
        { CL_ALLOCA_STACK;
          var uintD* MSDptr;
          var uintC n;
          BN_to_NDS(x, MSDptr=,n=,); // NDS zu x bilden
          // Es ist n>=bn_minlength,
          // und die ersten intDsize+1 Bit sind nicht alle gleich.
          not_loop_msp(MSDptr,n); // mit NOT komplementieren,
                          // wegen n>0 wird auch das Vorzeichenbit umgedreht
          // MSDptr/n/LSDptr ist immer noch eine NDS, da n>=bn_minlength
          // und die ersten intDsize+1 Bit nicht alle gleich sind.
          return NDS_to_I(MSDptr,n); // Ergebnis als Integer
    }   }

}  // namespace cln
