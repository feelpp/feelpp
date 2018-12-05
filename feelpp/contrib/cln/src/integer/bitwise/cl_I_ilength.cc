// integer_length().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

uintC integer_length (const cl_I& x)
{
	if (fixnump(x))
	  { var uintL bitcount = 0;
	    var uintV x_ = FN_to_V(x); // x als intVsize-Bit-Zahl
	    if (FN_V_minusp(x,(sintV)x_)) { x_ = ~ x_; } // falls <0, komplementieren
	    if (!(x_==0)) {
              #if (intVsize>32)
              integerlength64(x_,bitcount=);
              #else
              integerlength32(x_,bitcount=);
              #endif
            }
	    return bitcount; // 0 <= bitcount < intVsize.
	  }
          else
          { var const uintD* MSDptr;
            var uintC len;
            BN_to_NDS_nocopy(x, MSDptr=,len=,); // normalisierte DS zu x bilden.
            var uintC bitcount = intDsize*(len-1); // Anzahl Digits mal intDsize
            // MSDigit nehmen, testen, welches das höchste Bit ist, das vom
            // Vorzeichenbit abweicht:
            var uintD msd = mspref(MSDptr,0); // MSDigit
            if ((sintD)msd < 0) { msd = ~msd; } // falls negativ, invertieren
            // Position des höchsten Bits in msd suchen und entsprechend bit_count
            // erhöhen (um höchstens intDsize-1):
            if (!(msd == 0)) { integerlengthD(msd, bitcount += ); }
            return bitcount; // 0 <= bitcount < intDsize*2^intCsize.
          }
}

}  // namespace cln
