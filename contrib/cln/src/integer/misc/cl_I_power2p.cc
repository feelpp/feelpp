// power2p().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

uintC power2p (const cl_I& x) // x > 0
{
// Methode 1: Wenn ord2(x) = integer_length(x)-1.
// Methode 2: Wenn logand(x,x-1) = 0.
// Methode 3: Wenn das erste Digit /=0 eine Zweierpotenz ist und alle weiteren
//            Digits Null sind.
	if (fixnump(x))
	  { var uintV x_ = FN_to_UV(x);
	    if (!((x_ & (x_-1)) == 0)) return 0; // keine Zweierpotenz
            #if (intVsize>32)
            integerlength64(x_,return); // Zweierpotenz: n = integer_length(x)
            #else
	    integerlength32(x_,return); // Zweierpotenz: n = integer_length(x)
            #endif
	  }
	  else
	  { var const uintD* MSDptr;
	    var uintC len;
	    var const uintD* LSDptr;
	    BN_to_NDS_nocopy(x, MSDptr=,len=,LSDptr=); // normalisierte DS zu x bilden.
	    var uintD msd = mspref(MSDptr,0);
	    if (msd==0) { msshrink(MSDptr); msd = mspref(MSDptr,0); len--; }
	    // len = Anzahl der Digits ab MSDptr, len>0, msd = erstes Digit (/=0)
	    if (!((msd & (msd-1)) == 0)) return 0; // erstes Digit muß Zweierpotenz sein
	    if (DS_test_loop(MSDptr mspop 1,len-1,LSDptr)) return 0; // danach alles Nullen
	   {var uintL msdlen;
	    integerlengthD(msd, msdlen=);
	    return intDsize*(len-1) + msdlen; // integer_length(x) als Ergebnis
	  }}
}

}  // namespace cln
