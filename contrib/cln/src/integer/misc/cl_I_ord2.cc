// ord2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

// Methode 1a:
//   Sei n = ord2(x). Dann ist logxor(x,x-1) = 2^n + (2^n-1) = 2^(n+1)-1.
//   Also  (ord2 x) = (1- (integer-length (logxor x (1- x)))) .
// Methode 1b:
//   Sei n = ord2(x). Dann ist logand(x,-x) = 2^n.
//   Also  (ord2 x) = (1- (integer-length (logand x (- x)))) .
// Methode 1c:
//   Sei n = ord2(x). Dann ist lognot(logior(x,-x)) = 2^n-1.
//   Also  (ord2 x) = (integer-length (lognot (logior x (- x)))) .
// Methode 2:
//   Nullbits am Schluß von x abzählen:
//   (ord2 x) = intDsize * Anzahl der Nulldigits am Schluß
//              + Anzahl der Nullbits am Ende des letzten Digits /=0.

uintC ord2 (const cl_I& x) // x /= 0
{
	if (fixnump(x))
	  { var uintV x_ = FN_to_V(x); // x als intVsize-Bit-Zahl
            // This assumes cl_value_len <= intVsize.
            #if (intVsize>32)
	    ord2_64(x_,return);
            #else
            ord2_32(x_,return);
            #endif
	  }
	  else
	  { var uintC bitcount = 0;
	    var const uintD* ptr;
	    BN_to_NDS_nocopy(x, ,,ptr=); // normalisierte DS zu x bilden.
	    while (lspref(ptr,0) == 0) { lsshrink(ptr); bitcount += intDsize; } // Nulldigits abzählen
	    var uintD lsd = lspref(ptr,0); // letztes Digit /=0
	    ord2_D(lsd,bitcount +=); // dessen Nullbits abzählen
	    return bitcount;
          }
}

}  // namespace cln
