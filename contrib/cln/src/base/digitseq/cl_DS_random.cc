// Digit sequence level random number generator.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/random/cl_random_impl.h"


// Implementation.

#include "cln/random.h"
#include "base/digitseq/cl_DS.h"
#include "base/cl_low.h"

namespace cln {

void random_UDS (random_state& randomstate, uintD* ptr, uintC len)
{
	var uintC count;
	#if (intDsize==64)
	dotimesC(count,len,
	  { mspref(ptr,0) = random64(randomstate); ptr = ptr mspop 1; });
	#else // (intDsize<=32)
	dotimesC(count,floor(len,32/intDsize),
	  { var uint32 next = random32(randomstate); // weitere 32/intDsize Digits besorgen
	    set_32_Dptr(ptr,next); ptr = ptr mspop 32/intDsize;
	  });
	len = len % (32/intDsize); // Anzahl noch fehlender Digits
	if (len>0)
	  { var uint32 next = random32(randomstate); // weitere 32/intDsize Digits besorgen
	    set_max32_Dptr(intDsize*len,ptr,next);
	  }
	#endif
}

}  // namespace cln
