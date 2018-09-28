// cl_LF_to_SF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "float/sfloat/cl_SF.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_SF cl_LF_to_SF (const cl_LF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintE exp;
	var uintD* ptr;
	var uintC len;
	LF_decode(x, { return SF_0; }, sign=,exp=,ptr=,len=,);
	// intDsize*len-SF_mant_len-1 Bits der Mantisse wegrunden:
	// erste k := ceiling(SF_mant_len+2,intDsize) Digits nach mant holen:
	#if (intDsize==64)
	var uint64 mant = get_max64_Dptr(SF_mant_len+2,ptr);
	#else
	var uint32 mant = get_max32_Dptr(SF_mant_len+2,ptr);
	#endif
	ptr = ptr mspop ceiling(SF_mant_len+2,intDsize);
	var const int shiftcount = ceiling(SF_mant_len+2,intDsize)*intDsize-(SF_mant_len+1);
	if ( ((mant & bit(shiftcount-1)) ==0) // Bit 14 war 0 -> abrunden
	     || ( ((mant & (bit(shiftcount-1)-1)) ==0) // war 1, Bits 13..0 >0 -> aufrunden
	          && !test_loop_msp(ptr,len-ceiling(SF_mant_len+2,intDsize)) // weitere Bits /=0 -> aufrunden
	          // round-to-even
	          && ((mant & bit(shiftcount)) ==0)
	   )    )
	  // abrunden
	  { mant = mant >> shiftcount; }
	  else
	  // aufrunden
	  { mant = mant >> shiftcount;
	    mant = mant+1;
	    if (mant >= bit(SF_mant_len+1))
	      // Überlauf durchs Runden
	      { mant = mant>>1; exp = exp+1; } // Mantisse rechts schieben
	  }
	return encode_SF(sign,exp,mant);
}

}  // namespace cln
