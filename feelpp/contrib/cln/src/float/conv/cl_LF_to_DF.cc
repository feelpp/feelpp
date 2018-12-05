// cl_LF_to_DF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "float/dfloat/cl_DF.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_DF cl_LF_to_DF (const cl_LF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintE exp;
	var uintD* ptr;
	var uintC len;
	LF_decode(x, { return cl_DF_0; }, sign=,exp=,ptr=,len=,);
	// intDsize*len-DF_mant_len-1 Bits der Mantisse wegrunden:
	// erste k := ceiling(DF_mant_len+2,intDsize) Digits nach manthi,mantlo holen:
	var const int shiftcount = ceiling(DF_mant_len+2,intDsize)*intDsize-(DF_mant_len+1);
	#if (cl_word_size==64)
	var uint64 mant = get_max64_Dptr(DF_mant_len+2,ptr);
	ptr = ptr mspop ceiling(DF_mant_len+2,intDsize);
	if ( ((mant & bit(shiftcount-1)) ==0) // Bit 10 war 0 -> abrunden
	     || ( ((mant & (bit(shiftcount-1)-1)) ==0) // war 1, Bits 9..0 >0 -> aufrunden
	          && !test_loop_msp(ptr,len-ceiling(DF_mant_len+2,intDsize)) // weitere Bits /=0 -> aufrunden
	          // round-to-even
	          && ((mant & bit(shiftcount)) ==0)
	   )    )
	  // abrunden
	  { mant = mant >> shiftcount; }
	  else
	  // aufrunden
	  { mant = mant >> shiftcount;
	    mant = mant+1;
	    if (mant >= bit(DF_mant_len+1))
	      // Überlauf durchs Runden
	      { mant = mant>>1; exp = exp+1; } // Mantisse rechts schieben
	  }
	return encode_DF(sign,exp,mant);
	#else
	var uint32 manthi = get_max32_Dptr(DF_mant_len+2-32,ptr);
	var uint32 mantlo = get_32_Dptr(ptr mspop ceiling(DF_mant_len+2-32,intDsize));
	ptr = ptr mspop ceiling(DF_mant_len+2,intDsize);
	if ( ((mantlo & bit(shiftcount-1)) ==0) // Bit 10 war 0 -> abrunden
	     || ( ((mantlo & (bit(shiftcount-1)-1)) ==0) // war 1, Bits 9..0 >0 -> aufrunden
	          && !test_loop_msp(ptr,len-ceiling(DF_mant_len+2,intDsize)) // weitere Bits /=0 -> aufrunden
	          // round-to-even
	          && ((mantlo & bit(shiftcount)) ==0)
	   )    )
	  // abrunden
	  { mantlo = (manthi << (32-shiftcount)) | (mantlo >> shiftcount);
	    manthi = manthi >> shiftcount;
	  }
	  else
	  // aufrunden
	  { mantlo = (manthi << (32-shiftcount)) | (mantlo >> shiftcount);
	    manthi = manthi >> shiftcount;
	    mantlo = mantlo+1;
	    if (mantlo==0)
	      { manthi = manthi+1;
	        if (manthi >= bit(DF_mant_len+1-32))
	          // Überlauf durchs Runden
	          { manthi = manthi>>1; exp = exp+1; } // Mantisse rechts schieben
	  }   }
	return encode_DF(sign,exp,manthi,mantlo);
	#endif
}

}  // namespace cln
