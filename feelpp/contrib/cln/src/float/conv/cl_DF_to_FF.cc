// cl_DF_to_FF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "float/ffloat/cl_FF.h"

namespace cln {

const cl_FF cl_DF_to_FF (const cl_DF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	#if (cl_word_size==64)
	var uint64 mant;
	DF_decode(x, { return cl_FF_0; }, sign=,exp=,mant=);
	// 52-23=29 Bits wegrunden:
	var const int shiftcount = DF_mant_len-FF_mant_len;
	if ( ((mant & bit(shiftcount-1)) ==0) // Bit 28 war 0 -> abrunden
	     || ( ((mant & (bit(shiftcount-1)-1)) ==0) // war 1, Bits 27..0 >0 -> aufrunden
	          // round-to-even
	          && ((mant & bit(shiftcount)) ==0)
	   )    )
	  // abrunden
	  { mant = mant >> shiftcount; }
	  else
	  // aufrunden
	  { mant = mant >> shiftcount;
	    mant = mant+1;
	    if (mant >= bit(FF_mant_len+1))
	      // Überlauf durchs Runden
	      { mant = mant>>1; exp = exp+1; } // Mantisse rechts schieben
	  }
	return encode_FF(sign,exp,mant);
	#else
	var uint32 manthi;
	var uint32 mantlo;
	DF_decode2(x, { return cl_FF_0; }, sign=,exp=,manthi=,mantlo=);
	// 52-23=29 Bits wegrunden:
	var const int shiftcount = DF_mant_len-FF_mant_len;
	manthi = (manthi << (32-shiftcount)) | (mantlo >> shiftcount);
	if ( ((mantlo & bit(shiftcount-1)) ==0) // Bit 28 war 0 -> abrunden
	     || ( ((mantlo & (bit(shiftcount-1)-1)) ==0) // war 1, Bits 27..0 >0 -> aufrunden
	          // round-to-even
	          && ((mantlo & bit(shiftcount)) ==0)
	   )    )
	  // abrunden
	  {}
	  else
	  // aufrunden
	  { manthi = manthi+1;
	    if (manthi >= bit(FF_mant_len+1))
	      // Überlauf durchs Runden
	      { manthi = manthi>>1; exp = exp+1; } // Mantisse rechts schieben
	  }
	return encode_FF(sign,exp,manthi);
	#endif
}

}  // namespace cln
