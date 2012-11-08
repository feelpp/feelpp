// cl_DF_to_SF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "float/sfloat/cl_SF.h"

namespace cln {

const cl_SF cl_DF_to_SF (const cl_DF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	#if (cl_word_size==64)
	var uint64 mant;
	DF_decode(x, { return SF_0; }, sign=,exp=,mant=);
	// 52-16=36 Bits wegrunden:
	var const int shiftcount = DF_mant_len-SF_mant_len;
	if ( ((mant & bit(shiftcount-1)) ==0) // Bit 35 war 0 -> abrunden
	     || ( ((mant & (bit(shiftcount-1)-1)) ==0) // war 1, Bits 34..0 >0 -> aufrunden
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
	#else
	var uint32 manthi;
	var uint32 mantlo;
	DF_decode2(x, { return SF_0; }, sign=,exp=,manthi=,mantlo=);
	// 52-16=36 Bits wegrunden:
	var const int shiftcount = DF_mant_len-SF_mant_len-32;
	if ( ((manthi & bit(shiftcount-1)) ==0) // Bit 35 war 0 -> abrunden
	     || ( ((manthi & (bit(shiftcount-1)-1)) ==0) // war 1, Bits 34..0 >0 -> aufrunden
	          && (mantlo==0)
	          // round-to-even
	          && ((manthi & bit(shiftcount)) ==0)
	   )    )
	  // abrunden
	  { manthi = manthi >> shiftcount; }
	  else
	  // aufrunden
	  { manthi = manthi >> shiftcount;
	    manthi = manthi+1;
	    if (manthi >= bit(SF_mant_len+1))
	      // Überlauf durchs Runden
	      { manthi = manthi>>1; exp = exp+1; } // Mantisse rechts schieben
	  }
	return encode_SF(sign,exp,manthi);
	#endif
}

}  // namespace cln
