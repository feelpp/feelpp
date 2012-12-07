// cl_DF_to_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "float/ffloat/cl_FF.h"

namespace cln {

float float_approx (const cl_DF& x)
{
	union { ffloat eksplicit; float machine_float; } u;
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	#if (cl_word_size==64)
	var uint64 mant;
	DF_decode(x, { return 0.0; }, sign=,exp=,mant=);
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
	#else
	var uint32 manthi;
	var uint32 mantlo;
	DF_decode2(x, { return 0.0; }, sign=,exp=,manthi=,mantlo=);
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
	#define mant manthi
	#endif
	if (exp > (sintL)(FF_exp_high-FF_exp_mid))
	  { u.eksplicit = make_FF_word(sign,bit(FF_exp_len)-1,0); } // Infinity
	else
	if (exp < (sintL)(FF_exp_low-FF_exp_mid))
	  { u.eksplicit = make_FF_word(sign,0,0); } // 0.0
	else
	  { u.eksplicit = make_FF_word(sign,exp+FF_exp_mid,mant); }
	return u.machine_float;
}

}  // namespace cln
