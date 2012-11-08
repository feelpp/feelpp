// cl_LF_to_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "float/ffloat/cl_FF.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

float float_approx (const cl_LF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintE exp;
	var uintD* ptr;
	var uintC len;
	LF_decode(x, { return 0.0; }, sign=,exp=,ptr=,len=,);
	// intDsize*len-FF_mant_len-1 Bits der Mantisse wegrunden:
	// erste k := ceiling(FF_mant_len+2,intDsize) Digits nach mant holen:
	#if (intDsize==64)
	var uint64 mant = get_max64_Dptr(FF_mant_len+2,ptr);
	#else
	var uint32 mant = get_max32_Dptr(FF_mant_len+2,ptr);
	#endif
	ptr = ptr mspop ceiling(FF_mant_len+2,intDsize);
	var const int shiftcount = ceiling(FF_mant_len+2,intDsize)*intDsize-(FF_mant_len+1);
	if ( ((mant & bit(shiftcount-1)) ==0) // Bit 7 war 0 -> abrunden
	     || ( ((mant & (bit(shiftcount-1)-1)) ==0) // war 1, Bits 6..0 >0 -> aufrunden
	          && !test_loop_msp(ptr,len-ceiling(FF_mant_len+2,intDsize)) // weitere Bits /=0 -> aufrunden
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
	union { ffloat eksplicit; float machine_float; } u;
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
