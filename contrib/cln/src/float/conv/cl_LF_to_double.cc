// cl_LF_to_double().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "float/dfloat/cl_DF.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

double double_approx (const cl_LF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintE exp;
	var uintD* ptr;
	var uintC len;
	LF_decode(x, { return 0.0; }, sign=,exp=,ptr=,len=,);
	// intDsize*len-DF_mant_len-1 Bits der Mantisse wegrunden:
	// erste k := ceiling(DF_mant_len+2,intDsize) Digits nach manthi,mantlo holen:
	var const int shiftcount = ceiling(DF_mant_len+2,intDsize)*intDsize-(DF_mant_len+1);
	union { dfloat eksplicit; double machine_double; } u;
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
	if (exp > (sintL)(DF_exp_high-DF_exp_mid))
	  { u.eksplicit =
	        ((sint64)sign & bit(63))
	      | ((uint64)(bit(DF_exp_len)-1) << DF_mant_len); // Infinity
	  }
	else
	if (exp < (sintL)(DF_exp_low-DF_exp_mid))
	  { u.eksplicit = ((sint64)sign & bit(63)); } // 0.0
	else
	  { u.eksplicit =
	        ((sint64)sign & bit(63))                  /* Vorzeichen */
	      | ((uint64)(exp+DF_exp_mid) << DF_mant_len) /* Exponent   */
	      | ((uint64)mant & (bit(DF_mant_len)-1));    /* Mantisse   */
	  }
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
	if (exp > (sintL)(DF_exp_high-DF_exp_mid))
	  { u.eksplicit.semhi =
	        ((sint32)sign & bit(31))
	      | ((uint32)(bit(DF_exp_len)-1) << (DF_mant_len-32)); // Infinity
	    u.eksplicit.mlo = 0;
	  }
	else
	if (exp < (sintL)(DF_exp_low-DF_exp_mid))
	  { u.eksplicit.semhi = ((sint32)sign & bit(31)); // 0.0
	    u.eksplicit.mlo = 0;
	  }
	else
	  { u.eksplicit.semhi =
	        ((sint32)sign & bit(31))                       /* Vorzeichen */
	      | ((uint32)(exp+DF_exp_mid) << (DF_mant_len-32)) /* Exponent   */
	      | ((uint32)manthi & (bit(DF_mant_len-32)-1));    /* Mantisse   */
	    u.eksplicit.mlo = mantlo;
	  }
	#endif
	return u.machine_double;
}

}  // namespace cln
