// cl_FF_to_double().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "float/dfloat/cl_DF.h"

namespace cln {

double double_approx (const cl_FF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	FF_decode(x, { return 0.0; }, sign=,exp=,mant=);
	// Mantisse um 52-23=29 Nullbits erweitern:
	union { dfloat eksplicit; double machine_double; } u;
	#if (cl_word_size==64)
	if (((sintL)(FF_exp_high-FF_exp_mid) > (sintL)(DF_exp_high-DF_exp_mid))
	    && (exp > (sintL)(DF_exp_high-DF_exp_mid)))
	  { u.eksplicit =
	        ((sint64)sign & bit(63))
	      | ((uint64)(bit(DF_exp_len)-1) << DF_mant_len); // Infinity
	  }
	else
	if (((sintL)(FF_exp_low-FF_exp_mid) < (sintL)(DF_exp_low-DF_exp_mid))
	    && (exp < (sintL)(DF_exp_low-DF_exp_mid)))
	  { u.eksplicit = ((sint64)sign & bit(63)); } // 0.0
	else
	  { u.eksplicit =
	        ((sint64)sign & bit(63))                  /* Vorzeichen */
	      | ((uint64)(exp+DF_exp_mid) << DF_mant_len) /* Exponent   */
	      | (((uint64)mant<<(DF_mant_len-FF_mant_len)) & (bit(DF_mant_len)-1)); /* Mantisse   */
	  }
	#else
	if (((sintL)(FF_exp_high-FF_exp_mid) > (sintL)(DF_exp_high-DF_exp_mid))
	    && (exp > (sintL)(DF_exp_high-DF_exp_mid)))
	  { u.eksplicit.semhi =
	        ((sint32)sign & bit(31))
	      | ((uint32)(bit(DF_exp_len)-1) << (DF_mant_len-32)); // Infinity
	    u.eksplicit.mlo = 0;
	  }
	else
	if (((sintL)(FF_exp_low-FF_exp_mid) < (sintL)(DF_exp_low-DF_exp_mid))
	    && (exp < (sintL)(DF_exp_low-DF_exp_mid)))
	  { u.eksplicit.semhi = ((sint32)sign & bit(31)); // 0.0
	    u.eksplicit.mlo = 0;
	  }
	else
	  { u.eksplicit.semhi =
	        ((sint32)sign & bit(31))                       /* Vorzeichen */
	      | ((uint32)(exp+DF_exp_mid) << (DF_mant_len-32)) /* Exponent   */
	      | (((uint32)mant>>(32-(DF_mant_len-FF_mant_len))) & (bit(DF_mant_len-32)-1)); /* Mantisse   */
	    u.eksplicit.mlo = mant<<(DF_mant_len-FF_mant_len);
	  }
	#endif
	return u.machine_double;
}

}  // namespace cln
