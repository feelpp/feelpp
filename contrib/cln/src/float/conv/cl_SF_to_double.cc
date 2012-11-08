// cl_SF_to_double().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "float/dfloat/cl_DF.h"

namespace cln {

double double_approx (const cl_SF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	SF_decode(x, { return 0.0; }, sign=,exp=,mant=);
	// Mantisse um 52-16=36 Nullbits erweitern:
	union { dfloat eksplicit; double machine_double; } u;
	#if (cl_word_size==64)
	if (((sintL)(SF_exp_high-SF_exp_mid) > (sintL)(DF_exp_high-DF_exp_mid))
	    && (exp > (sintL)(DF_exp_high-DF_exp_mid)))
	  { u.eksplicit =
	        ((sint64)sign & bit(63))
	      | ((uint64)(bit(DF_exp_len)-1) << DF_mant_len); // Infinity
	  }
	else
	if (((sintL)(SF_exp_low-SF_exp_mid) < (sintL)(DF_exp_low-DF_exp_mid))
	    && (exp < (sintL)(DF_exp_low-DF_exp_mid)))
	  { u.eksplicit = ((sint64)sign & bit(63)); } // 0.0
	else
	  { u.eksplicit =
	        ((sint64)sign & bit(63))                  /* Vorzeichen */
	      | ((uint64)(exp+DF_exp_mid) << DF_mant_len) /* Exponent   */
	      | (((uint64)mant<<(DF_mant_len-SF_mant_len)) & (bit(DF_mant_len)-1)); /* Mantisse   */
	  }
	#else
	if (((sintL)(SF_exp_high-SF_exp_mid) > (sintL)(DF_exp_high-DF_exp_mid))
	    && (exp > (sintL)(DF_exp_high-DF_exp_mid)))
	  { u.eksplicit.semhi =
	        ((sint32)sign & bit(31))
	      | ((uint32)(bit(DF_exp_len)-1) << (DF_mant_len-32)); // Infinity
	    u.eksplicit.mlo = 0;
	  }
	else
	if (((sintL)(SF_exp_low-SF_exp_mid) < (sintL)(DF_exp_low-DF_exp_mid))
	    && (exp < (sintL)(DF_exp_low-DF_exp_mid)))
	  { u.eksplicit.semhi = ((sint32)sign & bit(31)); // 0.0
	    u.eksplicit.mlo = 0;
	  }
	else
	  { u.eksplicit.semhi =
	        ((sint32)sign & bit(31))                       /* Vorzeichen */
	      | ((uint32)(exp+DF_exp_mid) << (DF_mant_len-32)) /* Exponent   */
	      | (((uint32)mant<<(DF_mant_len-SF_mant_len-32)) & (bit(DF_mant_len-32)-1)); /* Mantisse   */
	    u.eksplicit.mlo = 0;
	  }
	#endif
	return u.machine_double;
}

}  // namespace cln
