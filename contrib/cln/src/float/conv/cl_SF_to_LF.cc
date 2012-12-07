// cl_SF_to_LF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_LF cl_SF_to_LF (const cl_SF& x, uintC len)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	#if (intDsize==64)
	var uint64 mant;
	#else
	var uint32 mant;
	#endif
	SF_decode(x, { return encode_LF0(len); }, sign=,exp=(sintL),mant=);
	// Long-Float allozieren,
	// Mantisse mit intDsize*len-SF_mant_len-1 Nullbits auffüllen:
	var Lfloat y = allocate_lfloat(len,exp+LF_exp_mid,sign);
	var uintD* ptr = arrayMSDptr(TheLfloat(y)->data,len);
	// erste k := ceiling(SF_mant_len+1,intDsize) Digits mit mant füllen:
	mant = mant << (ceiling(SF_mant_len+1,intDsize)*intDsize-(SF_mant_len+1));
	#if (intDsize==64)
	set_max64_Dptr(SF_mant_len+1,ptr,mant);
	#else
	set_max32_Dptr(SF_mant_len+1,ptr,mant);
	#endif
	clear_loop_msp(ptr mspop ceiling(SF_mant_len+1,intDsize),len-ceiling(SF_mant_len+1,intDsize));
	return y;
}

}  // namespace cln
