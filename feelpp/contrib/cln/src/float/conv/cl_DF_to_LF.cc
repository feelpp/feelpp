// cl_DF_to_LF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_LF cl_DF_to_LF (const cl_DF& x, uintC len)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	#if (cl_word_size==64)
	var uint64 mant;
	DF_decode(x, { return encode_LF0(len); }, sign=,exp=(sintL),mant=);
	#else
	var uint32 manthi;
	var uint32 mantlo;
	DF_decode2(x, { return encode_LF0(len); }, sign=,exp=(sintL),manthi=,mantlo=);
	#endif
	// Long-Float allozieren,
	// Mantisse mit intDsize*len-DF_mant_len-1 Nullbits auffüllen:
	var Lfloat y = allocate_lfloat(len,exp+LF_exp_mid,sign);
	var uintD* ptr = arrayMSDptr(TheLfloat(y)->data,len);
	// erste k := ceiling(DF_mant_len+1,intDsize) Digits mit mant füllen:
	var const int shiftcount = (ceiling(DF_mant_len+1,intDsize)*intDsize-(DF_mant_len+1));
	#if (cl_word_size==64)
	mant = mant<<shiftcount;
	set_max64_Dptr(DF_mant_len+1,ptr,mant);
	#else
	manthi = (manthi<<shiftcount) | (mantlo>>(32-shiftcount));
	mantlo = mantlo<<shiftcount;
	set_max32_Dptr(DF_mant_len+1-32,ptr,manthi);
	set_32_Dptr(ptr mspop ceiling(DF_mant_len+1-32,intDsize),mantlo);
	#endif
	clear_loop_msp(ptr mspop ceiling(DF_mant_len+1,intDsize),len-ceiling(DF_mant_len+1,intDsize));
	return y;
}

}  // namespace cln
