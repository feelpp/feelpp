// float_negative_epsilon().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"

// Implementation.

#include "float/cl_F.h"
#include "float/sfloat/cl_SF.h"
#include "float/ffloat/cl_FF.h"
#include "float/dfloat/cl_DF.h"
#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"

namespace cln {

static inline const cl_LF LF_negative_epsilon (uintC len)
{
	var Lfloat erg = allocate_lfloat(len,LF_exp_mid-intDsize*len,0);
	var uintD* ptr = &TheLfloat(erg)->data[0];
	#if CL_DS_BIG_ENDIAN_P
	  *ptr++ = bit(intDsize-1);
	  ptr = clear_loop_up(ptr,len-2);
	  *ptr = bit(0);
	#else
	  *ptr++ = bit(0);
	  ptr = clear_loop_up(ptr,len-2);
	  *ptr = bit(intDsize-1);
	#endif
	return erg;
}

const cl_F float_negative_epsilon (float_format_t f)
{
	// Bei Floats mit d Bits (incl. Hiddem Bit, also d = ?F_mant_len+1)
	// ist ?F_negative_epsilon = 2^(-d-1)*(1+2^(1-d)),
	// d.h. Mantisse 10...01, Vorzeichen +.

	static const cl_SF SF_negative_epsilon =
		make_SF(0,SF_exp_mid-SF_mant_len-1,bit(SF_mant_len)+1);

	static const cl_FF FF_negative_epsilon =
		encode_FF(0,-FF_mant_len-1,bit(FF_mant_len)+1);

	static const cl_DF DF_negative_epsilon =
	#if (cl_word_size==64)
		encode_DF(0,-DF_mant_len-1,bit(DF_mant_len)+1);
	#else
		encode_DF(0,-DF_mant_len-1,bit(DF_mant_len-32),1);
	#endif

	floatformatcase((uintC)f
	,	return SF_negative_epsilon;
	,	return FF_negative_epsilon;
	,	return DF_negative_epsilon;
	,	return LF_negative_epsilon(len);
	);
}

}  // namespace cln

