// most_positive_float().

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

static inline const cl_LF most_positive_LF (uintC len)
{
	var Lfloat erg = allocate_lfloat(len,LF_exp_high,0);
	fill_loop_up(&TheLfloat(erg)->data[0],len,~(uintD)0);
	return erg;
}

const cl_F most_positive_float (float_format_t f)
{
	// Exponent so groß wie möglich, Mantisse 1...1, Vorzeichen +.

	static const cl_SF most_positive_SF =
		make_SF(0,SF_exp_high,bit(SF_mant_len+1)-1);

	static const cl_FF most_positive_FF =
		encode_FF(0,FF_exp_high-FF_exp_mid,bit(FF_mant_len+1)-1);

	static const cl_DF most_positive_DF =
#if (cl_word_size==64)
		encode_DF(0,DF_exp_high-DF_exp_mid,bit(DF_mant_len+1)-1);
#else
		encode_DF(0,DF_exp_high-DF_exp_mid,bit(DF_mant_len-32+1)-1,bitm(32)-1);
#endif

	floatformatcase((uintC)f
	,	return most_positive_SF;
	,	return most_positive_FF;
	,	return most_positive_DF;
	,	return most_positive_LF(len);
	);
}

}  // namespace cln

