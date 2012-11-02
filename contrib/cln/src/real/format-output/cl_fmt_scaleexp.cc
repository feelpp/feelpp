// format_scale_exponent().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "real/format-output/cl_format.h"


// Implementation.

#include "cln/real.h"
#include "cln/integer.h"
#include "cln/float.h"
#include "float/cl_F.h"
#include "float/sfloat/cl_SF.h"
#include "float/ffloat/cl_FF.h"
#include "float/dfloat/cl_DF.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

// NOTE: This may introduce roundoff-errors, through the use of *, /, expt.
// But this doesn't matter since format_float_to_string() works with
// exact integers, starting with integer_decode_float().

// For a floating point format f, five characteristic numbers:
struct float_format_params {
	cl_F zero;	// cl_float(0,f)
	cl_F one;	// cl_float(1,f)
	cl_F ten;	// cl_float(10,f)
	cl_F tenth;	// cl_float(1/10,f)
	cl_F lg2;	// log(10,2), as needed (max. 32 bits)
// Constructor:
	float_format_params (cl_F a, cl_F b, cl_F c, cl_F d, cl_F e)
		: zero(a), one(b), ten(c), tenth(d), lg2(e) {}
};

static const float_format_params get_float_params (const cl_F& arg)
{
	static const cl_RA tenth = (cl_RA)"1/10";
	static const cl_SF SF_zero = cl_RA_to_SF(0);
	static const cl_SF SF_one = cl_RA_to_SF(1);
	static const cl_SF SF_ten = cl_RA_to_SF(10);
	static const cl_SF SF_tenth = cl_RA_to_SF(tenth);
	static const cl_FF FF_zero = cl_RA_to_FF(0);
	static const cl_FF FF_one = cl_RA_to_FF(1);
	static const cl_FF FF_ten = cl_RA_to_FF(10);
	static const cl_FF FF_tenth = cl_RA_to_FF(tenth);
	static const cl_DF DF_zero = cl_RA_to_DF(0);
	static const cl_DF DF_one = cl_RA_to_DF(1);
	static const cl_DF DF_ten = cl_RA_to_DF(10);
	static const cl_DF DF_tenth = cl_RA_to_DF(tenth);
	static const cl_SF SF_lg2 = (cl_SF)"0.30103";
	static const cl_DF DF_lg2 = (cl_DF)"0.30102999566";

	floattypecase(arg
	,	return float_format_params(SF_zero,SF_one,SF_ten,SF_tenth,SF_lg2);
	,	return float_format_params(FF_zero,FF_one,FF_ten,FF_tenth,SF_lg2);
	,	return float_format_params(DF_zero,DF_one,DF_ten,DF_tenth,SF_lg2);
	,	var uintC len = TheLfloat(arg)->len;
		return float_format_params(
			cl_I_to_LF(0,len),
			cl_I_to_LF(1,len),
			cl_I_to_LF(10,len),
			cl_RA_to_LF(tenth,len),
			DF_lg2 // lg2 wird mit 32 Bit Genauigkeit gebraucht
		       );
	);
}

const decoded_float format_scale_exponent (const cl_F& arg)
{
	// Get float format parameters.
	var const float_format_params params = get_float_params(arg);
	var const cl_F& zero = params.zero;
	var const cl_F& one = params.one;
	var const cl_F& ten = params.ten;
	var const cl_F& tenth = params.tenth;
	var const cl_F& lg2 = params.lg2;
	// Decode arg.
	if (zerop(arg))
		return decoded_float(zero,0,one);
	var cl_F abs_arg = abs(arg);
	var decoded_float decoded = decode_float(abs_arg);
	var cl_I& expon = decoded.exponent;
	var cl_I expon10a = truncate1(expon*lg2); // nicht round, um Überlauf zu vermeiden
	var cl_F signif10a = abs_arg / expt(ten,expon10a);
	// Maybe need to increment expon10.
	var cl_I expon10b = expon10a;
	var cl_F signif10b = signif10a;
	{
		var cl_F tenpow = ten;
		until (signif10b < one) {
			expon10b = expon10b + 1;
			signif10b = signif10a / tenpow;
			tenpow = tenpow * ten;
		}
	}
	// Maybe need to decrement expon10.
	var cl_I expon10c = expon10b;
	var cl_F signif10c = signif10b;
	{
		var cl_F tenpow = ten;
		until (signif10c >= tenth) {
			expon10c = expon10c - 1;
			signif10c = signif10b * tenpow;
			tenpow = tenpow * ten;
		}
	}
	return decoded_float(signif10c,expon10c,float_sign(arg));
}

}  // namespace cln

