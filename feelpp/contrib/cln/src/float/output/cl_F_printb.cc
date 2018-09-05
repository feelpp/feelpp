// print_float_binary().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float_io.h"


// Implementation.

#include "cln/float.h"
#include "float/cl_F.h"
#include "cln/integer_io.h"
#include "integer/cl_I.h"

namespace cln {

void print_float_binary (std::ostream& stream, const cl_F& z)
{
// Vorzeichen, Punkt, Mantisse (binär), (Zweiersystem-)Exponent (dezimal)
	cl_idecoded_float m_e_s = integer_decode_float(z);
	var cl_I& m = m_e_s.mantissa;
	var cl_I& s = m_e_s.sign;
	// Vorzeichen ausgeben, falls <0:
	if (eq(s,-1))
		fprintchar(stream,'-');
	// Mantisse binär(!) ausgeben:
	fprintchar(stream,'.');
	print_integer(stream,2,m);
	// Exponent-Marker ausgeben:
	{
		var char exp_marker;
		floattypecase(z
		,	exp_marker = 's';
		,	exp_marker = 'f';
		,	exp_marker = 'd';
		,	exp_marker = 'L';
		);
		fprintchar(stream,exp_marker);
	}
	// Exponenten dezimal ausgeben:
	print_integer(stream,10,cl_I(float_exponent(z)));
}

}  // namespace cln
