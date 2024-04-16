// print_rational().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational_io.h"


// Implementation.

#include "cln/output.h"
#include "cln/integer_io.h"
#include "cln/rational.h"
#include "rational/cl_RA.h"

namespace cln {

void print_rational (std::ostream& stream, const cl_print_rational_flags& flags, const cl_RA& z)
{
	var unsigned int base = flags.rational_base;
	if (flags.rational_readably)
		// Radix-Specifier ausgeben:
		switch (base) {
		case 2:
			fprintchar(stream,'#');
			fprintchar(stream,'b');
			break;
		case 8:
			fprintchar(stream,'#');
			fprintchar(stream,'o');
			break;
		case 16:
			fprintchar(stream,'#');
			fprintchar(stream,'x');
			break;
		case 10:
			if (integerp(z)) {
				DeclareType(cl_I,z);
				// Basis 10 bei Integers durch
				// nachgestellten Punkt kennzeichnen:
				print_integer(stream,base,z);
				fprintchar(stream,'.');
				return;
			}
			// fallthrough
		default:
			// Basis in #nR-Schreibweise ausgeben:
			fprintchar(stream,'#');
			print_integer(stream,10,base);
			fprintchar(stream,'r');
			break;
		}
	if (integerp(z)) {
		DeclareType(cl_I,z);
		// Integer in Basis base ausgeben:
		print_integer(stream,base,z);
	} else {
		DeclareType(cl_RT,z);
		// Ratio in Basis base ausgeben; Zähler / Nenner
		print_integer(stream,base,numerator(z));
		fprintchar(stream,'/');
		print_integer(stream,base,denominator(z));
	}
}

}  // namespace cln
