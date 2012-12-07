// print_integer().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer_io.h"


// Implementation.

#include "cln/output.h"

namespace cln {

void print_integer (std::ostream& stream, const cl_print_rational_flags& flags, const cl_I& z)
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
			// Basis 10 bei Integers durch
			// nachgestellten Punkt kennzeichnen:
			print_integer(stream,base,z);
			fprintchar(stream,'.');
			return;
		default:
			// Basis in #nR-Schreibweise ausgeben:
			fprintchar(stream,'#');
			print_integer(stream,10,base);
			fprintchar(stream,'r');
			break;
		}
	// Integer in Basis base ausgeben:
	print_integer(stream,base,z);
}

}  // namespace cln
