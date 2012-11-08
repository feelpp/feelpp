// print_rational().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational_io.h"


// Implementation.

#include "cln/integer_io.h"
#include "cln/rational.h"
#include "rational/cl_RA.h"

namespace cln {

void print_rational (std::ostream& stream, unsigned int base, const cl_RA& z)
{
	if (integerp(z)) {
		DeclareType(cl_I,z);
		print_integer(stream,base,z);
	} else {
		DeclareType(cl_RT,z);
		var const cl_I& num = numerator(z);
		var const cl_I& den = denominator(z);
		// Der Zähler trägt das Vorzeichen.
		print_integer(stream,base,num); // Zähler ausgeben
		fprintchar(stream,'/');
		print_integer(stream,base,den); // Nenner ausgeben
	}
}

}  // namespace cln
