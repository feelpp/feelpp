// format_ordinal().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "real/format-output/cl_format.h"


// Implementation.

#include "cln/integer.h"

namespace cln {

static const char * const cl_format_ordinal_ones [20] = {
	NULL,
	"first",
	"second",
	"third",
	"fourth",
	"fifth",
	"sixth",
	"seventh",
	"eighth",
	"ninth",
	"tenth",
	"eleventh",
	"twelfth",
	"thirteenth",
	"fourteenth",
	"fifteenth",
	"sixteenth",
	"seventeenth",
	"eighteenth",
	"nineteenth",
};

static const char * const cl_format_ordinal_tens [10] = {
	NULL,
	"tenth",
	"twentieth",
	"thirtieth",
	"fortieth",
	"fiftieth",
	"sixtieth",
	"seventieth",
	"eightieth",
	"ninetieth",
};

void format_ordinal (std::ostream& stream, const cl_I& argument)
{
	if (zerop(argument))
		fprint(stream,"zeroth");
	else {
		var cl_I arg = argument;
		if (minusp(arg)) {
			fprint(stream,"minus ");
			arg = -arg;
		}
		var cl_I_div_t div = floor2(arg,100);
		var const cl_I& hundreds = div.quotient;
		var uintL tens_and_ones = cl_I_to_UL(div.remainder);
		if (hundreds > 0)
			format_cardinal(stream,hundreds*100);
		if (tens_and_ones == 0)
			fprint(stream,"th");
		else {
			var uintL tens = floor(tens_and_ones,10);
			var uintL ones = tens_and_ones % 10;
			if (hundreds > 0)
				fprintchar(stream,' ');
			if (tens < 2)
				fprint(stream,cl_format_ordinal_ones[tens_and_ones]);
			elif (ones == 0)
				fprint(stream,cl_format_ordinal_tens[tens]);
			else {
				fprint(stream,cl_format_tens[tens]);
				fprintchar(stream,'-');
				fprint(stream,cl_format_ordinal_ones[ones]);
			}
		}
	}
}

}  // namespace cln
