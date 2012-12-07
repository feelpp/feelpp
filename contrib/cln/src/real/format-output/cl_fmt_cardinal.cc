// format_cardinal().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "real/format-output/cl_format.h"


// Implementation.

#include <sstream>
#include "cln/integer.h"
#include "cln/integer_io.h"
#include "cln/exception.h"

namespace cln {

static const char * const cl_format_ones [20] = {
	NULL,
	"one",
	"two",
	"three",
	"four",
	"five",
	"six",
	"seven",
	"eight",
	"nine",
	"ten",
	"eleven",
	"twelve",
	"thirteen",
	"fourteen",
	"fifteen",
	"sixteen",
	"seventeen",
	"eighteen",
	"nineteen",
};

// gibt eine ganze Zahl >0, <1000 im Klartext auf englisch auf den stream aus.
// (arg=0 -> gibt nichts aus.)
static void format_small_cardinal (std::ostream& stream, uintL arg)
{
	var uintL hundreds = floor(arg,100);
	var uintL tens_and_ones = arg % 100;
	if (hundreds > 0) {
		fprint(stream,cl_format_ones[hundreds]);
		fprint(stream," hundred");
	}
	if (tens_and_ones > 0) {
		if (hundreds > 0)
			fprint(stream," and ");
		var uintL tens = floor(tens_and_ones,10);
		var uintL ones = tens_and_ones % 10;
		if (tens < 2)
			fprint(stream,cl_format_ones[tens_and_ones]);
		else {
			fprint(stream,cl_format_tens[tens]);
			if (ones > 0) {
				fprintchar(stream,'-');
				fprint(stream,cl_format_ones[ones]);
			}
		}
	}
}

void format_cardinal (std::ostream& stream, const cl_I& argument)
{
	if (zerop(argument))
		fprint(stream,"zero");
	else {
		var cl_I arg = argument;
		if (minusp(arg)) {
			fprint(stream,"minus ");
			arg = -arg;
		}
		// amerikanisch (billion=10^9)
		static const char * const illions[] = {
			"",
			" thousand",
			" million",
			" billion",
			" trillion",
			" quadrillion",
			" quintillion",
			" sextillion",
			" septillion",
			" octillion",
			" nonillion",
			" decillion",
			" undecillion",
			" duodecillion",
			" tredecillion",
			" quattuordecillion",
			" quindecillion",
			" sexdecillion",
			" septendecillion",
			" octodecillion",
			" novemdecillion",
			" vigintillion",
			NULL
		};
		var uintL small_pieces [sizeof(illions)/sizeof(illions[0])];
		// Let the recursion begin.
		var const char * const * illion_ptr = &illions[0];
		var uintL * small_piece_ptr = &small_pieces[0];
		do {
			if (*illion_ptr == NULL) {
				std::ostringstream buf;
				fprint(buf, "format_cardinal: argument too large: ");
				fprint(buf, argument);
				throw runtime_exception(buf.str());
			}
			var cl_I_div_t div = floor2(arg,1000);
			var const cl_I& thousands = div.quotient;
			var uintL small = cl_I_to_UL(div.remainder);
			illion_ptr++;
			*small_piece_ptr++ = small;
			arg = thousands;
		} while (arg > 0);
		// Roll back the recursion.
		var bool first_piece = true;
		do {
			var uintL small = *--small_piece_ptr;
			var const char * illion = *--illion_ptr;
			if (small > 0) {
				if (!first_piece)
					fprint(stream,", ");
				format_small_cardinal(stream,small);
				fprint(stream,illion);
				first_piece = false;
			}
		} until (illion_ptr == &illions[0]);
	}
}

}  // namespace cln
