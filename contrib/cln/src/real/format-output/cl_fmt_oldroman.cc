// format_old_roman().

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

void format_old_roman (std::ostream& stream, const cl_I& arg)
{
	if (!(0 < arg && arg < 5000)) {
		std::ostringstream buf;
		fprint(buf, "format_old_roman: argument should be in the range 1 - 4999, not ");
		fprint(buf, arg);
		fprint(buf, "\n");
		throw runtime_exception(buf.str());
	}
	var uintL value = cl_I_to_UL(arg);
	struct roman { char symbol; uintL value; };
	static const roman scale[7] = {
		{ 'I',    1 },
		{ 'V',    5 },
		{ 'X',   10 },
		{ 'L',   50 },
		{ 'C',  100 },
		{ 'D',  500 },
		{ 'M', 1000 },
	};
	for (int i = 6; value > 0 /* && i >= 0 */ ; i--) {
		var const roman * p = &scale[i];
		var uintL multiplicity = floor(value,p->value);
		value = value % p->value;
		while (multiplicity > 0) {
			fprintchar(stream,p->symbol);
			multiplicity--;
		}
	}
}

}  // namespace cln
