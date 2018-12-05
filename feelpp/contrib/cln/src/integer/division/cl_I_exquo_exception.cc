// exquo_exception().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "cln/io.h"
#include "cln/integer_io.h"
#include <sstream>

namespace cln {

static inline const std::string
exquo_error_msg (const cl_I& x, const cl_I& y)
{
	std::ostringstream buf;
	fprint(buf, "Quotient ");
	fprint(buf, x);
	fprint(buf, " / ");
	fprint(buf, y);
	fprint(buf, " is not an integer.");
	return buf.str();
}

exquo_exception::exquo_exception (const cl_I& x, const cl_I& y)
	: runtime_exception(exquo_error_msg(x,y))
{}

}  // namespace cln
