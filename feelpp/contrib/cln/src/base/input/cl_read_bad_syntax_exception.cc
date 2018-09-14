// read_number_bad_syntax_exception().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/number_io.h"


// Implementation.

#include <sstream>
#include "cln/io.h"

namespace cln {

static inline const std::string
read_number_bad_syntax_msg (const char * string, const char * string_limit)
{
	std::ostringstream buf;
	fprint(buf, "Illegal number syntax: \"");
	for (const char * ptr = string; ptr != string_limit; ptr++)
		fprintchar(buf, *ptr);
	fprint(buf, "\"");
	return buf.str();
}

read_number_bad_syntax_exception::read_number_bad_syntax_exception (const char * string, const char * string_limit)
	: read_number_exception(read_number_bad_syntax_msg(string, string_limit))
{}

}  // namespace cln
