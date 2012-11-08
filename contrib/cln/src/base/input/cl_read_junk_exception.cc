// read_number_junk_exception().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/number_io.h"


// Implementation.

#include <sstream>
#include "cln/io.h"

namespace cln {

static inline const std::string
read_number_junk_msg (const char * string_rest, const char * string, const char * string_limit)
{
	std::ostringstream buf;
	fprint(buf, "Junk after number: ");
	{ for (const char * ptr = string; ptr != string_rest; ptr++)
		fprintchar(buf, *ptr);
	}
	fprint(buf, "\"");
	{ for (const char * ptr = string_rest; ptr != string_limit; ptr++)
		fprintchar(buf, *ptr);
	}
	fprint(buf, "\"");
	return buf.str();
}

read_number_junk_exception::read_number_junk_exception (const char * string_rest, const char * string, const char * string_limit)
	: read_number_exception(read_number_junk_msg(string_rest, string, string_limit))
{}

}  // namespace cln
