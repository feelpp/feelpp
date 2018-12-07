// cl_notreached_exception().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/cl_macros.h"


// Implementation.

#include "cln/io.h"
#include <sstream>

namespace cln {

static inline const std::string
notreached_error_msg (const char* filename, int lineno)
{
	std::ostringstream buf;
	fprint(buf, "Internal error: statement in file ");
	fprint(buf, filename);
	fprint(buf, ", line ");
	fprintdecimal(buf, lineno);
	fprint(buf, " has been reached!!\n");
	fprint(buf, "Please send the authors of the program a description how you produced this error!");
	return buf.str();
}

notreached_exception::notreached_exception (const char* filename, int lineno)
	: runtime_exception(notreached_error_msg(filename, lineno))
{}

}  // namespace cln
