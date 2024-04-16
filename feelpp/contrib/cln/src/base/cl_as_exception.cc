// as_exception().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/exception.h"


// Implementation.

#include "cln/io.h"
#include "base/cl_N.h"
#include <sstream>

namespace cln {

static inline const std::string
as_error_msg (const cl_number& obj, const char * typestring, const char * filename, int line)
{
	std::ostringstream buf;
	fprint(buf, "Type assertion failed: in file ");
	fprint(buf, filename);
	fprint(buf, ", line ");
	fprintdecimal(buf, line);
	fprint(buf, ", not ");
	fprint(buf, typestring);
	fprint(buf, ": ");
#if 0 // This brings in a dependency from the complex and float printer and all the float stuff.
	fprint(buf, obj);
#else
	fprint(buf, "@0x");
	fprinthexadecimal(buf, (uintptr_t)(void*)&obj);
	fprint(buf, ": 0x");
	fprinthexadecimal(buf, (uintptr_t)obj.word);
#endif
	return buf.str();
}

as_exception::as_exception (const cl_number& obj, const char * typestring, const char * filename, int line)
	: runtime_exception(as_error_msg(obj, typestring, filename, line))
{}

}  // namespace cln
