// cl_string debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/string.h"
#include "cln/io.h"
#include <cctype>

namespace cln {

static void dprint (cl_heap* pointer)
{
	var const cl_string& obj = *(const cl_string*)&pointer;
	fprint(cl_debugout, "(cl_string) \"");
	var unsigned long l = obj.size();
	for (var unsigned long i = 0; i < l; i++) {
		var unsigned char c = obj[i];
		if (c >= 0x20) {
			if (c == '"' || c == '\\')
				fprintchar(cl_debugout, '\\');
			fprintchar(cl_debugout, c);
		} else
		switch (c) {
			case '\n': fprint(cl_debugout, "\\n"); break;
			case '\t': fprint(cl_debugout, "\\t"); break;
			case '\b': fprint(cl_debugout, "\\b"); break;
			case '\r': fprint(cl_debugout, "\\r"); break;
			case '\f': fprint(cl_debugout, "\\f"); break;
			case '\v': fprint(cl_debugout, "\\v"); break;
			default:
				static const char hexdigits[] = "0123456789abcdef";
				fprintchar(cl_debugout, '\\');
				fprintchar(cl_debugout, 'x');
				fprintchar(cl_debugout, hexdigits[(c>>4)&0x0f]);
				fprintchar(cl_debugout, hexdigits[c&0x0f]);
				break;
		}
	}
	fprint(cl_debugout, "\"");
}
AT_INITIALIZATION(dprint_string)
{ cl_register_type_printer(cl_class_string,dprint); }

// This dummy links in this module when <cln/string.h> requires it.
int cl_string_debug_module;

}  // namespace cln
