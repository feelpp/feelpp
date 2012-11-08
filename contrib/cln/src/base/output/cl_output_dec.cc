// fprintdecimal().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/io.h"


// Implementation.

namespace cln {

// We don't use `stream << x' or `stream << dec << x', because an ostream
// carries so many attributes, and we don't want to modifies these attributes.

void fprintdecimal (std::ostream& stream, unsigned long x)
{
	#define bufsize 20
	var char buf[bufsize+1];
	var char* bufptr = &buf[bufsize];
	*bufptr = '\0';
	do {
		unsigned long q = x / 10;
		unsigned long r = x % 10;
		*--bufptr = '0'+r;
		x = q;
	} while (x > 0);
	fprint(stream,bufptr);
	#undef bufsize
}

void fprintdecimal (std::ostream& stream, long x)
{
	if (x >= 0)
		fprintdecimal(stream,(unsigned long)x);
	else {
		fprintchar(stream,'-');
		fprintdecimal(stream,(unsigned long)(-1-x)+1);
	}
}

}  // namespace cln
