// fprinthexadecimal().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/io.h"


// Implementation.

namespace cln {

static void fprinthexadecimal_impl (std::ostream& stream, uintptr_t x)
{
	#define bufsize (sizeof(uintptr_t)*2)
	var char buf[bufsize+1];
	var char* bufptr = &buf[bufsize];
	*bufptr = '\0';
	do {
		uintptr_t q = x / 16;
		uintptr_t r = x % 16;
		*--bufptr = (r<10 ? '0'+r : 'A'-10+r);
		x = q;
	} while (x > 0);
	fprint(stream,bufptr);
	#undef bufsize
}

static void fprinthexadecimal_impl (std::ostream& stream, intptr_t x)
{
	if (x >= 0)
		fprinthexadecimal(stream,(uintptr_t)x);
	else {
		fprintchar(stream,'-');
		fprinthexadecimal(stream,(uintptr_t)(-1-x)+1);
	}
}

void fprinthexadecimal (std::ostream& stream, unsigned int x)
{
        fprinthexadecimal_impl(stream,(uintptr_t)x);
}
void fprinthexadecimal (std::ostream& stream, int x)
{
        fprinthexadecimal_impl(stream,(intptr_t)x);
}

void fprinthexadecimal (std::ostream& stream, unsigned long x)
{
        fprinthexadecimal_impl(stream,(uintptr_t)x);
}
void fprinthexadecimal (std::ostream& stream, long x)
{
        fprinthexadecimal_impl(stream,(intptr_t)x);
}

void fprinthexadecimal (std::ostream& stream, unsigned long long x)
{
#if long_long_bitsize <= pointer_bitsize
	fprinthexadecimal_impl(stream,(uintptr_t)x);
#else
	#define bufsize (sizeof(unsigned long long)*2)
	var char buf[bufsize+1];
	var char* bufptr = &buf[bufsize];
	*bufptr = '\0';
	do {
		unsigned long long q = x / 16;
		unsigned long long r = x % 16;
		*--bufptr = (r<10 ? '0'+r : 'A'-10+r);
		x = q;
	} while (x > 0);
	fprint(stream,bufptr);
	#undef bufsize
#endif
}

void fprinthexadecimal (std::ostream& stream, long long x)
{
#if long_long_bitsize <= pointer_bitsize
	fprinthexadecimal_impl(stream,(intptr_t)x);
#else
	if (x >= 0)
		fprinthexadecimal(stream,(unsigned long long)x);
	else {
		fprintchar(stream,'-');
		fprinthexadecimal(stream,(unsigned long long)(-1-x)+1);
	}
#endif
}

}  // namespace cln
