// cl_sstring().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/string/cl_sstring.h"


// Implementation.

#include "cln/malloc.h"

namespace cln {

char * cl_sstring (const char * ptr, uintC len)
{
	var char * string = (char *) malloc_hook(len+1);
	{
		var const char* ptr1 = ptr;
		var char* ptr2 = string;
		var uintC count;
		for (count = len; count > 0; count--)
			*ptr2++ = *ptr1++;
		*ptr2++ = '\0';
	}
	return string;
}

}  // namespace cln
