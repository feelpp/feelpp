// class cl_spushstring.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/string/cl_spushstring.h"


// Implementation.

#include <cstring> // declares memcpy()

namespace cln {

void cl_spushstring::push (char c)
{
	if (index >= alloc) {
		var uintL newalloc = 2*alloc;
		var char* newbuffer = (char *) malloc_hook(newalloc);
		memcpy(newbuffer,buffer,alloc);
		free_hook(buffer);
		buffer = newbuffer;
		alloc = newalloc;
	}
	// Now index < alloc.
	buffer[index++] = c;
}

}  // namespace cln
