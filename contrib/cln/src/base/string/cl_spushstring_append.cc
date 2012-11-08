// class cl_spushstring.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/string/cl_spushstring.h"


// Implementation.

#include <cstring> // declares memcpy()

namespace cln {

void cl_spushstring::append (const char * ptr, uintL len)
{
	if (index + len > alloc) {
		var uintL newalloc = index+2*len;
		if (newalloc < 2*alloc) { newalloc = 2*alloc; }
		var char* newbuffer = (char *) malloc_hook(newalloc);
		memcpy(newbuffer,buffer,alloc);
		free_hook(buffer);
		buffer = newbuffer;
		alloc = newalloc;
	}
	// Now index+len <= alloc.
	for (uintL count = len; count > 0; count--)
		buffer[index++] = *ptr++;
}

}  // namespace cln
