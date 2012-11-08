// cl_make_heap_string().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/string.h"


// Implementation.

#include "cln/malloc.h"
#include "base/cl_offsetof.h"

namespace cln {

cl_heap_string* cl_make_heap_string (const char * s)
{
	var unsigned long len = ::strlen(s);
	var cl_heap_string* str = (cl_heap_string*) malloc_hook(offsetofa(cl_heap_string,data)+sizeof(char)*(len+1));
	str->refcount = 1;
	str->type = &cl_class_string;
	str->length = len;
	{
		var const char* ptr1 = s;
		var char* ptr2 = &str->data[0];
		var uintL count;
		for (count = len; count > 0; count--)
			*ptr2++ = *ptr1++;
		*ptr2++ = '\0';
	}
	return str;
}

}  // namespace cln
