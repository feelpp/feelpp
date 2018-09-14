// cl_make_heap_string().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/string.h"


// Implementation.

#include "cln/malloc.h"
#include "base/cl_offsetof.h"

namespace cln {

cl_heap_string* cl_make_heap_string (unsigned long len)
{
	var cl_heap_string* str = (cl_heap_string*) malloc_hook(offsetofa(cl_heap_string,data)+sizeof(char)*(len+1));
	str->refcount = 1;
	str->type = &cl_class_string;
	str->length = len;
	return str;	/* Have to fill data[0..len] yourself. */
}

}  // namespace cln
