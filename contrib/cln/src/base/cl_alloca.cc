// cl_alloc_alloca_header(), cl_free_alloca_header().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/cl_alloca.h"


// Implementation.

#include "cln/malloc.h"
#include "base/cl_offsetof.h"

namespace cln {

cl_alloca_header* cl_alloc_alloca_header (size_t size)
{
	var cl_alloca_header* pointer =
	  (cl_alloca_header*)malloc_hook(size+offsetofa(cl_alloca_header,usable_memory));
	pointer->next = NULL;
	return pointer;
}

void cl_free_alloca_header (cl_alloca_header* pointer)
{
	do {
		cl_alloca_header* next = pointer->next;
		free_hook(pointer);
		pointer = next;
	} while (pointer != NULL);
}

}  // namespace cln

