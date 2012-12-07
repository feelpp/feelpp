// malloc_hook, free_hook.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/malloc.h"


// Implementation.

#include <cstdlib>
#include "cln/io.h"
#include "cln/exception.h"

#ifndef malloc
  extern "C" void* malloc (size_t size);
#endif
#ifndef free
  extern "C" void free (void* ptr);
#endif

namespace cln {

// Just like malloc() but never return NULL pointers.
static void* xmalloc (size_t size)
{
	void* ptr = malloc(size);
	if (ptr)
		return ptr;
	throw runtime_exception("Out of virtual memory.");
}

void* (*malloc_hook) (size_t size) = xmalloc;
void (*free_hook) (void* ptr)      = free;

}  // namespace cln
