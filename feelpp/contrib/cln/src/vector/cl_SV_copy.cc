// copy().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#define CL_SV_NO_RANGECHECKS
#include "cln/SV.h"


// Implementation.

#include "cln/malloc.h"

namespace cln {

const cl_SV_any copy (const cl_SV_any& src)
{
	std::size_t len = src.size();
	cl_heap_SV_any* hv = (cl_heap_SV_any*) malloc_hook(sizeof(cl_heap_SV_any)+sizeof(cl_gcobject)*len);
	hv->refcount = 1;
	hv->type = src.pointer_type();
	new (&hv->v) cl_SV_inner<cl_gcobject> (len);
	for (std::size_t i = 0; i < len; i++)
		init1(cl_gcobject, hv->v[i]) (src[i]);
	return hv;
}

}  // namespace cln
