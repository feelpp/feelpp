// cl_make_heap_SV_number().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/SV_number.h"

// Implementation.

namespace cln {

static void cl_svector_number_destructor (cl_heap* pointer)
{
	(*(cl_heap_SV_number*)pointer).~cl_heap_SV_number();
}

// XXX: ugh, this needs to be non-static (and non-const) so redefining
// debugging function is possible (see cl_SV_number_debug.cc) 
cl_class& cl_class_svector_number()
{
	static cl_class instance = {
		cl_svector_number_destructor,
		0
	};
	return instance;
}

cl_heap_SV_number* cl_make_heap_SV_number_uninit (std::size_t len)
{
	cl_heap_SV_number* hv = (cl_heap_SV_number*) malloc_hook(sizeof(cl_heap_SV_number)+sizeof(cl_number)*len);
	hv->refcount = 1;
	hv->type = &cl_class_svector_number();
	new (&hv->v) cl_SV_inner<cl_number> (len);
	// Have to fill hv->v[i] (0 <= i < len) yourself.
	return hv;
}

cl_heap_SV_number* cl_make_heap_SV_number (std::size_t len)
{
	cl_heap_SV_number* hv = (cl_heap_SV_number*) malloc_hook(sizeof(cl_heap_SV_number)+sizeof(cl_number)*len);
	hv->refcount = 1;
	hv->type = &cl_class_svector_number();
	new (&hv->v) cl_SV_inner<cl_number> (len);
	for (std::size_t i = 0; i < len; i++)
		init1(cl_number, hv->v[i]) (0);
	return hv;
}

// An empty vector.
const cl_SV_number cl_null_SV_number = cl_null_SV_number;

int cl_SV_number_init_helper::count = 0;

cl_SV_number_init_helper::cl_SV_number_init_helper()
{
	if (count++ == 0)
		new ((void *)&cl_null_SV_number) cl_SV_number((std::size_t)0);
}

cl_SV_number_init_helper::~cl_SV_number_init_helper()
{
	if (--count == 0) {
		// Nothing to clean up
	}
}

}  // namespace cln

