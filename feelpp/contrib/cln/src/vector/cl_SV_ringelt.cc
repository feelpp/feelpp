// cl_make_heap_SV_ringelt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/SV_ringelt.h"


// Implementation.

namespace cln {

static void cl_svector_ringelt_destructor (cl_heap* pointer)
{
	(*(cl_heap_SV_ringelt*)pointer).~cl_heap_SV_ringelt();
}

// XXX: this ought to be static const, but it would be impossible to
// set the printing function (see cl_SV_ringelt_debug.cc)
cl_class& cl_class_svector_ringelt()
{
	static cl_class instance = {
		cl_svector_ringelt_destructor,
		0
	};
	return instance;
}

cl_heap_SV_ringelt* cl_make_heap_SV_ringelt_uninit (std::size_t len)
{
	cl_heap_SV_ringelt* hv = (cl_heap_SV_ringelt*) malloc_hook(sizeof(cl_heap_SV_ringelt)+sizeof(_cl_ring_element)*len);
	hv->refcount = 1;
	hv->type = &cl_class_svector_ringelt();
	new (&hv->v) cl_SV_inner<_cl_ring_element> (len);
	// Have to fill hv->v[i] (0 <= i < len) yourself.
	return hv;
}

cl_heap_SV_ringelt* cl_make_heap_SV_ringelt (std::size_t len)
{
	cl_heap_SV_ringelt* hv = (cl_heap_SV_ringelt*) malloc_hook(sizeof(cl_heap_SV_ringelt)+sizeof(_cl_ring_element)*len);
	hv->refcount = 1;
	hv->type = &cl_class_svector_ringelt();
	new (&hv->v) cl_SV_inner<_cl_ring_element> (len);
	for (std::size_t i = 0; i < len; i++)
		init1(_cl_ring_element, hv->v[i]) ();
	return hv;
}

// An empty vector.
const cl_SV_ringelt cl_null_SV_ringelt = cl_null_SV_ringelt;

int cl_SV_ringelt_init_helper::count = 0;

cl_SV_ringelt_init_helper::cl_SV_ringelt_init_helper()
{
	if (count++ == 0)
		new ((void *)&cl_null_SV_ringelt) cl_SV_ringelt((std::size_t)0);
}

cl_SV_ringelt_init_helper::~cl_SV_ringelt_init_helper()
{
	if (--count == 0) {
		// Nothing to clean up here
	}
}

}  // namespace cln

