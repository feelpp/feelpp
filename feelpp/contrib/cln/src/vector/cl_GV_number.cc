// cl_make_heap_GV_number().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/GV_number.h"


// Implementation.

#include "cln/exception.h"
#include "base/cl_offsetof.h"

namespace cln {

static void cl_gvector_number_destructor (cl_heap* pointer)
{
	(*(cl_heap_GV_number*)pointer).~cl_heap_GV_number();
}

// XXX: this ought to be const, but it would be impossible to register
// the printing function (in cl_GV_number_debug.cc)
cl_class& cl_class_gvector_number()
{
	static cl_class instance = {
		cl_gvector_number_destructor,
		0
	};
	return instance;
}

static inline cl_heap_GV_number * outcast (cl_GV_inner<cl_number>* vec)
{
	return (cl_heap_GV_number *)((char *) vec - offsetof(cl_heap_GV_number,v));
}
static inline const cl_heap_GV_number * outcast (const cl_GV_inner<cl_number>* vec)
{
	return (const cl_heap_GV_number *)((const char *) vec - offsetof(cl_heap_GV_number,v));
}


// Vectors of numbers.

struct cl_heap_GV_number_general : public cl_heap_GV_number {
	cl_number data[1];
	// Standard allocation disabled.
	void* operator new (size_t size) = delete;
	// Standard deallocation disabled.
	void operator delete (void* ptr) = delete;
	// No default constructor.
	cl_heap_GV_number_general ();
};

static const cl_number general_element (const cl_GV_inner<cl_number>* vec, std::size_t index)
{
	return ((const cl_heap_GV_number_general *) outcast(vec))->data[index];
}

static void general_set_element (cl_GV_inner<cl_number>* vec, std::size_t index, const cl_number& x)
{
	((cl_heap_GV_number_general *) outcast(vec))->data[index] = x;
}

static void general_do_delete (cl_GV_inner<cl_number>* vec)
{
	cl_heap_GV_number_general* hv = (cl_heap_GV_number_general *) outcast(vec);
	std::size_t len = hv->v.size();
	for (std::size_t i = 0; i < len; i++)
		hv->data[i].~cl_number();
}

static void general_copy_elements (const cl_GV_inner<cl_number>* srcvec, std::size_t srcindex, cl_GV_inner<cl_number>* destvec, std::size_t destindex, std::size_t count)
{
	if (count > 0) {
		const cl_heap_GV_number_general* srcv =
		  (const cl_heap_GV_number_general *) outcast(srcvec);
		cl_heap_GV_number_general* destv =
		  (cl_heap_GV_number_general *) outcast(destvec);
		std::size_t srclen = srcv->v.size();
		std::size_t destlen = destv->v.size();
		if (!(srcindex <= srcindex+count && srcindex+count <= srclen))
			throw runtime_exception();
		if (!(destindex <= destindex+count && destindex+count <= destlen))
			throw runtime_exception();
		do {
			destv->data[destindex++] = srcv->data[srcindex++];
		} while (--count > 0);
	}
}


cl_heap_GV_number* cl_make_heap_GV_number (std::size_t len)
{
	static cl_GV_vectorops<cl_number> general_vectorops = {
		general_element,
		general_set_element,
		general_do_delete,
		general_copy_elements
	};

	cl_heap_GV_number_general* hv = (cl_heap_GV_number_general*) malloc_hook(offsetofa(cl_heap_GV_number_general,data)+sizeof(cl_number)*len);
	hv->refcount = 1;
	hv->type = &cl_class_gvector_number();
	new (&hv->v) cl_GV_inner<cl_number> (len,&general_vectorops);
	for (std::size_t i = 0; i < len; i++)
		init1(cl_number, hv->data[i]) ();
	return hv;
}

// An empty vector.
const cl_GV_number cl_null_GV_number = cl_null_GV_number;

int cl_GV_number_init_helper::count = 0;

cl_GV_number_init_helper::cl_GV_number_init_helper()
{
	if (count++ == 0)
		new ((void *)&cl_null_GV_number) cl_GV_number((std::size_t)0);
}

cl_GV_number_init_helper::~cl_GV_number_init_helper()
{
	if (--count == 0) {
		// Nothing to clean up
	}
}

}  // namespace cln

