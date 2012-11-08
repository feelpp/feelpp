// Null ring.

// General includes.
#include "base/cl_sysdep.h"


// Specification.
#include "cln/null_ring.h"


// Implementation.

#include "cln/integer_class.h"

namespace cln {

static const _cl_ring_element null_op0 (cl_heap_ring* R)
{
	return _cl_ring_element(R, (cl_I)0);
}

static const _cl_ring_element null_op1 (cl_heap_ring* R, const _cl_ring_element& x)
{
	unused x;
	return _cl_ring_element(R, (cl_I)0);
}

static const _cl_ring_element null_op2 (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	unused x;
	unused y;
	return _cl_ring_element(R, (cl_I)0);
}

static void null_fprint (cl_heap_ring* R, std::ostream& stream, const _cl_ring_element& x)
{
	unused R;
	unused x;
	fprint(stream,"0");
}

static bool null_equal (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	unused R;
	unused x;
	unused y;
	return true;
}

#define null_zero null_op0
static bool null_zerop (cl_heap_ring* R, const _cl_ring_element& x)
{
	unused R;
	unused x;
	return true;
}
#define null_plus null_op2
#define null_minus null_op2
#define null_uminus null_op1

#define null_one null_op0
static const _cl_ring_element null_canonhom (cl_heap_ring* R, const cl_I& x)
{
	unused x;
	return _cl_ring_element(R, (cl_I)0);
}
#define null_mul null_op2
#define null_square null_op1
static const _cl_ring_element null_expt_pos (cl_heap_ring* R, const _cl_ring_element& x, const cl_I& y)
{
	unused x;
	unused y;
	return _cl_ring_element(R, (cl_I)0);
}

static cl_ring_setops& null_setops()
{
	static cl_ring_setops ops = {
		null_fprint,
		null_equal
	};
	return ops;
}

static cl_ring_addops& null_addops()
{
	static cl_ring_addops ops = {
		null_zero,
		null_zerop,
		null_plus,
		null_minus,
		null_uminus
	};
	return ops;
}

static cl_ring_mulops& null_mulops()
{
	static cl_ring_mulops ops = {
		 null_one,
		 null_canonhom,
		 null_mul,
		 null_square,
		 null_expt_pos
	};
	return ops;
}

static const cl_class& cl_class_null_ring();

class cl_heap_null_ring : public cl_heap_ring {
	SUBCLASS_cl_heap_ring()
public:
	// Constructor.
	cl_heap_null_ring ()
		: cl_heap_ring (&null_setops(),&null_addops(),&null_mulops())
		{ type = &cl_class_null_ring(); }
	// Destructor.
	~cl_heap_null_ring () {}
};

static void cl_null_ring_destructor (cl_heap* pointer)
{
	(*(cl_heap_null_ring*)pointer).~cl_heap_null_ring();
}

static void cl_null_ring_dprint (cl_heap* pointer)
{
	unused pointer;
	fprint(cl_debugout, "(cl_null_ring) cl_0_ring");
}

static const cl_class& cl_class_null_ring()
{
	static const cl_class cl_class_null_ring_instance = {
		 cl_null_ring_destructor,
		 cl_class_flags_number_ring,
		 cl_null_ring_dprint
	};
	return cl_class_null_ring_instance;
}

cl_heap_null_ring* cl_heap_null_ring_instance;
const cl_null_ring cl_0_ring = cl_0_ring;

inline cl_null_ring::cl_null_ring ()
	: cl_ring(cl_heap_null_ring_instance) {}

int cl_0_ring_init_helper::count = 0;

cl_0_ring_init_helper::cl_0_ring_init_helper()
{
	if (count++ == 0) {
		cl_heap_null_ring_instance = new cl_heap_null_ring();
		new ((void *)&cl_0_ring) cl_null_ring();
	 }
}

cl_0_ring_init_helper::~cl_0_ring_init_helper()
{
	if (--count == 0) {
		 delete cl_heap_null_ring_instance;
	}
}

}  // namespace cln

