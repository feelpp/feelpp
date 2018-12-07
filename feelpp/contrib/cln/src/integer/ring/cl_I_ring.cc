// Ring of integers.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer_ring.h"


// Implementation.

#include "cln/integer.h"
#include "cln/integer_io.h"
#define zerop zerop_inline
#include "integer/cl_I.h"
#undef zerop

namespace cln {

static void I_fprint (cl_heap_ring* R, std::ostream& stream, const _cl_ring_element& x)
{
	unused R;
	fprint(stream,The(cl_I)(x));
}

static bool I_equal (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	unused R;
	return equal(The(cl_I)(x),The(cl_I)(y));
}

static const _cl_ring_element I_zero (cl_heap_ring* R)
{
	return _cl_ring_element(R, (cl_I)0);
}

static bool CL_FLATTEN I_zerop (cl_heap_ring* R, const _cl_ring_element& x)
{
	unused R;
	return zerop_inline(The(cl_I)(x));
}

static const _cl_ring_element I_plus (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	return _cl_ring_element(R, The(cl_I)(x) + The(cl_I)(y));
}

static const _cl_ring_element I_minus (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	return _cl_ring_element(R, The(cl_I)(x) - The(cl_I)(y));
}

static const _cl_ring_element I_uminus (cl_heap_ring* R, const _cl_ring_element& x)
{
	return _cl_ring_element(R, - The(cl_I)(x));
}

static const _cl_ring_element I_one (cl_heap_ring* R)
{
	return _cl_ring_element(R, (cl_I)1);
}

static const _cl_ring_element I_canonhom (cl_heap_ring* R, const cl_I& x)
{
	return _cl_ring_element(R, (cl_I)x);
}

static const _cl_ring_element I_mul (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	return _cl_ring_element(R, The(cl_I)(x) * The(cl_I)(y));
}

static const _cl_ring_element I_square (cl_heap_ring* R, const _cl_ring_element& x)
{
	return _cl_ring_element(R, square(The(cl_I)(x)));
}

static const _cl_ring_element I_expt_pos (cl_heap_ring* R, const _cl_ring_element& x, const cl_I& y)
{
	return _cl_ring_element(R, expt_pos(The(cl_I)(x),y));
}

static bool cl_I_p (const cl_number& x)
{
	return (!x.pointer_p()
		? x.nonpointer_tag() == cl_FN_tag
		: x.pointer_type() == &cl_class_bignum);
}

static cl_ring_setops I_setops = {
	I_fprint,
	I_equal
};
static cl_ring_addops I_addops = {
	I_zero,
	I_zerop,
	I_plus,
	I_minus,
	I_uminus
};
static cl_ring_mulops I_mulops = {
	I_one,
	I_canonhom,
	I_mul,
	I_square,
	I_expt_pos
};

static cl_number_ring_ops<cl_I> I_ops = {
	cl_I_p,
	equal,
	zerop,
	operator+,
	operator-,
	operator-,
	operator*,
	square,
	expt_pos
};

class cl_heap_integer_ring : public cl_heap_number_ring {
	SUBCLASS_cl_heap_ring()
public:
	// Constructor.
	cl_heap_integer_ring ()
		: cl_heap_number_ring (&I_setops,&I_addops,&I_mulops,
		                       (cl_number_ring_ops<cl_number>*) &I_ops)
		{ type = &cl_class_integer_ring; }
	// Destructor.
	~cl_heap_integer_ring () {}
};

static void cl_integer_ring_destructor (cl_heap* pointer)
{
	(*(cl_heap_integer_ring*)pointer).~cl_heap_integer_ring();
}

static void cl_integer_ring_dprint (cl_heap* pointer)
{
	unused pointer;
	fprint(cl_debugout, "(cl_integer_ring) cl_I_ring");
}

cl_class cl_class_integer_ring;
static cl_heap_integer_ring* cl_heap_integer_ring_instance;

// Constructor.
template <>
inline cl_integer_ring::cl_specialized_number_ring ()
	: cl_number_ring(cl_heap_integer_ring_instance) {}

const cl_integer_ring cl_I_ring = cl_I_ring;

int cl_I_ring_init_helper::count = 0;

cl_I_ring_init_helper::cl_I_ring_init_helper()
{
	if (count++ == 0) {
		cl_class_integer_ring.destruct = cl_integer_ring_destructor;
		cl_class_integer_ring.flags = cl_class_flags_number_ring;
		cl_class_integer_ring.dprint = cl_integer_ring_dprint;
		cl_heap_integer_ring_instance = new cl_heap_integer_ring();
		new ((void *)&cl_I_ring) cl_integer_ring();
	}
}

cl_I_ring_init_helper::~cl_I_ring_init_helper()
{
	if (--count == 0) {
		delete cl_heap_integer_ring_instance;
	}
}

}  // namespace cln

