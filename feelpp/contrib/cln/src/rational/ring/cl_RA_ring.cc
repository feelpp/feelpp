// Ring of rational numbers.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational_ring.h"


// Implementation.

#include "cln/rational.h"
#include "cln/rational_io.h"
#define zerop zerop_inline
#include "rational/cl_RA.h"
#undef zerop

namespace cln {

static void RA_fprint (cl_heap_ring* R, std::ostream& stream, const _cl_ring_element& x)
{
	unused R;
	fprint(stream,The(cl_RA)(x));
}

static bool RA_equal (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	unused R;
	return equal(The(cl_RA)(x),The(cl_RA)(y));
}

static const _cl_ring_element RA_zero (cl_heap_ring* R)
{
	return _cl_ring_element(R, (cl_RA)0);
}

static bool CL_FLATTEN RA_zerop (cl_heap_ring* R, const _cl_ring_element& x)
{
	unused R;
	return zerop_inline(The(cl_RA)(x));
}

static const _cl_ring_element RA_plus (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	return _cl_ring_element(R, The(cl_RA)(x) + The(cl_RA)(y));
}

static const _cl_ring_element RA_minus (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	return _cl_ring_element(R, The(cl_RA)(x) - The(cl_RA)(y));
}

static const _cl_ring_element RA_uminus (cl_heap_ring* R, const _cl_ring_element& x)
{
	return _cl_ring_element(R, - The(cl_RA)(x));
}

static const _cl_ring_element RA_one (cl_heap_ring* R)
{
	return _cl_ring_element(R, (cl_RA)1);
}

static const _cl_ring_element RA_canonhom (cl_heap_ring* R, const cl_I& x)
{
	return _cl_ring_element(R, (cl_RA)x);
}

static const _cl_ring_element RA_mul (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	return _cl_ring_element(R, The(cl_RA)(x) * The(cl_RA)(y));
}

static const _cl_ring_element RA_square (cl_heap_ring* R, const _cl_ring_element& x)
{
	return _cl_ring_element(R, square(The(cl_RA)(x)));
}

static const _cl_ring_element RA_expt_pos (cl_heap_ring* R, const _cl_ring_element& x, const cl_I& y)
{
	return _cl_ring_element(R, expt_pos(The(cl_RA)(x),y));
}

static bool cl_RA_p (const cl_number& x)
{
	return (!x.pointer_p()
		? x.nonpointer_tag() == cl_FN_tag
		: (x.pointer_type()->flags & cl_class_flags_subclass_rational) != 0);
}

static cl_ring_setops RA_setops = {
	RA_fprint,
	RA_equal
};
static cl_ring_addops RA_addops = {
	RA_zero,
	RA_zerop,
	RA_plus,
	RA_minus,
	RA_uminus
};
static cl_ring_mulops RA_mulops = {
	RA_one,
	RA_canonhom,
	RA_mul,
	RA_square,
	RA_expt_pos
};

static cl_number_ring_ops<cl_RA> RA_ops = {
	cl_RA_p,
	equal,
	zerop,
	operator+,
	operator-,
	operator-,
	operator*,
	square,
	expt_pos
};

class cl_heap_rational_ring : public cl_heap_number_ring {
	SUBCLASS_cl_heap_ring()
public:
	// Constructor.
	cl_heap_rational_ring ()
		: cl_heap_number_ring (&RA_setops,&RA_addops,&RA_mulops,
		                       (cl_number_ring_ops<cl_number>*) &RA_ops)
		{ type = &cl_class_rational_ring; }
	// Destructor.
	~cl_heap_rational_ring () {}
};

static void cl_rational_ring_destructor (cl_heap* pointer)
{
	(*(cl_heap_rational_ring*)pointer).~cl_heap_rational_ring();
}

static void cl_rational_ring_dprint (cl_heap* pointer)
{
	unused pointer;
	fprint(cl_debugout, "(cl_rational_ring) cl_RA_ring");
}

cl_class cl_class_rational_ring;
static cl_heap_rational_ring* cl_heap_rational_ring_instance;

// Constructor.
template <>
inline cl_rational_ring::cl_specialized_number_ring ()
	: cl_number_ring(cl_heap_rational_ring_instance) { }

const cl_rational_ring cl_RA_ring = cl_RA_ring;

int cl_RA_ring_init_helper::count = 0;

cl_RA_ring_init_helper::cl_RA_ring_init_helper()
{
	if (count++ == 0) {
		cl_class_rational_ring.destruct = cl_rational_ring_destructor;
		cl_class_rational_ring.flags = cl_class_flags_number_ring;
		cl_class_rational_ring.dprint = cl_rational_ring_dprint;
		cl_heap_rational_ring_instance = new cl_heap_rational_ring();
		new ((void *)&cl_RA_ring) cl_rational_ring();
	}
}

cl_RA_ring_init_helper::~cl_RA_ring_init_helper()
{
	if (--count == 0) {
		delete cl_heap_rational_ring_instance;
	}
}

}  // namespace cln
