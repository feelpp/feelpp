// Ring of complex numbers.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex_ring.h"


// Implementation.

#include "cln/complex.h"
#include "cln/complex_io.h"
#include "complex/cl_C.h"

namespace cln {

static void N_fprint (cl_heap_ring* R, std::ostream& stream, const _cl_ring_element& x)
{
	unused R;
	fprint(stream,The(cl_N)(x));
}

static bool N_equal (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	unused R;
	return equal(The(cl_N)(x),The(cl_N)(y));
}

static const _cl_ring_element N_zero (cl_heap_ring* R)
{
	return _cl_ring_element(R, (cl_N)0);
}

static bool N_zerop (cl_heap_ring* R, const _cl_ring_element& x)
{
	unused R;
	// Here we return true only if x is the *exact* zero. Because we
	// don't want the degree of polynomials to depend on rounding errors.
	// For all ring theoretic purposes, we treat 0.0, 0+0.0i etc. as if
	// they were zero divisors.
	return exact_zerop(The(cl_N)(x));
}

static const _cl_ring_element N_plus (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	return _cl_ring_element(R, The(cl_N)(x) + The(cl_N)(y));
}

static const _cl_ring_element N_minus (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	return _cl_ring_element(R, The(cl_N)(x) - The(cl_N)(y));
}

static const _cl_ring_element N_uminus (cl_heap_ring* R, const _cl_ring_element& x)
{
	return _cl_ring_element(R, - The(cl_N)(x));
}

static const _cl_ring_element N_one (cl_heap_ring* R)
{
	return _cl_ring_element(R, (cl_N)1);
}

static const _cl_ring_element N_canonhom (cl_heap_ring* R, const cl_I& x)
{
	return _cl_ring_element(R, (cl_N)x);
}

static const _cl_ring_element N_mul (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	return _cl_ring_element(R, The(cl_N)(x) * The(cl_N)(y));
}

static const _cl_ring_element N_square (cl_heap_ring* R, const _cl_ring_element& x)
{
	return _cl_ring_element(R, square(The(cl_N)(x)));
}

static const _cl_ring_element N_expt_pos (cl_heap_ring* R, const _cl_ring_element& x, const cl_I& y)
{
	return _cl_ring_element(R, expt(The(cl_N)(x),y));
}

static bool cl_N_p (const cl_number& x)
{
	return (!x.pointer_p()
		|| (x.pointer_type()->flags & cl_class_flags_subclass_complex) != 0);
}

static cl_ring_setops N_setops = {
	N_fprint,
	N_equal
};
static cl_ring_addops N_addops = {
	N_zero,
	N_zerop,
	N_plus,
	N_minus,
	N_uminus
};
static cl_ring_mulops N_mulops = {
	N_one,
	N_canonhom,
	N_mul,
	N_square,
	N_expt_pos
};

static cl_number_ring_ops<cl_N> N_ops = {
	cl_N_p,
	equal,
	exact_zerop,
	operator+,
	operator-,
	operator-,
	operator*,
	square,
	expt
};

class cl_heap_complex_ring : public cl_heap_number_ring {
	SUBCLASS_cl_heap_ring()
public:
	// Constructor.
	cl_heap_complex_ring ()
		: cl_heap_number_ring (&N_setops,&N_addops,&N_mulops,
		                       (cl_number_ring_ops<cl_number>*) &N_ops)
		{ type = &cl_class_complex_ring; }
	// Destructor.
	~cl_heap_complex_ring () {}
};

static void cl_complex_ring_destructor (cl_heap* pointer)
{
	(*(cl_heap_complex_ring*)pointer).~cl_heap_complex_ring();
}

static void cl_complex_ring_dprint (cl_heap* pointer)
{
	unused pointer;
	fprint(cl_debugout, "(cl_complex_ring) cl_C_ring");
}

cl_class cl_class_complex_ring;
static cl_heap_complex_ring* cl_heap_complex_ring_instance;
const cl_complex_ring cl_C_ring = cl_C_ring;

// Constructor.
template <>
inline cl_complex_ring::cl_specialized_number_ring ()
	: cl_number_ring(cl_heap_complex_ring_instance) { }

int cl_C_ring_init_helper::count = 0;

cl_C_ring_init_helper::cl_C_ring_init_helper()
{
	if (count++ == 0) {
		cl_class_complex_ring.destruct = cl_complex_ring_destructor;
		cl_class_complex_ring.flags = cl_class_flags_number_ring;
		cl_class_complex_ring.dprint = cl_complex_ring_dprint;
		cl_heap_complex_ring_instance = new cl_heap_complex_ring();
		new ((void *)&cl_C_ring) cl_complex_ring();
	}
}

cl_C_ring_init_helper::~cl_C_ring_init_helper()
{
	if (--count == 0) {
		delete cl_heap_complex_ring_instance;
	}
}

}  // namespace cln

