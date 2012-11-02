// Dummy ring.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/univpoly.h"


// Implementation.

#include "cln/io.h"

namespace cln {

static const _cl_UP dummy_op0 (cl_heap_univpoly_ring* R)
{
	unused R;
	throw uninitialized_ring_exception();
}

static const _cl_UP dummy_op1 (cl_heap_univpoly_ring* R, const _cl_UP& x)
{
	unused R;
	throw uninitialized_exception(x);
}

static const _cl_UP dummy_op2 (cl_heap_univpoly_ring* R, const _cl_UP& x, const _cl_UP& y)
{
	unused R;
	throw uninitialized_exception(x, y);
}

static void dummy_fprint (cl_heap_univpoly_ring* R, std::ostream& stream, const _cl_UP& x)
{
	unused R;
	unused stream;
	throw uninitialized_exception(x);
}
static bool dummy_equal (cl_heap_univpoly_ring* R, const _cl_UP& x, const _cl_UP& y)
{
	unused R;
	throw uninitialized_exception(x, y);
}

#define dummy_zero dummy_op0
static bool dummy_zerop (cl_heap_univpoly_ring* R, const _cl_UP& x)
{
	unused R;
	throw uninitialized_exception(x);
}
#define dummy_plus dummy_op2
#define dummy_minus dummy_op2
#define dummy_uminus dummy_op1

#define dummy_one dummy_op0
static const _cl_UP dummy_canonhom (cl_heap_univpoly_ring* R, const cl_I& x)
{
	unused R;
	(void)&x; // unused x;
	throw uninitialized_ring_exception();
}
#define dummy_mul dummy_op2
#define dummy_square dummy_op1
static const _cl_UP dummy_expt_pos (cl_heap_univpoly_ring* R, const _cl_UP& x, const cl_I& y)
{
	unused R;
	(void)&y; // unused y;
	throw uninitialized_exception(x);
}

static const _cl_UP dummy_scalmul (cl_heap_univpoly_ring* R, const cl_ring_element& x, const _cl_UP& y)
{
	unused R;
	unused x;
	throw uninitialized_exception(y);
}

static sintL dummy_degree (cl_heap_univpoly_ring* R, const _cl_UP& x)
{
	unused R;
	throw uninitialized_exception(x);
}
static sintL dummy_ldegree (cl_heap_univpoly_ring* R, const _cl_UP& x)
{
	unused R;
	throw uninitialized_exception(x);
}
static const _cl_UP dummy_monomial (cl_heap_univpoly_ring* R, const cl_ring_element& x, uintL e)
{
	unused R;
	unused x;
	unused e;
	throw uninitialized_ring_exception();
}
static const cl_ring_element dummy_coeff (cl_heap_univpoly_ring* R, const _cl_UP& x, uintL index)
{
	unused R;
	unused index;
	throw uninitialized_exception(x);
}
static const _cl_UP dummy_create (cl_heap_univpoly_ring* R, sintL deg)
{
	unused R;
	unused deg;
	throw uninitialized_ring_exception();
}
static void dummy_set_coeff (cl_heap_univpoly_ring* R, _cl_UP& x, uintL index, const cl_ring_element& y)
{
	unused R;
	unused index;
	unused y;
	throw uninitialized_exception(x);
}
static void dummy_finalize (cl_heap_univpoly_ring* R, _cl_UP& x)
{
	unused R;
	throw uninitialized_exception(x);
}
static const cl_ring_element dummy_eval (cl_heap_univpoly_ring* R, const _cl_UP& x, const cl_ring_element& y)
{
	unused R;
	unused y;
	throw uninitialized_exception(x);
}

static cl_univpoly_setops dummy_setops = {
	dummy_fprint,
	dummy_equal
};
static cl_univpoly_addops dummy_addops = {
	dummy_zero,
	dummy_zerop,
	dummy_plus,
	dummy_minus,
	dummy_uminus
};
static cl_univpoly_mulops dummy_mulops = {
	dummy_one,
	dummy_canonhom,
	dummy_mul,
	dummy_square,
	dummy_expt_pos
};
static cl_univpoly_modulops dummy_modulops = {
	dummy_scalmul
};
static cl_univpoly_polyops dummy_polyops = {
	dummy_degree,
	dummy_ldegree,
	dummy_monomial,
	dummy_coeff,
	dummy_create,
	dummy_set_coeff,
	dummy_finalize,
	dummy_eval
};

class cl_heap_no_univpoly_ring : public cl_heap_univpoly_ring {
	SUBCLASS_cl_heap_univpoly_ring()
public:
	// Constructor.
	cl_heap_no_univpoly_ring ()
		: cl_heap_univpoly_ring (cl_no_ring,&dummy_setops,&dummy_addops,&dummy_mulops,&dummy_modulops,&dummy_polyops)
		{ type = &cl_class_no_univpoly_ring; }
	// Destructor.
	~cl_heap_no_univpoly_ring () {}
};

static void cl_no_univpoly_ring_destructor (cl_heap* pointer)
{
	(*(cl_heap_no_univpoly_ring*)pointer).~cl_heap_no_univpoly_ring();
}

cl_class cl_class_no_univpoly_ring;
static cl_heap_no_univpoly_ring* cl_heap_no_univpoly_ring_instance;
const cl_univpoly_ring cl_no_univpoly_ring = cl_no_univpoly_ring;

int cl_UP_no_ring_init_helper::count = 0;

cl_UP_no_ring_init_helper::cl_UP_no_ring_init_helper()
{
	if (count++ == 0) {
		cl_class_no_univpoly_ring.destruct = cl_no_univpoly_ring_destructor;
		cl_class_no_univpoly_ring.flags = 0;
		cl_heap_no_univpoly_ring_instance = new cl_heap_no_univpoly_ring();
		new ((void *)&cl_no_univpoly_ring) cl_univpoly_ring(cl_heap_no_univpoly_ring_instance);
	}
}

cl_UP_no_ring_init_helper::~cl_UP_no_ring_init_helper()
{
	if (--count == 0) {
		delete cl_heap_no_univpoly_ring_instance;
	}
}

}  // namespace cln

