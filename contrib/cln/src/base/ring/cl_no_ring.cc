// Dummy ring.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ring.h"


// Implementation.

#include <sstream>
#include "cln/io.h"

namespace cln {

uninitialized_ring_exception::uninitialized_ring_exception ()
	: runtime_exception("Uninitialized ring operation called.")
{}

static inline const std::string
uninitialized_error_msg (const _cl_ring_element& obj)
{
	std::ostringstream buf;
	fprint(buf, "Uninitialized ring element @0x");
	fprinthexadecimal(buf, (unsigned long)(void*)&obj);
	fprint(buf, ": 0x");
        fprinthexadecimal(buf, (unsigned long)obj.rep.word);
	return buf.str();
}

static inline const std::string
uninitialized_error_msg (const _cl_ring_element& obj_x, const _cl_ring_element& obj_y)
{
	std::ostringstream buf;
	fprint(buf, "Uninitialized ring elements @0x");
	fprinthexadecimal(buf, (unsigned long)(void*)&obj_x);
	fprint(buf, ": 0x");
        fprinthexadecimal(buf, (unsigned long)obj_x.rep.word);
	fprint(buf, ", @0x");
	fprinthexadecimal(buf, (unsigned long)(void*)&obj_y);
	fprint(buf, ": 0x");
        fprinthexadecimal(buf, (unsigned long)obj_y.rep.word);
	return buf.str();
}

uninitialized_exception::uninitialized_exception (const _cl_ring_element& obj)
	: runtime_exception(uninitialized_error_msg(obj))
{}

uninitialized_exception::uninitialized_exception (const _cl_ring_element& obj_x, const _cl_ring_element& obj_y)
	: runtime_exception(uninitialized_error_msg(obj_x, obj_y))
{}


static const _cl_ring_element dummy_op0 (cl_heap_ring* R)
{
	unused R;
	throw uninitialized_ring_exception();
}

static const _cl_ring_element dummy_op1 (cl_heap_ring* R, const _cl_ring_element& x)
{
	unused R;
	throw uninitialized_exception(x);
}

static const _cl_ring_element dummy_op2 (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	unused R;
	throw uninitialized_exception(x, y);
}

static void dummy_fprint (cl_heap_ring* R, std::ostream& stream, const _cl_ring_element& x)
{
	unused R;
	unused stream;
	throw uninitialized_exception(x);
}
static bool dummy_equal (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y)
{
	unused R;
	throw uninitialized_exception(x, y);
}

#define dummy_zero dummy_op0
static bool dummy_zerop (cl_heap_ring* R, const _cl_ring_element& x)
{
	unused R;
	throw uninitialized_exception(x);
}
#define dummy_plus dummy_op2
#define dummy_minus dummy_op2
#define dummy_uminus dummy_op1

#define dummy_one dummy_op0
static const _cl_ring_element dummy_canonhom (cl_heap_ring* R, const cl_I& x)
{
	unused R;
	(void)&x; // unused x;
	throw uninitialized_ring_exception();
}
#define dummy_mul dummy_op2
#define dummy_square dummy_op1
static const _cl_ring_element dummy_expt_pos (cl_heap_ring* R, const _cl_ring_element& x, const cl_I& y)
{
	unused R;
	(void)&y; // unused y;
	throw uninitialized_exception(x);
}

static cl_ring_setops dummy_setops = {
	dummy_fprint,
	dummy_equal
};
static cl_ring_addops dummy_addops = {
	dummy_zero,
	dummy_zerop,
	dummy_plus,
	dummy_minus,
	dummy_uminus
};
static cl_ring_mulops dummy_mulops = {
	dummy_one,
	dummy_canonhom,
	dummy_mul,
	dummy_square,
	dummy_expt_pos
};

class cl_heap_no_ring : public cl_heap_ring {
	SUBCLASS_cl_heap_ring()
public:
	// Constructor.
	cl_heap_no_ring ()
		: cl_heap_ring (&dummy_setops,&dummy_addops,&dummy_mulops)
		{ type = &cl_class_no_ring; }
	// Destructor.
	~cl_heap_no_ring () {}
};

static void cl_no_ring_destructor (cl_heap* pointer)
{
	(*(cl_heap_no_ring*)pointer).~cl_heap_no_ring();
}

static void cl_no_ring_dprint (cl_heap* pointer)
{
	unused pointer;
	fprint(cl_debugout, "(cl_ring) cl_no_ring");
}

cl_class cl_class_no_ring;

static cl_heap_no_ring* cl_heap_no_ring_instance;
// const cl_ring cl_no_ring = cl_ring (new cl_heap_no_ring());
const cl_ring cl_no_ring = cl_no_ring;


int cl_no_ring_init_helper::count = 0;

cl_no_ring_init_helper::cl_no_ring_init_helper()
{
	if (count++ == 0) {
		cl_class_no_ring.destruct = cl_no_ring_destructor;
		cl_class_no_ring.flags = 0;
		cl_class_no_ring.dprint = cl_no_ring_dprint;

		cl_heap_no_ring_instance = new cl_heap_no_ring();
		new((void*)&cl_no_ring) cl_ring(cl_heap_no_ring_instance);
	}
}

cl_no_ring_init_helper::~cl_no_ring_init_helper()
{
	if (--count == 0)
		delete cl_heap_no_ring_instance;
}

}  // namespace cln

