// Univariate Polynomial operations.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#define CL_GV_NO_RANGECHECKS
#define CL_SV_NO_RANGECHECKS
#include "cln/univpoly.h"
#include "polynomial/cl_UP.h"


// Implementation.

#include "cln/output.h"

namespace cln {

cl_symbol cl_univpoly_varname_key = (cl_symbol)(cl_string)"variable name";

// Prepare for looking into a polynomial.
  #define DeclarePoly(type,x)  \
    const type& __tmp_##x = *(const type*) &(x).rep;			\
    const type& x = __tmp_##x;
  #define DeclareMutablePoly(type,x)  \
    type& __tmp_##x = *(type*) &(x).rep;				\
    type& x = __tmp_##x;

}  // namespace cln

// Four different implementations of the polynomial operations, for efficiency:
#include "cl_UP_number.h"  // polynomials over number rings
#include "cl_UP_MI.h"      // polynomials over modular integer rings
#include "cl_UP_GF2.h"     // polynomials over the modular integer ring GF(2)
#include "cl_UP_gen.h"     // polynomials over all other rings

namespace cln {

static void cl_univpoly_ring_destructor (cl_heap* pointer)
{
	(*(cl_heap_univpoly_ring*)pointer).~cl_heap_univpoly_ring();
}

cl_class cl_class_univpoly_ring;

int cl_UP_init_helper::count = 0;

cl_UP_init_helper::cl_UP_init_helper() 
{
	if (count++ == 0) {
		cl_class_univpoly_ring.destruct = cl_univpoly_ring_destructor;
		cl_class_univpoly_ring.flags = cl_class_flags_univpoly_ring;
	}
}

cl_UP_init_helper::~cl_UP_init_helper() 
{
	if (--count == 0) {
		// nothing to clean up
	}
}

cl_heap_univpoly_ring::cl_heap_univpoly_ring (const cl_ring& r, cl_univpoly_setops* setopv, cl_univpoly_addops* addopv, cl_univpoly_mulops* mulopv, cl_univpoly_modulops* modulopv, cl_univpoly_polyops* polyopv)
	: setops (setopv), addops (addopv), mulops (mulopv), modulops (modulopv), polyops (polyopv),
	  _basering (r)
{
	refcount = 0; // will be incremented by the `cl_univpoly_ring' constructor
	type = &cl_class_univpoly_ring;
}


// Create a new univariate polynomial ring.

cl_heap_univpoly_ring* cl_make_univpoly_ring (const cl_ring& r)
{
	if (r.pointer_type()->flags & cl_class_flags_number_ring)
		return new cl_heap_num_univpoly_ring(r);
	else if (r.pointer_type()->flags & cl_class_flags_modint_ring) {
		if (((cl_heap_modint_ring*)r.heappointer)->modulus == 2)
			return new cl_heap_gf2_univpoly_ring(r);
		else
			return new cl_heap_modint_univpoly_ring(r);
	} else
		return new cl_heap_gen_univpoly_ring(r);
}

}  // namespace cln

