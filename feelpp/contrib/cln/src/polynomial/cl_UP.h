// Internals of Univariate polynomials.

#ifndef _CL_UP_H
#define _CL_UP_H

#include "cln/univpoly.h"
#include "cln/output.h"

namespace cln {

extern cl_heap_univpoly_ring* cl_make_univpoly_ring (const cl_ring& r);

struct cl_varname_property : public cl_property {
	SUBCLASS_cl_property();
public:
	cl_symbol varname;
	// Constructor.
	cl_varname_property (const cl_symbol& k, const cl_symbol& v) : cl_property (k), varname (v) {}
};

// The property list key used to look up the varname.
extern cl_symbol cl_univpoly_varname_key;

static inline const cl_string get_varname (cl_heap_univpoly_ring* UPR)
{
	cl_property* p = UPR->get_property(cl_univpoly_varname_key);
	if (p)
		return ((cl_varname_property*)p)->varname;
	else
		return default_print_flags.univpoly_varname;
}

}  // namespace cln

#endif /* _CL_UP_H */
