// class cl_property_list.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/proplist.h"


// Implementation.

#include "cln/exception.h"

namespace cln {

// This tells the compiler to put the `cl_property' vtable into this file.
void cl_property::dummy () {}

void cl_property_list::add_property (cl_property* new_property)
{
	if (new_property->next)
		throw runtime_exception();
	new_property->next = list;
	list = new_property;
}

}  // namespace cln
