// Conditions.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/condition.h"

// Implementation.

namespace cln {

// This tells the compiler to put the `cl_condition' vtable into this file.
void cl_condition::dummy () {}

// The destructor must be defined although it is virtual and abstract.
cl_condition::~cl_condition () {}

}  // namespace cln
