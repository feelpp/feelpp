// Global variables in CLN

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

namespace cln {

bool cl_inhibit_floating_point_underflow = false;

float_format_t default_float_format      = float_format_ffloat;

}  // namespace cln
