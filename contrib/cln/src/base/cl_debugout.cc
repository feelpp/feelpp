// Debugging stream.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/io.h"


// Implementation.

// Just assume that the debugger runs on /dev/tty, independently of
// cin, cout, cerr.

#include <fstream>

namespace cln {

std::ostream * cl_debugout_stream = new std::ofstream ("/dev/tty");

}  // namespace cln
