// cl_immediate_classes.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/object.h"


// Implementation.

namespace cln {

const struct cl_class * cl_immediate_classes [1<<cl_tag_len];
// Zero-initialized.

}  // namespace cln
