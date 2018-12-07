// Built-in integer ring.

#ifndef _CL_INTEGER_RING_H
#define _CL_INTEGER_RING_H

#include "cln/ring.h"
#include "cln/integer_class.h"

namespace cln {

typedef cl_specialized_number_ring<cl_I> cl_integer_ring;
extern const cl_integer_ring cl_I_ring;		// math. Z
extern cl_class cl_class_integer_ring;

class cl_I_ring_init_helper
{
	static int count;
public:
	cl_I_ring_init_helper();
	~cl_I_ring_init_helper();
};
static cl_I_ring_init_helper cl_I_ring_init_helper_instance;

}  // namespace cln

#endif /* _CL_INTEGER_RING_H */
