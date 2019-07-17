// Built-in complex number ring.

#ifndef _CL_COMPLEX_RING_H
#define _CL_COMPLEX_RING_H

#include "cln/ring.h"
#include "cln/complex_class.h"

namespace cln {

typedef cl_specialized_number_ring<cl_N> cl_complex_ring;
extern const cl_complex_ring cl_C_ring;		// math. C
extern cl_class cl_class_complex_ring;

class cl_C_ring_init_helper
{
	static int count;
public:
	cl_C_ring_init_helper();
	~cl_C_ring_init_helper();
};
static cl_C_ring_init_helper cl_C_ring_init_helper_instance;

}  // namespace cln

#endif /* _CL_COMPLEX_RING_H */
