// Built-in rational number ring.

#ifndef _CL_RATIONAL_RING_H
#define _CL_RATIONAL_RING_H

#include "cln/ring.h"
#include "cln/rational_class.h"

namespace cln {

typedef cl_specialized_number_ring<cl_RA> cl_rational_ring;
extern const cl_rational_ring cl_RA_ring;	// math. Q
extern cl_class cl_class_rational_ring;

class cl_RA_ring_init_helper
{
	static int count;
public:
	cl_RA_ring_init_helper();
	~cl_RA_ring_init_helper();
};
static cl_RA_ring_init_helper cl_RA_ring_init_helper_instance;

}  // namespace cln

#endif /* _CL_RATIONAL_RING_H */
