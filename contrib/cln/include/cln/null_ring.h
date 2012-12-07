// Built-in null ring.

#ifndef _CL_NULL_RING_H
#define _CL_NULL_RING_H

#include "cln/ring.h"

namespace cln {

class cl_null_ring : public cl_ring { public: cl_null_ring (); };
extern const cl_null_ring cl_0_ring;		// math. {0}

class cl_0_ring_init_helper
{
	static int count;
public:
	cl_0_ring_init_helper();
	~cl_0_ring_init_helper();
};
static cl_0_ring_init_helper cl_0_ring_init_helper_instance;

}  // namespace cln

#endif /* _CL_NULL_RING_H */
