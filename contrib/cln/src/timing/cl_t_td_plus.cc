// operator+ (const cl_time_duration&, const cl_time_duration&)

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/timing.h"


// Implementation.

namespace cln {

const cl_time_duration operator+ (const cl_time_duration& a, const cl_time_duration& b)
{
	var uintL sum_sec = a.tv_sec + b.tv_sec;
	var uintL sum_nsec = a.tv_nsec + b.tv_nsec;
	if (sum_nsec >= 1000000000) {
		sum_nsec -= 1000000000;
		sum_sec += 1;
	}
	return cl_time_duration(sum_sec,sum_nsec);
}

}  // namespace cln
