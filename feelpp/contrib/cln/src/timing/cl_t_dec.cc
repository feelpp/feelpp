// operator- (const cl_timespec&, const cl_time_duration&)

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/timing.h"


// Implementation.

namespace cln {

const cl_timespec operator- (const cl_timespec& a, const cl_time_duration& b)
{
	var uintL sec = a.tv_sec - b.tv_sec;
	var sintL nsec = a.tv_nsec - b.tv_nsec;
	if (nsec < 0) {
		nsec += 1000000000;
		sec -= 1;
	}
	return cl_timespec(sec,nsec);
}

}  // namespace cln
