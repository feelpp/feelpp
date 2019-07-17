// operator- (const cl_time_duration&, const cl_time_duration&)

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/timing.h"


// Implementation.

namespace cln {

const cl_time_duration operator- (const cl_time_duration& a, const cl_time_duration& b)
{
	var sintL sec = a.tv_sec - b.tv_sec;
	var sintL nsec = a.tv_nsec - b.tv_nsec;
	if (nsec < 0) {
		nsec += 1000000000;
		sec -= 1;
	}
	if (sec < 0) {
		sec = 0; nsec = 0;
	}
	return cl_time_duration(sec,nsec);
}

}  // namespace cln
