// cl_timing_report().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/timing.h"


// Implementation.

namespace cln {

// Round to 3 decimal places.
#define CL_HZ 1000
#define CL_HZ_NSECS (1000000000/CL_HZ)

void cl_timing_report (std::ostream& stream, const cl_time_consumption& t)
{
	var uintL real_sec = t.realtime.tv_sec;
	var uintL real_msec = (t.realtime.tv_nsec + (CL_HZ_NSECS-1)/2) / CL_HZ_NSECS;
	if (real_msec >= CL_HZ) { real_msec -= CL_HZ; real_sec += 1; }
	var uintL user_sec = t.usertime.tv_sec;
	var uintL user_msec = (t.usertime.tv_nsec + (CL_HZ_NSECS-1)/2) / CL_HZ_NSECS;
	if (user_msec >= CL_HZ) { user_msec -= CL_HZ; user_sec += 1; }
	var char oldfill = stream.fill();
	var int oldwidth = stream.width();
	stream << "real time: ";
	stream.width(4); stream << real_sec; stream << ".";
	stream.fill('0'); stream.width(3); stream << real_msec;
	stream.fill(oldfill);
	stream << " s, ";
	stream << "run time: ";
	stream.width(4); stream << user_sec; stream << ".";
	stream.fill('0'); stream.width(3); stream << user_msec;
	stream.fill(oldfill);
	stream << " s";
	stream.width(oldwidth);
}

}  // namespace cln
