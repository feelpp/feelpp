// constructor cl_timing(cl_time_consumption&).

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/timing.h"


// Implementation.

namespace cln {

static void report_accu (const cl_timing& t)
{
	var const cl_time_consumption usage_end = cl_current_time_consumption();
	var const cl_time_consumption& usage_start = t.tmp;
	var cl_time_consumption usage;
	usage.realtime = usage_end.realtime - usage_start.realtime;
	usage.usertime = usage_end.usertime - usage_start.usertime;

	var cl_time_consumption& accumulator = *(cl_time_consumption*)(t.report_destination);
	accumulator.realtime = accumulator.realtime + usage.realtime;
	accumulator.usertime = accumulator.usertime + usage.usertime;
}

cl_timing::cl_timing (cl_time_consumption& accumulator)
{
	report_fn = report_accu; report_destination = &accumulator;
	tmp = cl_current_time_consumption();
}

}  // namespace cln
