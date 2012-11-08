// constructor cl_timing(std::ostream&).

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/timing.h"


// Implementation.

namespace cln {

static void report_stream (const cl_timing& t)
{
	var const cl_time_consumption usage_end = cl_current_time_consumption();
	var const cl_time_consumption& usage_start = t.tmp;
	var cl_time_consumption usage;
	usage.realtime = usage_end.realtime - usage_start.realtime;
	usage.usertime = usage_end.usertime - usage_start.usertime;

	var std::ostream& destination = *(std::ostream*) t.report_destination;
	if (t.comment)
		fprint(destination,t.comment);
	cl_timing_report(destination,usage);
	fprint(destination,"\n");
}

cl_timing::cl_timing (std::ostream& destination)
{
	report_fn = report_stream;
	report_destination = &destination;
	comment = NULL;
	tmp = cl_current_time_consumption();
}

cl_timing::cl_timing (const char * msg, std::ostream& destination)
{
	report_fn = report_stream;
	report_destination = &destination;
	comment = msg;
	tmp = cl_current_time_consumption();
}

}  // namespace cln
