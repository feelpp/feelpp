// Timing tools.

#ifndef _CL_TIMING_H
#define _CL_TIMING_H

#include "cln/config.h"
#include "cln/intparam.h"
#include "cln/types.h"

#include "cln/io.h"

namespace cln {

struct cl_timespec {
	uintL tv_sec;	// seconds since 1970-01-01
	sintL tv_nsec;	// nanoseconds, >= 0, < 1000000000
	// Constructors.
	cl_timespec () {}
	cl_timespec (uintL sec, sintL nsec)
		: tv_sec (sec), tv_nsec (nsec) {}
};

struct cl_time_duration {
	uintL tv_sec;	// seconds
	uintL tv_nsec;	// nanoseconds
	// Constructors.
	cl_time_duration () {}
	cl_time_duration (uintL sec)
		: tv_sec (sec), tv_nsec (0) {}
	cl_time_duration (uintL sec, uintL nsec)
		: tv_sec (sec), tv_nsec (nsec) {}
};

struct cl_time_consumption {
	cl_time_duration realtime;	// elapsed time
	cl_time_duration usertime;	// system's notion of user time/run time
};

extern const cl_time_duration operator- (const cl_timespec&, const cl_timespec&);
extern const cl_timespec operator+ (const cl_timespec&, const cl_time_duration&);
extern const cl_timespec operator- (const cl_timespec&, const cl_time_duration&);
extern const cl_time_duration operator+ (const cl_time_duration&, const cl_time_duration&);
extern const cl_time_duration operator- (const cl_time_duration&, const cl_time_duration&);

extern const cl_timespec cl_current_time ();
extern const cl_time_consumption cl_current_time_consumption ();

// Report a time consumption.
// (Should better be a virtual member function of `cl_time_consumption').
extern void cl_timing_report (std::ostream&, const cl_time_consumption&);

struct cl_timing {
	// Constructor, starts the time interval.
	cl_timing (cl_time_consumption& accumulator);
	cl_timing (std::ostream& destination = std::cerr);
	cl_timing (const char *, std::ostream& destination = std::cerr);
	// Destructor, closes the time interval and does a report.
	~cl_timing ();	
//private:
	cl_time_consumption tmp;
	void (*report_fn) (const cl_timing&);
	void* report_destination;
	const char * comment;
};

// Macro for timing.
// Usage:
//     { CL_TIMING; computation(); }
// or  { CL_TIMING(accumulator); computation(); }
// or  { CL_TIMING(cout); computation(); }
// The timing interval starts immediately and ends at the closing brace.
#define CL_TIMING  CL_TIMING1(__LINE__)
#define CL_TIMING1(line)  CL_TIMING2(line)
#define CL_TIMING2(line)  cl_timing cl_timing_dummy_##line

}  // namespace cln

#endif /* _CL_TIMING_H */
