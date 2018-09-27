// cl_current_time_consumption().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/timing.h"


// Implementation.

#include "timing/cl_t_config.h"

#if defined(HAVE_GETRUSAGE)
  #include <sys/types.h>
  #include <sys/time.h>
  #include <sys/resource.h>
  extern "C" int getrusage (RUSAGE_WHO_T who, struct rusage * rusage);
#elif defined(HAVE_SYS_TIMES_H)
  #include <sys/types.h>
  #include <sys/param.h> // defines HZ, unit for times() is 1/HZ seconds
  #include <sys/times.h>
  extern "C" clock_t times (struct tms * buffer);
#endif
#ifdef HAVE_PERROR_DECL
  #include <cerrno>
  #include <cstdio>
#else
  extern "C" int perror (const char *);
#endif

namespace cln {

const cl_time_consumption cl_current_time_consumption ()
{
	var cl_time_consumption result;

	var cl_timespec time = cl_current_time();
	result.realtime.tv_sec  = time.tv_sec;
	result.realtime.tv_nsec = time.tv_nsec;

#if defined(HAVE_GETRUSAGE)
	var struct rusage usage;
	if (getrusage(RUSAGE_SELF,&usage) == 0) {
		// use ru_utime only, ignore ru_stime.
		result.usertime.tv_sec  = usage.ru_utime.tv_sec;
		result.usertime.tv_nsec = usage.ru_utime.tv_usec * (1000000000/1000000);
	} else {
		perror("getrusage");
		result.usertime.tv_sec = 0; result.usertime.tv_nsec = 0;
	}
#elif defined(HAVE_SYS_TIMES_H)
	var struct tms usage;
	if (times(&usage) != (clock_t)(-1)) {
		// use tms_utime only, ignore tms_stime.
		var uintL used_time = usage.tms_utime;
		result.usertime.tv_sec  = used_time / HZ;
		result.usertime.tv_nsec = (used_time % HZ) * ((2*1000000000+HZ)/(2*HZ));
	} else {
		// ignore error ??
		result.usertime.tv_sec = 0; result.usertime.tv_nsec = 0;
	}
#else
	result.usertime = result.realtime;
#endif

	return result;
}

}  // namespace cln
