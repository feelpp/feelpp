// cl_current_time().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/timing.h"


// Implementation.

#include "timing/cl_t_config.h"


#if defined(HAVE_GETTIMEOFDAY)
  #include <sys/time.h>
  #ifdef GETTIMEOFDAY_DOTS
    extern "C" int gettimeofday (struct timeval * tp, ...);
  #else
    extern "C" int gettimeofday (struct timeval * tp, GETTIMEOFDAY_TZP_T tzp);
  #endif
#else
  #include <ctime>
#endif
#ifdef HAVE_PERROR_DECL
  #include <cerrno>
  #include <cstdio>
#else
  extern "C" int perror (const char *);
#endif

namespace cln {

const cl_timespec cl_current_time ()
{
#if defined(HAVE_GETTIMEOFDAY)
	var struct timeval tv;
	if (gettimeofday(&tv,NULL) != 0) {
		perror("gettimeofday");
		tv.tv_sec = 0; tv.tv_usec = 0;
	}
	return cl_timespec(tv.tv_sec,
			   tv.tv_usec * (1000000000/1000000)
			  );
#else
	return cl_timespec(time(NULL),0);
#endif
}

}  // namespace cln
