// cl_current_time().

#if defined(_WIN32) && !defined(__CYGWIN__)
#include <windows.h> // For GetSystemTimeAsFileTime(), must be included first, sorry.
#endif

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/timing.h"


// Implementation.

#include "timing/cl_t_config.h"


#if defined(HAVE_GETTIMEOFDAY)
  #include <sys/time.h>
#elif defined(_WIN32) && !defined(__CYGWIN__)
  /* <windows.h> included above. */
#else
  #include <ctime>
#endif
#include <cerrno>
#include <cstdio>

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
#elif defined(_WIN32) && !defined(__CYGWIN__)

	/* GetSystemTimePreciseAsFileTime was introduced only in Windows 8.  */
	typedef void (WINAPI * GetSystemTimePreciseAsFileTimeFuncType) (FILETIME *lpTime);
	static GetSystemTimePreciseAsFileTimeFuncType GetSystemTimePreciseAsFileTimeFunc = NULL;
	static BOOL initialized = FALSE;

	if (!initialized) {
		HMODULE kernel32 = LoadLibrary ("kernel32.dll");
		if (kernel32 != NULL) {
			GetSystemTimePreciseAsFileTimeFunc =
				(GetSystemTimePreciseAsFileTimeFuncType) (void *) GetProcAddress (kernel32, "GetSystemTimePreciseAsFileTime");
		}
		initialized = TRUE;
	}

	FILETIME current_time;

	if (GetSystemTimePreciseAsFileTimeFunc != NULL)
		GetSystemTimePreciseAsFileTimeFunc (&current_time);
	else
		GetSystemTimeAsFileTime (&current_time);

	/* Convert from FILETIME to 'struct timeval'.  */
	/* FILETIME: <https://docs.microsoft.com/en-us/windows/desktop/api/minwinbase/ns-minwinbase-filetime> */
	ULONGLONG since_1601 =
		((ULONGLONG) current_time.dwHighDateTime << 32)
		| (ULONGLONG) current_time.dwLowDateTime;
	/* Between 1601-01-01 and 1970-01-01 there were 280 normal years and 89 leap
	   years, in total 134774 days.  */
	ULONGLONG since_1970 =
		since_1601 - (ULONGLONG) 134774 * (ULONGLONG) 86400 * (ULONGLONG) 10000000;
	ULONGLONG microseconds_since_1970 = since_1970 / (ULONGLONG) 10;
	return cl_timespec(microseconds_since_1970 / (ULONGLONG) 1000000,
			   (microseconds_since_1970 % (ULONGLONG) 1000000) * (1000000000/1000000)
			  );

#else
	return cl_timespec(time(NULL),0);
#endif
}

}  // namespace cln
