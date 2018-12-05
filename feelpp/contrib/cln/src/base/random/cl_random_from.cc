// random_state constructor.


#if defined(_WIN32)
#include <windows.h> // For GetCurrentProcessId(), must be included first, sorry.
#endif

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/random.h"


// Implementation.

#include "base/cl_base_config.h"
#include "base/cl_low.h"
#include <cstdlib>  // declares rand()

#if defined(unix) || defined(__unix) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__DragonFly__) || defined(_AIX) || defined(sinix) || (defined(__MACH__) && defined(__APPLE__)) || (defined(__CYGWIN__) && defined(__GNUC__)) || defined(__BEOS__)

#include <sys/types.h>
#include <unistd.h> // declares getpid()

#if defined(HAVE_GETTIMEOFDAY)

#include <sys/time.h>
#ifdef GETTIMEOFDAY_DOTS
  extern "C" int gettimeofday (struct timeval * tp, ...);
#else
  extern "C" int gettimeofday (struct timeval * tp, GETTIMEOFDAY_TZP_T tzp);
#endif

namespace cln {
inline uint32 get_seed (void)
{
	var struct timeval tv;
	gettimeofday(&tv,0);
	return highlow32(tv.tv_sec,tv.tv_usec); // 16+16 zufällige Bits
}
}  // namespace cln

#elif defined(HAVE_TIMES_CLOCK)

#include <ctime>
#ifndef CLK_TCK
#include <sys/time.h>
#endif
#include <sys/times.h>
extern "C" clock_t times (struct tms * buffer);

namespace cln {
inline uint32 get_seed (void)
{
	var struct tms tmsbuf;
	var uint32 seed_lo = times(&tmsbuf);
	return seed_lo + tmsbuf.tms_utime + tmsbuf.tms_stime;
}
}  // namespace cln

#endif

#elif defined(_WIN32)
#include <sys/time.h>
#include <sys/timeb.h>

namespace cln {
inline uint32 get_seed (void)
{
	struct timeb timebuf;
	ftime(&timebuf);
	return highlow32(timebuf.time, (long)(timebuf.millitm)*1000);
}
}  // namespace cln

#endif

namespace cln {

// Counter, to avoid that two random-states created immediately one after
// the other contain the same seed.
static uint32 counter = 0;

random_state::random_state ()
{
	var uint32 seed_hi;
	var uint32 seed_lo;
#if defined(unix) || defined(__unix) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__DragonFly__) || defined(_AIX) || defined(sinix) || (defined(__MACH__) && defined(__APPLE__)) || (defined(__CYGWIN__) && defined(__GNUC__)) || defined(__BEOS__)
	seed_lo = get_seed();
	seed_hi = (rand() // zufällige 31 Bit (bei UNIX_BSD) bzw. 16 Bit (bei UNIX_SYSV)
                          << 8) ^ (uintL)(getpid()); // ca. 8 Bit von der Process ID
#elif defined(__OpenBSD__)
	seed_lo = arc4random();
	seed_hi = arc4random();
#elif defined(_WIN32)
	seed_lo = get_seed();
	seed_hi = (rand() << 8) ^ (uintL)(GetCurrentProcessId());
#elif defined(__atarist)
	seed_lo = highlow32(GEMDOS_GetDate(),GEMDOS_GetTime()); // 16+16 zufällige Bits
	seed_hi = XBIOS_Random(); // 24 Bit zufällig vom XBIOS, vorne 8 Nullbits
#elif defined(amiga) || defined(AMIGA)
	seed_lo = get_real_time(); // Uhrzeit
	seed_hi = FindTask(NULL); // Pointer auf eigene Task
#elif defined(__MSDOS__) || defined(__EMX__) || defined(__riscos)
	// Keine Zufallszahlen, keine PID, nichts Zufälliges da.
	seed_lo = get_real_time(); // Uhrzeit, 100 Hz
	seed_hi = time(NULL);
#else
#error "Must implement random_state constructor!"
#endif
	seed_hi ^= counter++ << 5;
	seed.hi = seed_hi;
	seed.lo = seed_lo;
}

}  // namespace cln
