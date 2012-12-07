// Low level: multiplication.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/cl_low.h"


// Implementation.

#ifdef NEED_VAR_mulu32_high
uint32 mulu32_high;
#endif

#ifdef NEED_FUNCTION_mulu32_
uint32 mulu32_high;
namespace cln {
uint32 mulu32_ (uint32 x, uint32 y)
{
	var uint16 x1 = high16(x);
	var uint16 x0 = low16(x);
	var uint16 y1 = high16(y);
	var uint16 y0 = low16(y);
	var uint32 hi = mulu16(x1,y1); // obere Portion
	var uint32 lo = mulu16(x0,y0); // untere Portion
	{var uint32 mid = mulu16(x0,y1); // 1. mittlere Portion
	 hi += high16(mid); mid = highlow32_0(low16(mid));
	 lo += mid; if (lo < mid) { hi += 1; } // 64-Bit-Addition
	}
	{var uint32 mid = mulu16(x1,y0); // 2. mittlere Portion
	 hi += high16(mid); mid = highlow32_0(low16(mid));
	 lo += mid; if (lo < mid) { hi += 1; } // 64-Bit-Addition
	}
	mulu32_high = hi; return lo;
}
}  // namespace cln
#endif

#ifdef NEED_FUNCTION_mulu32_w
namespace cln {
uint64 mulu32_w (uint32 arg1, uint32 arg2)
{
	var uint32 lo = mulu32_(arg1,arg2);
	var uint32 hi = mulu32_high;
	return highlow64(hi,lo);
}
}  // namespace cln
#endif


#ifdef NEED_VAR_mulu64_high
uint64 mulu64_high;
#endif

#ifdef NEED_FUNCTION_mulu64_
uint64 mulu64_high;
namespace cln {
extern "C" uint64 mulu64_ (uint64 x, uint64 y);
uint64 mulu64_ (uint64 x, uint64 y)
{
	var uint32 x1 = high32(x);
	var uint32 x0 = low32(x);
	var uint32 y1 = high32(y);
	var uint32 y0 = low32(y);
	var uint64 hi = mulu32_w(x1,y1); // obere Portion
	var uint64 lo = mulu32_w(x0,y0); // untere Portion
	{var uint64 mid = mulu32_w(x0,y1); // 1. mittlere Portion
	 hi += high32(mid); mid = highlow64_0(low32(mid));
	 lo += mid; if (lo < mid) { hi += 1; } // 128-Bit-Addition
	}
	{var uint64 mid = mulu32_w(x1,y0); // 2. mittlere Portion
	 hi += high32(mid); mid = highlow64_0(low32(mid));
	 lo += mid; if (lo < mid) { hi += 1; } // 128-Bit-Addition
	}
	mulu64_high = hi; return lo;
}
}  // namespace cln
#endif

