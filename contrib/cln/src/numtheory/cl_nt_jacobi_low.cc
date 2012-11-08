// jacobi().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/numtheory.h"


// Implementation.

#include "cln/exception.h"
#include "base/cl_xmacros.h"

namespace cln {

// Assume 0 <= a < b.
inline int jacobi_aux (uintV a, uintV b)
{
	var int v = 1;
	for (;;) {
		// (a/b) * v is invariant.
		if (b == 1)
			// b=1 implies (a/b) = 1.
			return v;
		if (a == 0)
			// b>1 and a=0 imply (a/b) = 0.
			return 0;
		if (a > (b >> 1)) {
			// a > b/2, so (a/b) = (-1/b) * ((b-a)/b),
			// and (-1/b) = -1 if b==3 mod 4.
			a = b-a;
			switch (b % 4) {
				case 1: break;
				case 3: v = -v; break;
				default: throw runtime_exception();
			}
			continue;
		}
		if ((a & 1) == 0) {
			// b>1 and a=2a', so (a/b) = (2/b) * (a'/b),
			// and (2/b) = -1 if b==3,5 mod 8.
			a = a>>1;
			switch (b % 8) {
				case 1: case 7: break;
				case 3: case 5: v = -v; break;
				default: throw runtime_exception();
			}
			continue;
		}
		// a and b odd, 0 < a < b/2 < b, so apply quadratic reciprocity
		// law  (a/b) = (-1)^((a-1)/2)((b-1)/2) * (b/a).
		if ((a & b & 3) == 3)
			v = -v;
		swap(uintV, a,b);
		// Now a > 2*b, set a := a mod b.
		if ((a >> 3) >= b)
			a = a % b;
		else
			{ a = a-b; do { a = a-b; } while (a >= b); }
	}
}

int jacobi (sintV a, sintV b)
{
	// Check b > 0, b odd.
	if (!(b > 0))
		throw runtime_exception();
	if ((b & 1) == 0)
		throw runtime_exception();
	// Ensure 0 <= a < b.
	if (a >= 0)
		a = (uintV)a % (uintV)b;
	else
		a = b-1-((uintV)(~a) % (uintV)b);
	return jacobi_aux(a,b);
}

}  // namespace cln
