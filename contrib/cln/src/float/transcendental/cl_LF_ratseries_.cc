// eval_rational_series().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_LF_tran.h"


// Implementation.

#include "cln/lfloat.h"
#include "cln/integer.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

// Subroutine.
// Evaluates S = sum(N1 <= n < N2, a(n)/b(n) * (p(N1)...p(n))/(q(N1)...q(n)))
// and returns P = p(N1)...p(N2-1), Q = q(N1)...q(N2-1), B = B(N1)...B(N2-1)
// and T = B*Q*S (all integers). On entry N1 < N2.
// P will not be computed if a NULL pointer is passed.

static inline void eval__series_aux (uintC N1, uintC N2,
                                     const cl__series& args,
                                     cl_I* T)
{
	unused args;
	*T = N2-N1;
}

const cl_LF eval_rational_series (uintC N, const cl__series& args, uintC len)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_I T;
	eval__series_aux(0,N,args,&T);
	return cl_I_to_LF(T,len);
}
// Bit complexity (if p(n), q(n), a(n), b(n) have length O(log(n))):
// O(log(N)^2*M(N)).

}  // namespace cln
