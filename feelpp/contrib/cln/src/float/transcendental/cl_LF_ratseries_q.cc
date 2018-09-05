// eval_rational_series<bool>().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_LF_tran.h"


// Implementation.

#include "cln/lfloat.h"
#include "cln/integer.h"
#include "cln/exception.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

// Subroutine.
// Evaluates S = sum(N1 <= n < N2, a(n)/b(n) * (p(N1)...p(n))/(q(N1)...q(n)))
// and returns P = p(N1)...p(N2-1), Q = q(N1)...q(N2-1), B = B(N1)...B(N2-1)
// and T = B*Q*S (all integers). On entry N1 < N2.
// P will not be computed if a NULL pointer is passed.

static void eval_q_series_aux (uintC N1, uintC N2,
                               const cl_q_series& args,
                               cl_I* Q, cl_I* T)
{
	switch (N2 - N1) {
	case 0:
		throw runtime_exception(); break;
	case 1:
		*Q = args.qv[N1];
		*T = 1;
		break;
	case 2: {
		*Q = args.qv[N1] * args.qv[N1+1];
		*T = args.qv[N1+1]
		   + 1;
		break;
		}
	case 3: {
		var cl_I q12 = args.qv[N1+1] * args.qv[N1+2];
		*Q = args.qv[N1] * q12;
		*T = q12
		   + args.qv[N1+2]
		   + 1;
		break;
		}
	case 4: {
		var cl_I q23 = args.qv[N1+2] * args.qv[N1+3];
		var cl_I q123 = args.qv[N1+1] * q23;
		*Q = args.qv[N1] * q123;
		*T = q123
		   + q23
		   + args.qv[N1+3]
		   + 1;
		break;
		}
	default: {
		var uintC Nm = (N1+N2)/2; // midpoint
		// Compute left part.
		var cl_I LQ, LT;
		eval_q_series_aux(N1,Nm,args,&LQ,&LT);
		// Compute right part.
		var cl_I RQ, RT;
		eval_q_series_aux(Nm,N2,args,&RQ,&RT);
		// Put together partial results.
		*Q = LQ*RQ;
		// S = LS + 1/LQ * RS, so T = RQ*LT + RT.
		*T = RQ*LT + RT;
		break;
		}
	}
}

template<>
const cl_LF eval_rational_series<false> (uintC N, const cl_q_series& args, uintC len)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_I Q, T;
	eval_q_series_aux(0,N,args,&Q,&T);
	return cl_I_to_LF(T,len) / cl_I_to_LF(Q,len);
}

static void eval_q_series_aux (uintC N1, uintC N2,
                               cl_q_series_stream& args,
                               cl_I* Q, cl_I* T)
{
	switch (N2 - N1) {
	case 0:
		throw runtime_exception(); break;
	case 1: {
		var cl_q_series_term v0 = args.next(); // [N1]
		*Q = v0.q;
		*T = 1;
		break;
		}
	case 2: {
		var cl_q_series_term v0 = args.next(); // [N1]
		var cl_q_series_term v1 = args.next(); // [N1+1]
		*Q = v0.q * v1.q;
		*T = v1.q + 1;
		break;
		}
	case 3: {
		var cl_q_series_term v0 = args.next(); // [N1]
		var cl_q_series_term v1 = args.next(); // [N1+1]
		var cl_q_series_term v2 = args.next(); // [N1+2]
		var cl_I q12 = v1.q * v2.q;
		*Q = v0.q * q12;
		*T = q12 + v2.q + 1;
		break;
		}
	case 4: {
		var cl_q_series_term v0 = args.next(); // [N1]
		var cl_q_series_term v1 = args.next(); // [N1+1]
		var cl_q_series_term v2 = args.next(); // [N1+2]
		var cl_q_series_term v3 = args.next(); // [N1+3]
		var cl_I q23 = v2.q * v3.q;
		var cl_I q123 = v1.q * q23;
		*Q = v0.q * q123;
		*T = q123 + q23 + v3.q + 1;
		break;
		}
	default: {
		var uintC Nm = (N1+N2)/2; // midpoint
		// Compute left part.
		var cl_I LQ, LT;
		eval_q_series_aux(N1,Nm,args,&LQ,&LT);
		// Compute right part.
		var cl_I RQ, RT;
		eval_q_series_aux(Nm,N2,args,&RQ,&RT);
		// Put together partial results.
		*Q = LQ*RQ;
		// S = LS + 1/LQ * RS, so T = RQ*LT + RT.
		*T = RQ*LT + RT;
		break;
		}
	}
}

template<>
const cl_LF eval_rational_series<false> (uintC N, cl_q_series_stream& args, uintC len)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_I Q, T;
	eval_q_series_aux(0,N,args,&Q,&T);
	return cl_I_to_LF(T,len) / cl_I_to_LF(Q,len);
}
// Bit complexity (if p(n), q(n), a(n), b(n) have length O(log(n))):
// O(log(N)^2*M(N)).

}  // namespace cln
