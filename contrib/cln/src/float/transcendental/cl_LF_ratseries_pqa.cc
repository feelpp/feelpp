// eval_rational_series<bool>().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_LF_tran.h"


// Implementation.

#include "cln/lfloat.h"
#include "cln/integer.h"
#include "cln/real.h"
#include "cln/exception.h"
#include "float/lfloat/cl_LF.h"
#include "base/cl_alloca.h"

namespace cln {

// Subroutine.
// Evaluates S = sum(N1 <= n < N2, a(n)/b(n) * (p(N1)...p(n))/(q(N1)...q(n)))
// and returns P = p(N1)...p(N2-1), Q = q(N1)...q(N2-1), B = B(N1)...B(N2-1)
// and T = B*Q*S (all integers). On entry N1 < N2.
// P will not be computed if a NULL pointer is passed.

static void eval_pqa_series_aux (uintC N1, uintC N2,
                                 const cl_pqa_series& args,
                                 cl_I* P, cl_I* Q, cl_I* T)
{
	switch (N2 - N1) {
	case 0:
		throw runtime_exception(); break;
	case 1:
		if (P) { *P = args.pv[N1]; }
		*Q = args.qv[N1];
		*T = args.av[N1] * args.pv[N1];
		break;
	case 2: {
		var cl_I p01 = args.pv[N1] * args.pv[N1+1];
		if (P) { *P = p01; }
		*Q = args.qv[N1] * args.qv[N1+1];
		*T = args.qv[N1+1] * args.av[N1] * args.pv[N1]
		   + args.av[N1+1] * p01;
		break;
		}
	case 3: {
		var cl_I p01 = args.pv[N1] * args.pv[N1+1];
		var cl_I p012 = p01 * args.pv[N1+2];
		if (P) { *P = p012; }
		var cl_I q12 = args.qv[N1+1] * args.qv[N1+2];
		*Q = args.qv[N1] * q12;
		*T = q12 * args.av[N1] * args.pv[N1]
		   + args.qv[N1+2] * args.av[N1+1] * p01
		   + args.av[N1+2] * p012;
		break;
		}
	case 4: {
		var cl_I p01 = args.pv[N1] * args.pv[N1+1];
		var cl_I p012 = p01 * args.pv[N1+2];
		var cl_I p0123 = p012 * args.pv[N1+3];
		if (P) { *P = p0123; }
		var cl_I q23 = args.qv[N1+2] * args.qv[N1+3];
		var cl_I q123 = args.qv[N1+1] * q23;
		*Q = args.qv[N1] * q123;
		*T = q123 * args.av[N1] * args.pv[N1]
		   + q23 * args.av[N1+1] * p01
		   + args.qv[N1+3] * args.av[N1+2] * p012
		   + args.av[N1+3] * p0123;
		break;
		}
	default: {
		var uintC Nm = (N1+N2)/2; // midpoint
		// Compute left part.
		var cl_I LP, LQ, LT;
		eval_pqa_series_aux(N1,Nm,args,&LP,&LQ,&LT);
		// Compute right part.
		var cl_I RP, RQ, RT;
		eval_pqa_series_aux(Nm,N2,args,(P?&RP:(cl_I*)0),&RQ,&RT);
		// Put together partial results.
		if (P) { *P = LP*RP; }
		*Q = LQ*RQ;
		// S = LS + LP/LQ * RS, so T = RQ*LT + LP*RT.
		*T = RQ*LT + LP*RT;
		break;
		}
	}
}

template<>
const cl_LF eval_rational_series<false> (uintC N, const cl_pqa_series& args, uintC len)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_I Q, T;
	eval_pqa_series_aux(0,N,args,NULL,&Q,&T);
	return cl_I_to_LF(T,len) / cl_I_to_LF(Q,len);
}

static void eval_pqsa_series_aux (uintC N1, uintC N2,
                                  const cl_pqa_series& args, const uintC* qsv,
                                  cl_I* P, cl_I* Q, uintC* QS, cl_I* T)
{
	switch (N2 - N1) {
	case 0:
		throw runtime_exception(); break;
	case 1:
		if (P) { *P = args.pv[N1]; }
		*Q = args.qv[N1];
		*QS = qsv[N1];
		*T = args.av[N1] * args.pv[N1];
		break;
	case 2: {
		var cl_I p01 = args.pv[N1] * args.pv[N1+1];
		if (P) { *P = p01; }
		*Q = args.qv[N1] * args.qv[N1+1];
		*QS = qsv[N1] + qsv[N1+1];
		*T = ((args.qv[N1+1] * args.av[N1] * args.pv[N1]) << qsv[N1+1])
		   + args.av[N1+1] * p01;
		break;
		}
	case 3: {
		var cl_I p01 = args.pv[N1] * args.pv[N1+1];
		var cl_I p012 = p01 * args.pv[N1+2];
		if (P) { *P = p012; }
		var cl_I q12 = args.qv[N1+1] * args.qv[N1+2];
		*Q = args.qv[N1] * q12;
		*QS = qsv[N1] + qsv[N1+1] + qsv[N1+2];
		*T = ((q12 * args.av[N1] * args.pv[N1]) << (qsv[N1+1] + qsv[N1+2]))
		   + ((args.qv[N1+2] * args.av[N1+1] * p01) << qsv[N1+2])
		   + args.av[N1+2] * p012;
		break;
		}
	case 4: {
		var cl_I p01 = args.pv[N1] * args.pv[N1+1];
		var cl_I p012 = p01 * args.pv[N1+2];
		var cl_I p0123 = p012 * args.pv[N1+3];
		if (P) { *P = p0123; }
		var cl_I q23 = args.qv[N1+2] * args.qv[N1+3];
		var cl_I q123 = args.qv[N1+1] * q23;
		*Q = args.qv[N1] * q123;
		*QS = qsv[N1] + qsv[N1+1] + qsv[N1+2] + qsv[N1+3];
		*T = ((((((q123 * args.av[N1] * args.pv[N1]) << qsv[N1+1])
		         + q23 * args.av[N1+1] * p01) << qsv[N1+2])
		       + args.qv[N1+3] * args.av[N1+2] * p012) << qsv[N1+3])
		   + args.av[N1+3] * p0123;
		break;
		}
	default: {
		var uintC Nm = (N1+N2)/2; // midpoint
		// Compute left part.
		var cl_I LP, LQ, LT;
		var uintC LQS;
		eval_pqsa_series_aux(N1,Nm,args,qsv,&LP,&LQ,&LQS,&LT);
		// Compute right part.
		var cl_I RP, RQ, RT;
		var uintC RQS;
		eval_pqsa_series_aux(Nm,N2,args,qsv,(P?&RP:(cl_I*)0),&RQ,&RQS,&RT);
		// Put together partial results.
		if (P) { *P = LP*RP; }
		*Q = LQ*RQ;
		*QS = LQS+RQS;
		// S = LS + LP/LQ * RS, so T = RQ*LT + LP*RT.
		*T = ((RQ*LT) << RQS) + LP*RT;
		break;
		}
	}
}

template<>
const cl_LF eval_rational_series<true> (uintC N, const cl_pqa_series& args, uintC len)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_I Q, T;
	// Precomputation of the shift counts:
	// Split qv[n] into qv[n]*2^qsv[n].
	CL_ALLOCA_STACK;
	var uintC* qsv = (uintC*) cl_alloca(N*sizeof(uintC));
	var cl_I* qp = args.qv;
	var uintC* qsp = qsv;
	for (var uintC n = 0; n < N; n++, qp++, qsp++) {
		*qsp = pullout_shiftcount(*qp);
	}
	// Main computation.
	var uintC QS;
	eval_pqsa_series_aux(0,N,args,qsv,NULL,&Q,&QS,&T);
	return cl_I_to_LF(T,len) / scale_float(cl_I_to_LF(Q,len),QS);
}

static void eval_pqa_series_aux (uintC N1, uintC N2,
                                 cl_pqa_series_stream& args,
                                 cl_I* P, cl_I* Q, cl_I* T)
{
	switch (N2 - N1) {
	case 0:
		throw runtime_exception(); break;
	case 1: {
		var cl_pqa_series_term v0 = args.next(); // [N1]
		if (P) { *P = v0.p; }
		*Q = v0.q;
		*T = v0.a * v0.p;
		break;
		}
	case 2: {
		var cl_pqa_series_term v0 = args.next(); // [N1]
		var cl_pqa_series_term v1 = args.next(); // [N1+1]
		var cl_I p01 = v0.p * v1.p;
		if (P) { *P = p01; }
		*Q = v0.q * v1.q;
		*T = v1.q * v0.a * v0.p
		   + v1.a * p01;
		break;
		}
	case 3: {
		var cl_pqa_series_term v0 = args.next(); // [N1]
		var cl_pqa_series_term v1 = args.next(); // [N1+1]
		var cl_pqa_series_term v2 = args.next(); // [N1+2]
		var cl_I p01 = v0.p * v1.p;
		var cl_I p012 = p01 * v2.p;
		if (P) { *P = p012; }
		var cl_I q12 = v1.q * v2.q;
		*Q = v0.q * q12;
		*T = q12 * v0.a * v0.p
		   + v2.q * v1.a * p01
		   + v2.a * p012;
		break;
		}
	case 4: {
		var cl_pqa_series_term v0 = args.next(); // [N1]
		var cl_pqa_series_term v1 = args.next(); // [N1+1]
		var cl_pqa_series_term v2 = args.next(); // [N1+2]
		var cl_pqa_series_term v3 = args.next(); // [N1+3]
		var cl_I p01 = v0.p * v1.p;
		var cl_I p012 = p01 * v2.p;
		var cl_I p0123 = p012 * v3.p;
		if (P) { *P = p0123; }
		var cl_I q23 = v2.q * v3.q;
		var cl_I q123 = v1.q * q23;
		*Q = v0.q * q123;
		*T = q123 * v0.a * v0.p
		   + q23 * v1.a * p01
		   + v3.q * v2.a * p012
		   + v3.a * p0123;
		break;
		}
	default: {
		var uintC Nm = (N1+N2)/2; // midpoint
		// Compute left part.
		var cl_I LP, LQ, LT;
		eval_pqa_series_aux(N1,Nm,args,&LP,&LQ,&LT);
		// Compute right part.
		var cl_I RP, RQ, RT;
		eval_pqa_series_aux(Nm,N2,args,(P?&RP:(cl_I*)0),&RQ,&RT);
		// Put together partial results.
		if (P) { *P = LP*RP; }
		*Q = LQ*RQ;
		// S = LS + LP/LQ * RS, so T = RQ*LT + LP*RT.
		*T = RQ*LT + LP*RT;
		break;
		}
	}
}

template<>
const cl_LF eval_rational_series<false> (uintC N, cl_pqa_series_stream& args, uintC len)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_I Q, T;
	eval_pqa_series_aux(0,N,args,NULL,&Q,&T);
	return cl_I_to_LF(T,len) / cl_I_to_LF(Q,len);
}

static void eval_pqa_series_aux (uintC N1, uintC N2,
                                 cl_pqa_series_stream& args,
                                 cl_R* P, cl_R* Q, cl_R* T,
                                 uintC trunclen)
{
	switch (N2 - N1) {
	case 0:
		throw runtime_exception(); break;
	case 1: {
		var cl_pqa_series_term v0 = args.next(); // [N1]
		if (P) { *P = v0.p; }
		*Q = v0.q;
		*T = v0.a * v0.p;
		break;
		}
	case 2: {
		var cl_pqa_series_term v0 = args.next(); // [N1]
		var cl_pqa_series_term v1 = args.next(); // [N1+1]
		var cl_I p01 = v0.p * v1.p;
		if (P) { *P = p01; }
		*Q = v0.q * v1.q;
		*T = v1.q * v0.a * v0.p
		   + v1.a * p01;
		break;
		}
	case 3: {
		var cl_pqa_series_term v0 = args.next(); // [N1]
		var cl_pqa_series_term v1 = args.next(); // [N1+1]
		var cl_pqa_series_term v2 = args.next(); // [N1+2]
		var cl_I p01 = v0.p * v1.p;
		var cl_I p012 = p01 * v2.p;
		if (P) { *P = p012; }
		var cl_I q12 = v1.q * v2.q;
		*Q = v0.q * q12;
		*T = q12 * v0.a * v0.p
		   + v2.q * v1.a * p01
		   + v2.a * p012;
		break;
		}
	case 4: {
		var cl_pqa_series_term v0 = args.next(); // [N1]
		var cl_pqa_series_term v1 = args.next(); // [N1+1]
		var cl_pqa_series_term v2 = args.next(); // [N1+2]
		var cl_pqa_series_term v3 = args.next(); // [N1+3]
		var cl_I p01 = v0.p * v1.p;
		var cl_I p012 = p01 * v2.p;
		var cl_I p0123 = p012 * v3.p;
		if (P) { *P = p0123; }
		var cl_I q23 = v2.q * v3.q;
		var cl_I q123 = v1.q * q23;
		*Q = v0.q * q123;
		*T = q123 * v0.a * v0.p
		   + q23 * v1.a * p01
		   + v3.q * v2.a * p012
		   + v3.a * p0123;
		break;
		}
	default: {
		var uintC Nm = (N1+N2)/2; // midpoint
		// Compute left part.
		var cl_R LP, LQ, LT;
		eval_pqa_series_aux(N1,Nm,args,&LP,&LQ,&LT,trunclen);
		// Compute right part.
		var cl_R RP, RQ, RT;
		eval_pqa_series_aux(Nm,N2,args,(P?&RP:(cl_R*)0),&RQ,&RT,trunclen);
		// Put together partial results.
		if (P) {
			*P = LP*RP;
			truncate_precision(*P,trunclen);
		}
		*Q = LQ*RQ;
		truncate_precision(*Q,trunclen);
		// S = LS + LP/LQ * RS, so T = RQ*LT + LP*RT.
		*T = RQ*LT + LP*RT;
		truncate_precision(*T,trunclen);
		break;
		}
	}
}

template<>
const cl_LF eval_rational_series<false> (uintC N, cl_pqa_series_stream& args, uintC len, uintC trunclen)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_R Q, T;
	eval_pqa_series_aux(0,N,args,NULL,&Q,&T,trunclen);
	return cl_R_to_LF(T,len) / cl_R_to_LF(Q,len);
}
// Bit complexity (if p(n), q(n), a(n), b(n) have length O(log(n))):
// O(log(N)^2*M(N)).

}  // namespace cln
