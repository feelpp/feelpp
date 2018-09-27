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

static void eval_pqb_series_aux (uintC N1, uintC N2,
                                 const cl_pqb_series& args,
                                 cl_I* P, cl_I* Q, cl_I* B, cl_I* T)
{
	switch (N2 - N1) {
	case 0:
		throw runtime_exception(); break;
	case 1:
		if (P) { *P = args.pv[N1]; }
		*Q = args.qv[N1];
		*B = args.bv[N1];
		*T = args.pv[N1];
		break;
	case 2: {
		var cl_I p01 = args.pv[N1] * args.pv[N1+1];
		if (P) { *P = p01; }
		*Q = args.qv[N1] * args.qv[N1+1];
		*B = args.bv[N1] * args.bv[N1+1];
		*T = args.bv[N1+1] * args.qv[N1+1] * args.pv[N1]
		   + args.bv[N1] * p01;
		break;
		}
	case 3: {
		var cl_I p01 = args.pv[N1] * args.pv[N1+1];
		var cl_I p012 = p01 * args.pv[N1+2];
		if (P) { *P = p012; }
		var cl_I q12 = args.qv[N1+1] * args.qv[N1+2];
		*Q = args.qv[N1] * q12;
		var cl_I b12 = args.bv[N1+1] * args.bv[N1+2];
		*B = args.bv[N1] * b12;
		*T = b12 * q12 * args.pv[N1]
		   + args.bv[N1] * (args.bv[N1+2] * args.qv[N1+2] * p01
		                    + args.bv[N1+1] * p012);
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
		var cl_I b01 = args.bv[N1] * args.bv[N1+1];
		var cl_I b23 = args.bv[N1+2] * args.bv[N1+3];
		*B = b01 * b23;
		*T = b23 * (args.bv[N1+1] * q123 * args.pv[N1]
		            + args.bv[N1] * q23 * p01)
		   + b01 * (args.bv[N1+3] * args.qv[N1+3] * p012
		            + args.bv[N1+2] * p0123);
		break;
		}
	default: {
		var uintC Nm = (N1+N2)/2; // midpoint
		// Compute left part.
		var cl_I LP, LQ, LB, LT;
		eval_pqb_series_aux(N1,Nm,args,&LP,&LQ,&LB,&LT);
		// Compute right part.
		var cl_I RP, RQ, RB, RT;
		eval_pqb_series_aux(Nm,N2,args,(P?&RP:(cl_I*)0),&RQ,&RB,&RT);
		// Put together partial results.
		if (P) { *P = LP*RP; }
		*Q = LQ*RQ;
		*B = LB*RB;
		// S = LS + LP/LQ * RS, so T = RB*RQ*LT + LB*LP*RT.
		*T = RB*RQ*LT + LB*LP*RT;
		break;
		}
	}
}

template<>
const cl_LF eval_rational_series<false> (uintC N, const cl_pqb_series& args, uintC len)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_I Q, B, T;
	eval_pqb_series_aux(0,N,args,NULL,&Q,&B,&T);
	return cl_I_to_LF(T,len) / cl_I_to_LF(B*Q,len);
}

static void eval_pqsb_series_aux (uintC N1, uintC N2,
                                  const cl_pqb_series& args, const uintC* qsv,
                                  cl_I* P, cl_I* Q, uintC* QS, cl_I* B, cl_I* T)
{
	switch (N2 - N1) {
	case 0:
		throw runtime_exception(); break;
	case 1:
		if (P) { *P = args.pv[N1]; }
		*Q = args.qv[N1];
		*QS = qsv[N1];
		*B = args.bv[N1];
		*T = args.pv[N1];
		break;
	case 2: {
		var cl_I p01 = args.pv[N1] * args.pv[N1+1];
		if (P) { *P = p01; }
		*Q = args.qv[N1] * args.qv[N1+1];
		*QS = qsv[N1] + qsv[N1+1];
		*B = args.bv[N1] * args.bv[N1+1];
		*T = ((args.bv[N1+1] * args.qv[N1+1] * args.pv[N1]) << qsv[N1+1])
		   + args.bv[N1] * p01;
		break;
		}
	case 3: {
		var cl_I p01 = args.pv[N1] * args.pv[N1+1];
		var cl_I p012 = p01 * args.pv[N1+2];
		if (P) { *P = p012; }
		var cl_I q12 = args.qv[N1+1] * args.qv[N1+2];
		*Q = args.qv[N1] * q12;
		*QS = qsv[N1] + qsv[N1+1] + qsv[N1+2];
		var cl_I b12 = args.bv[N1+1] * args.bv[N1+2];
		*B = args.bv[N1] * b12;
		*T = ((b12 * q12 * args.pv[N1]) << (qsv[N1+1] + qsv[N1+2]))
		   + args.bv[N1] * (((args.bv[N1+2] * args.qv[N1+2] * p01) << qsv[N1+2])
		                    + args.bv[N1+1] * p012);
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
		var cl_I b01 = args.bv[N1] * args.bv[N1+1];
		var cl_I b23 = args.bv[N1+2] * args.bv[N1+3];
		*B = b01 * b23;
		*T = ((b23 * (((args.bv[N1+1] * q123 * args.pv[N1]) << qsv[N1+1])
		              + args.bv[N1] * q23 * p01)) << (qsv[N1+2] + qsv[N1+3]))
		   + b01 * (((args.bv[N1+3] * args.qv[N1+3] * p012) << qsv[N1+3])
		            + args.bv[N1+2] * p0123);
		break;
		}
	default: {
		var uintC Nm = (N1+N2)/2; // midpoint
		// Compute left part.
		var cl_I LP, LQ, LB, LT;
		var uintC LQS;
		eval_pqsb_series_aux(N1,Nm,args,qsv,&LP,&LQ,&LQS,&LB,&LT);
		// Compute right part.
		var cl_I RP, RQ, RB, RT;
		var uintC RQS;
		eval_pqsb_series_aux(Nm,N2,args,qsv,(P?&RP:(cl_I*)0),&RQ,&RQS,&RB,&RT);
		// Put together partial results.
		if (P) { *P = LP*RP; }
		*Q = LQ*RQ;
		*QS = LQS+RQS;
		*B = LB*RB;
		// S = LS + LP/LQ * RS, so T = RB*RQ*LT + LB*LP*RT.
		*T = ((RB*RQ*LT) << RQS) + LB*LP*RT;
		break;
		}
	}
}

template<>
const cl_LF eval_rational_series<true> (uintC N, const cl_pqb_series& args, uintC len)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_I Q, B, T;
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
	eval_pqsb_series_aux(0,N,args,qsv,NULL,&Q,&QS,&B,&T);
	return cl_I_to_LF(T,len) / scale_float(cl_I_to_LF(B*Q,len),QS);
}

static void eval_pqb_series_aux (uintC N1, uintC N2,
                                 cl_pqb_series_stream& args,
                                 cl_I* P, cl_I* Q, cl_I* B, cl_I* T)
{
	switch (N2 - N1) {
	case 0:
		throw runtime_exception(); break;
	case 1: {
		var cl_pqb_series_term v0 = args.next(); // [N1]
		if (P) { *P = v0.p; }
		*Q = v0.q;
		*B = v0.b;
		*T = v0.p;
		break;
		}
	case 2: {
		var cl_pqb_series_term v0 = args.next(); // [N1]
		var cl_pqb_series_term v1 = args.next(); // [N1+1]
		var cl_I p01 = v0.p * v1.p;
		if (P) { *P = p01; }
		*Q = v0.q * v1.q;
		*B = v0.b * v1.b;
		*T = v1.b * v1.q * v0.p
		   + v0.b * p01;
		break;
		}
	case 3: {
		var cl_pqb_series_term v0 = args.next(); // [N1]
		var cl_pqb_series_term v1 = args.next(); // [N1+1]
		var cl_pqb_series_term v2 = args.next(); // [N1+2]
		var cl_I p01 = v0.p * v1.p;
		var cl_I p012 = p01 * v2.p;
		if (P) { *P = p012; }
		var cl_I q12 = v1.q * v2.q;
		*Q = v0.q * q12;
		var cl_I b12 = v1.b * v2.b;
		*B = v0.b * b12;
		*T = b12 * q12 * v0.p
		   + v0.b * (v2.b * v2.q * p01
		             + v1.b * p012);
		break;
		}
	case 4: {
		var cl_pqb_series_term v0 = args.next(); // [N1]
		var cl_pqb_series_term v1 = args.next(); // [N1+1]
		var cl_pqb_series_term v2 = args.next(); // [N1+2]
		var cl_pqb_series_term v3 = args.next(); // [N1+3]
		var cl_I p01 = v0.p * v1.p;
		var cl_I p012 = p01 * v2.p;
		var cl_I p0123 = p012 * v3.p;
		if (P) { *P = p0123; }
		var cl_I q23 = v2.q * v3.q;
		var cl_I q123 = v1.q * q23;
		*Q = v0.q * q123;
		var cl_I b01 = v0.b * v1.b;
		var cl_I b23 = v2.b * v3.b;
		*B = b01 * b23;
		*T = b23 * (v1.b * q123 * v0.p
		            + v0.b * q23 * p01)
		   + b01 * (v3.b * v3.q * p012
		            + v2.b * p0123);
		break;
		}
	default: {
		var uintC Nm = (N1+N2)/2; // midpoint
		// Compute left part.
		var cl_I LP, LQ, LB, LT;
		eval_pqb_series_aux(N1,Nm,args,&LP,&LQ,&LB,&LT);
		// Compute right part.
		var cl_I RP, RQ, RB, RT;
		eval_pqb_series_aux(Nm,N2,args,(P?&RP:(cl_I*)0),&RQ,&RB,&RT);
		// Put together partial results.
		if (P) { *P = LP*RP; }
		*Q = LQ*RQ;
		*B = LB*RB;
		// S = LS + LP/LQ * RS, so T = RB*RQ*LT + LB*LP*RT.
		*T = RB*RQ*LT + LB*LP*RT;
		break;
		}
	}
}

template<>
const cl_LF eval_rational_series<false> (uintC N, cl_pqb_series_stream& args, uintC len)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_I Q, B, T;
	eval_pqb_series_aux(0,N,args,NULL,&Q,&B,&T);
	return cl_I_to_LF(T,len) / cl_I_to_LF(B*Q,len);
}

static void eval_pqb_series_aux (uintC N1, uintC N2,
                                 cl_pqb_series_stream& args,
                                 cl_R* P, cl_R* Q, cl_R* B, cl_R* T,
                                 uintC trunclen)
{
	switch (N2 - N1) {
	case 0:
		throw runtime_exception(); break;
	case 1: {
		var cl_pqb_series_term v0 = args.next(); // [N1]
		if (P) { *P = v0.p; }
		*Q = v0.q;
		*B = v0.b;
		*T = v0.p;
		break;
		}
	case 2: {
		var cl_pqb_series_term v0 = args.next(); // [N1]
		var cl_pqb_series_term v1 = args.next(); // [N1+1]
		var cl_I p01 = v0.p * v1.p;
		if (P) { *P = p01; }
		*Q = v0.q * v1.q;
		*B = v0.b * v1.b;
		*T = v1.b * v1.q * v0.p
		   + v0.b * p01;
		break;
		}
	case 3: {
		var cl_pqb_series_term v0 = args.next(); // [N1]
		var cl_pqb_series_term v1 = args.next(); // [N1+1]
		var cl_pqb_series_term v2 = args.next(); // [N1+2]
		var cl_I p01 = v0.p * v1.p;
		var cl_I p012 = p01 * v2.p;
		if (P) { *P = p012; }
		var cl_I q12 = v1.q * v2.q;
		*Q = v0.q * q12;
		var cl_I b12 = v1.b * v2.b;
		*B = v0.b * b12;
		*T = b12 * q12 * v0.p
		   + v0.b * (v2.b * v2.q * p01
		             + v1.b * p012);
		break;
		}
	case 4: {
		var cl_pqb_series_term v0 = args.next(); // [N1]
		var cl_pqb_series_term v1 = args.next(); // [N1+1]
		var cl_pqb_series_term v2 = args.next(); // [N1+2]
		var cl_pqb_series_term v3 = args.next(); // [N1+3]
		var cl_I p01 = v0.p * v1.p;
		var cl_I p012 = p01 * v2.p;
		var cl_I p0123 = p012 * v3.p;
		if (P) { *P = p0123; }
		var cl_I q23 = v2.q * v3.q;
		var cl_I q123 = v1.q * q23;
		*Q = v0.q * q123;
		var cl_I b01 = v0.b * v1.b;
		var cl_I b23 = v2.b * v3.b;
		*B = b01 * b23;
		*T = b23 * (v1.b * q123 * v0.p
		            + v0.b * q23 * p01)
		   + b01 * (v3.b * v3.q * p012
		            + v2.b * p0123);
		break;
		}
	default: {
		var uintC Nm = (N1+N2)/2; // midpoint
		// Compute left part.
		var cl_R LP, LQ, LB, LT;
		eval_pqb_series_aux(N1,Nm,args,&LP,&LQ,&LB,&LT,trunclen);
		// Compute right part.
		var cl_R RP, RQ, RB, RT;
		eval_pqb_series_aux(Nm,N2,args,(P?&RP:(cl_I*)0),&RQ,&RB,&RT,trunclen);
		// Put together partial results.
		if (P) {
			*P = LP*RP;
			truncate_precision(*P,trunclen);
		}
		*Q = LQ*RQ;
		truncate_precision(*Q,trunclen);
		*B = LB*RB;
		truncate_precision(*B,trunclen);
		// S = LS + LP/LQ * RS, so T = RB*RQ*LT + LB*LP*RT.
		*T = RB*RQ*LT + LB*LP*RT;
		truncate_precision(*T,trunclen);
		break;
		}
	}
}

template<>
const cl_LF eval_rational_series<false> (uintC N, cl_pqb_series_stream& args, uintC len, uintC trunclen)
{
	if (N==0)
		return cl_I_to_LF(0,len);
	var cl_R Q, B, T;
	eval_pqb_series_aux(0,N,args,NULL,&Q,&B,&T,trunclen);
	return cl_R_to_LF(T,len) / cl_R_to_LF(B*Q,len);
}

// Bit complexity (if p(n), q(n), a(n), b(n) have length O(log(n))):
// O(log(N)^2*M(N)).

}  // namespace cln
