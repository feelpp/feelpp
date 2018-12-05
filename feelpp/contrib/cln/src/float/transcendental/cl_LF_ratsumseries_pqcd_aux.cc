// eval_pqcd_series_aux().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_LF_tran.h"


// Implementation.

#include "cln/integer.h"
#include "cln/real.h"
#include "cln/exception.h"

namespace cln {

void eval_pqcd_series_aux (uintC N, cl_pqcd_series_term* args, cl_pqcd_series_result<cl_I>& Z, bool rightmost)
{
	// N = N2-N1
	switch (N) {
	case 0:
		throw runtime_exception(); break;
	case 1:
		if (!rightmost) { Z.P = args[0].p; }
		Z.Q = args[0].q;
		Z.T = args[0].p;
		if (!rightmost) { Z.C = args[0].c; }
		Z.D = args[0].d;
		Z.V = args[0].c * args[0].p;
		break;
	case 2: {
		var cl_I p01 = args[0].p * args[1].p;
		if (!rightmost) { Z.P = p01; }
		Z.Q = args[0].q * args[1].q;
		var cl_I p0q1 = args[0].p * args[1].q + p01;
		Z.T = p0q1;
		var cl_I c0d1 = args[0].c * args[1].d;
		var cl_I c1d0 = args[1].c * args[0].d;
		if (!rightmost) { Z.C = c0d1 + c1d0; }
		Z.D = args[0].d * args[1].d;
		Z.V = c0d1 * p0q1 + c1d0 * p01;
		break;
		}
	case 3: {
		var cl_I p01 = args[0].p * args[1].p;
		var cl_I p012 = p01 * args[2].p;
		if (!rightmost) { Z.P = p012; }
		Z.Q = args[0].q * args[1].q * args[2].q;
		var cl_I p0q1 = args[0].p * args[1].q + p01;
		Z.T = args[2].q * p0q1 + p012;
		var cl_I c0d1 = args[0].c * args[1].d;
		var cl_I c1d0 = args[1].c * args[0].d;
		var cl_I d01 = args[0].d * args[1].d;
		if (!rightmost) { Z.C = (c0d1 + c1d0) * args[2].d + args[2].c * d01; }
		Z.D = d01 * args[2].d;
		Z.V = args[2].d * (args[2].q * (c0d1 * p0q1 + c1d0 * p01) + (c0d1 + c1d0) * p012) + args[2].c * d01 * p012;
		break;
		}
	default: {
		var uintC Nm = N/2; // midpoint
		// Compute left part.
		var cl_pqcd_series_result<cl_I> L;
		eval_pqcd_series_aux(Nm,args+0,L,false);
		// Compute right part.
		var cl_pqcd_series_result<cl_I> R;
		eval_pqcd_series_aux(N-Nm,args+Nm,R,rightmost);
		// Put together partial results.
		if (!rightmost) { Z.P = L.P * R.P; }
		Z.Q = L.Q * R.Q;
		// Z.S = L.S + L.P/L.Q*R.S;
		var cl_I tmp = L.P * R.T;
		Z.T = R.Q * L.T + tmp;
		if (!rightmost) { Z.C = L.C * R.D + L.D * R.C; }
		Z.D = L.D * R.D;
		// Z.U = L.U + L.C/L.D * L.P/L.Q * R.S + L.P/L.Q * R.U;
		// Z.V = R.D * R.Q * L.V + R.D * L.C * L.P * R.T + L.D * L.P * R.V;
		Z.V = R.D * (R.Q * L.V + L.C * tmp) + L.D * L.P * R.V;
		break;
		}
	}
}

void eval_pqcd_series_aux (uintC N, cl_pqcd_series_stream& args, cl_pqcd_series_result<cl_I>& Z, bool rightmost)
{
	// N = N2-N1
	switch (N) {
	case 0:
		throw runtime_exception(); break;
	case 1: {
		var cl_pqcd_series_term v0 = args.next(); // [N1]
		if (!rightmost) { Z.P = v0.p; }
		Z.Q = v0.q;
		Z.T = v0.p;
		if (!rightmost) { Z.C = v0.c; }
		Z.D = v0.d;
		Z.V = v0.c * v0.p;
		break;
		}
	case 2: {
		var cl_pqcd_series_term v0 = args.next(); // [N1]
		var cl_pqcd_series_term v1 = args.next(); // [N1+1]
		var cl_I p01 = v0.p * v1.p;
		if (!rightmost) { Z.P = p01; }
		Z.Q = v0.q * v1.q;
		var cl_I p0q1 = v0.p * v1.q + p01;
		Z.T = p0q1;
		var cl_I c0d1 = v0.c * v1.d;
		var cl_I c1d0 = v1.c * v0.d;
		if (!rightmost) { Z.C = c0d1 + c1d0; }
		Z.D = v0.d * v1.d;
		Z.V = c0d1 * p0q1 + c1d0 * p01;
		break;
		}
	case 3: {
		var cl_pqcd_series_term v0 = args.next(); // [N1]
		var cl_pqcd_series_term v1 = args.next(); // [N1+1]
		var cl_pqcd_series_term v2 = args.next(); // [N1+2]
		var cl_I p01 = v0.p * v1.p;
		var cl_I p012 = p01 * v2.p;
		if (!rightmost) { Z.P = p012; }
		Z.Q = v0.q * v1.q * v2.q;
		var cl_I p0q1 = v0.p * v1.q + p01;
		Z.T = v2.q * p0q1 + p012;
		var cl_I c0d1 = v0.c * v1.d;
		var cl_I c1d0 = v1.c * v0.d;
		var cl_I d01 = v0.d * v1.d;
		if (!rightmost) { Z.C = (c0d1 + c1d0) * v2.d + v2.c * d01; }
		Z.D = d01 * v2.d;
		Z.V = v2.d * (v2.q * (c0d1 * p0q1 + c1d0 * p01) + (c0d1 + c1d0) * p012) + v2.c * d01 * p012;
		break;
		}
	default: {
		var uintC Nm = N/2; // midpoint
		// Compute left part.
		var cl_pqcd_series_result<cl_I> L;
		eval_pqcd_series_aux(Nm,args,L,false);
		// Compute right part.
		var cl_pqcd_series_result<cl_I> R;
		eval_pqcd_series_aux(N-Nm,args,R,rightmost);
		// Put together partial results.
		if (!rightmost) { Z.P = L.P * R.P; }
		Z.Q = L.Q * R.Q;
		// Z.S = L.S + L.P/L.Q*R.S;
		var cl_I tmp = L.P * R.T;
		Z.T = R.Q * L.T + tmp;
		if (!rightmost) { Z.C = L.C * R.D + L.D * R.C; }
		Z.D = L.D * R.D;
		// Z.U = L.U + L.C/L.D * L.P/L.Q * R.S + L.P/L.Q * R.U;
		// Z.V = R.D * R.Q * L.V + R.D * L.C * L.P * R.T + L.D * L.P * R.V;
		Z.V = R.D * (R.Q * L.V + L.C * tmp) + L.D * L.P * R.V;
		break;
		}
	}
}

void eval_pqcd_series_aux (uintC N, cl_pqcd_series_stream& args, cl_pqcd_series_result<cl_R>& Z, uintC trunclen, bool rightmost)
{
	// N = N2-N1
	switch (N) {
	case 0:
		throw runtime_exception(); break;
	case 1: {
		var cl_pqcd_series_term v0 = args.next(); // [N1]
		if (!rightmost) { Z.P = v0.p; }
		Z.Q = v0.q;
		Z.T = v0.p;
		if (!rightmost) { Z.C = v0.c; }
		Z.D = v0.d;
		Z.V = v0.c * v0.p;
		break;
		}
	case 2: {
		var cl_pqcd_series_term v0 = args.next(); // [N1]
		var cl_pqcd_series_term v1 = args.next(); // [N1+1]
		var cl_I p01 = v0.p * v1.p;
		if (!rightmost) { Z.P = p01; }
		Z.Q = v0.q * v1.q;
		var cl_I p0q1 = v0.p * v1.q + p01;
		Z.T = p0q1;
		var cl_I c0d1 = v0.c * v1.d;
		var cl_I c1d0 = v1.c * v0.d;
		if (!rightmost) { Z.C = c0d1 + c1d0; }
		Z.D = v0.d * v1.d;
		Z.V = c0d1 * p0q1 + c1d0 * p01;
		break;
		}
	case 3: {
		var cl_pqcd_series_term v0 = args.next(); // [N1]
		var cl_pqcd_series_term v1 = args.next(); // [N1+1]
		var cl_pqcd_series_term v2 = args.next(); // [N1+2]
		var cl_I p01 = v0.p * v1.p;
		var cl_I p012 = p01 * v2.p;
		if (!rightmost) { Z.P = p012; }
		Z.Q = v0.q * v1.q * v2.q;
		var cl_I p0q1 = v0.p * v1.q + p01;
		Z.T = v2.q * p0q1 + p012;
		var cl_I c0d1 = v0.c * v1.d;
		var cl_I c1d0 = v1.c * v0.d;
		var cl_I d01 = v0.d * v1.d;
		if (!rightmost) { Z.C = (c0d1 + c1d0) * v2.d + v2.c * d01; }
		Z.D = d01 * v2.d;
		Z.V = v2.d * (v2.q * (c0d1 * p0q1 + c1d0 * p01) + (c0d1 + c1d0) * p012) + v2.c * d01 * p012;
		break;
		}
	default: {
		var uintC Nm = N/2; // midpoint
		// Compute left part.
		var cl_pqcd_series_result<cl_R> L;
		eval_pqcd_series_aux(Nm,args,L,trunclen,false);
		// Compute right part.
		var cl_pqcd_series_result<cl_R> R;
		eval_pqcd_series_aux(N-Nm,args,R,trunclen,rightmost);
		// Put together partial results.
		if (!rightmost) {
			Z.P = L.P * R.P;
			truncate_precision(Z.P,trunclen);
		}
		Z.Q = L.Q * R.Q;
		truncate_precision(Z.Q,trunclen);
		// Z.S = L.S + L.P/L.Q*R.S;
		var cl_R tmp = L.P * R.T;
		Z.T = R.Q * L.T + tmp;
		truncate_precision(Z.T,trunclen);
		if (!rightmost) {
			Z.C = L.C * R.D + L.D * R.C;
			truncate_precision(Z.C,trunclen);
		}
		Z.D = L.D * R.D;
		truncate_precision(Z.D,trunclen);
		// Z.U = L.U + L.C/L.D * L.P/L.Q * R.S + L.P/L.Q * R.U;
		// Z.V = R.D * R.Q * L.V + R.D * L.C * L.P * R.T + L.D * L.P * R.V;
		Z.V = R.D * (R.Q * L.V + L.C * tmp) + L.D * L.P * R.V;
		truncate_precision(Z.V,trunclen);
		break;
		}
	}
}

}  // namespace cln
