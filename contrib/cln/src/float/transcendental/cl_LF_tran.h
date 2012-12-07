// cl_LF internals, transcendental functions

#ifndef _CL_LF_TRAN_H
#define _CL_LF_TRAN_H

#include "cln/integer.h"
#include "cln/integer_ring.h"
#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

// Subroutine for evaluating
// sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
// where all the entries are small integers (ideally polynomials in n).
// Some of the factors (a,b,p,q) may be omitted. They are then understood to
// be 1. This is fast because it groups factors together before multiplying.
// Result will be a cl_LF with len digits.

// There are various alternative implementations of the same algorithm that
// differ in the way the series is represented and, as a consequence, in memory
// consumption.
//
// 1st implementation (series is precomputed entirely)
// Arguments:
//   Vectors p[0..N-1], q[0..N-1], a[0..N-1], b[0..N-1], N.
//   Some of the vectors (a,b,p,q) can be a NULL pointer, all of its entries
//   are then understood to be 1.
//   If given, a vector qs[0..N-1] which the evaluation routine may use to
//   split off q[n] into q[n]*2^qs[n]. qs may be NULL, in that case no shift
//   optimizations will be used. (They are worth it only if a significant
//   amount of multiplication work can be saved by shifts.)
//
// 2nd implemenation (series is computed on demand, as a stream)
// In this alternate implementation the series is not represented as a couple
// of arrays, but as a method returning each tuple (p(n),q(n),a(n),b(n))
// in turn. This is preferrable if the a(n) are big, in order to avoid too
// much memory usage at the same time.
// The next() function is called N times and is expected to return
// (p(n),q(n),a(n),b(n)) for n=0..N-1 in that order.
//
// 3rd implemenation (series is computed on demand and truncated early)
// This is like the second implementation, but it coerces the integer factors
// to cl_LF of a given length (trunclen) as soon as the integer factor's size
// exceeds the size to store the cl_LF. For this to make sense, trunclen must
// not be smaller than len. In practice, this can shave off substantially from
// the memory consumption but it also bears a potential for rounding errors.
// A minimum trunclen that guarantees correctness must be evaluated on a
// case-by-case basis.
//
// As a variation, it is sometimes advantageous to factor out from q[0..N-1]
// powers of two and put them back in again at the end of the computation by
// left-shifting. These are distinguished by a boolean template parameter.
// (Some combinations might not be implemented yet.)


// In each of the special cases below, none of (a,b,p,q) can be NULL.

struct cl_pqab_series {
	const cl_I* pv;
	      cl_I* qv;
	const cl_I* av;
	const cl_I* bv;
};
struct cl_pqab_series_term {
	cl_I p;
	cl_I q;
	cl_I a;
	cl_I b;
};
struct cl_pqab_series_stream {
	cl_pqab_series_term (*nextfn)(cl_pqab_series_stream&);
	cl_pqab_series_term next () { return nextfn(*this); }
	// Constructor.
	cl_pqab_series_stream (cl_pqab_series_term (*n)(cl_pqab_series_stream&)) : nextfn (n) {}
};
template<bool shift_q> const cl_LF eval_rational_series (uintC N, const cl_pqab_series& args, uintC len);
template<bool shift_q> const cl_LF eval_rational_series (uintC N, cl_pqab_series_stream& args, uintC len);
template<bool shift_q> const cl_LF eval_rational_series (uintC N, cl_pqab_series_stream& args, uintC len, uintC trunclen);

struct cl_pqb_series {
	const cl_I* pv;
	      cl_I* qv;
	const cl_I* bv;
};
struct cl_pqb_series_term {
	cl_I p;
	cl_I q;
	cl_I b;
};
struct cl_pqb_series_stream {
	cl_pqb_series_term (*nextfn)(cl_pqb_series_stream&);
	cl_pqb_series_term next () { return nextfn(*this); }
	// Constructor.
	cl_pqb_series_stream (cl_pqb_series_term (*n)(cl_pqb_series_stream&)) : nextfn (n) {}
};
template<bool shift_q> const cl_LF eval_rational_series (uintC N, const cl_pqb_series& args, uintC len);
template<bool shift_q> const cl_LF eval_rational_series (uintC N, cl_pqb_series_stream& args, uintC len);
template<bool shift_q> const cl_LF eval_rational_series (uintC N, cl_pqb_series_stream& args, uintC len, uintC trunclen);

struct cl_pqa_series {
	const cl_I* pv;
	      cl_I* qv;
	const cl_I* av;
};
struct cl_pqa_series_term {
	cl_I p;
	cl_I q;
	cl_I a;
};
struct cl_pqa_series_stream {
	cl_pqa_series_term (*nextfn)(cl_pqa_series_stream&);
	cl_pqa_series_term next () { return nextfn(*this); }
	// Constructor.
	cl_pqa_series_stream (cl_pqa_series_term (*n)(cl_pqa_series_stream&)) : nextfn (n) {}
};
template<bool shift_q> const cl_LF eval_rational_series (uintC N, const cl_pqa_series& args, uintC len);
template<bool shift_q> const cl_LF eval_rational_series (uintC N, cl_pqa_series_stream& args, uintC len);
template<bool shift_q> const cl_LF eval_rational_series (uintC N, cl_pqa_series_stream& args, uintC len, uintC trunclen);

struct cl_pq_series {
	const cl_I* pv;
	      cl_I* qv;
};
struct cl_pq_series_term {
	cl_I p;
	cl_I q;
};
struct cl_pq_series_stream {
	cl_pq_series_term (*nextfn)(cl_pq_series_stream&);
	cl_pq_series_term next () { return nextfn(*this); }
	// Constructor.
	cl_pq_series_stream (cl_pq_series_term (*n)(cl_pq_series_stream&)) : nextfn (n) {}
};
template<bool shift_q> const cl_LF eval_rational_series (uintC N, const cl_pq_series& args, uintC len);
template<bool shift_q> const cl_LF eval_rational_series (uintC N, cl_pq_series_stream& args, uintC len);
template<bool shift_q> const cl_LF eval_rational_series (uintC N, cl_pq_series_stream& args, uintC len, uintC trunclen);

struct cl_pab_series {
	const cl_I* pv;
	const cl_I* av;
	const cl_I* bv;
};
const cl_LF eval_rational_series (uintC N, const cl_pab_series& args, uintC len);

struct cl_pb_series {
	const cl_I* pv;
	const cl_I* bv;
};
const cl_LF eval_rational_series (uintC N, const cl_pb_series& args, uintC len);

struct cl_pa_series {
	const cl_I* pv;
	const cl_I* av;
};
const cl_LF eval_rational_series (uintC N, const cl_pa_series& args, uintC len);

struct cl_p_series {
	const cl_I* pv;
};
const cl_LF eval_rational_series (uintC N, const cl_p_series& args, uintC len);

struct cl_qab_series {
	      cl_I* qv;
	const cl_I* av;
	const cl_I* bv;
};
template<bool shift_q> const cl_LF eval_rational_series (uintC N, const cl_qab_series& args, uintC len);

struct cl_qb_series {
	      cl_I* qv;
	const cl_I* bv;
};
struct cl_qb_series_term {
	cl_I q;
	cl_I b;
};
struct cl_qb_series_stream {
	cl_qb_series_term (*nextfn)(cl_qb_series_stream&);
	cl_qb_series_term next () { return nextfn(*this); }
	// Constructor.
	cl_qb_series_stream (cl_qb_series_term (*n)(cl_qb_series_stream&)) : nextfn (n) {}
};
template<bool shift_q> const cl_LF eval_rational_series (uintC N, const cl_qb_series& args, uintC len);
template<bool shift_q> const cl_LF eval_rational_series (uintC N, cl_qb_series_stream& args, uintC len);

struct cl_qa_series {
	      cl_I* qv;
	const cl_I* av;
};
template<bool shift_q> const cl_LF eval_rational_series (uintC N, const cl_qa_series& args, uintC len);

struct cl_q_series {
	      cl_I* qv;
};
struct cl_q_series_term {
	cl_I q;
};
struct cl_q_series_stream {
	cl_q_series_term (*nextfn)(cl_q_series_stream&);
	cl_q_series_term next () { return nextfn(*this); }
	// Constructor.
	cl_q_series_stream (cl_q_series_term (*n)(cl_q_series_stream&)) : nextfn (n) {}
};
template<bool shift_q> const cl_LF eval_rational_series (uintC N, const cl_q_series& args, uintC len);
template<bool shift_q> const cl_LF eval_rational_series (uintC N, cl_q_series_stream& args, uintC len);

struct cl_ab_series {
	const cl_I* av;
	const cl_I* bv;
};
const cl_LF eval_rational_series (uintC N, const cl_ab_series& args, uintC len);

struct cl_b_series {
	const cl_I* bv;
};
const cl_LF eval_rational_series (uintC N, const cl_b_series& args, uintC len);

struct cl_a_series {
	const cl_I* av;
};
const cl_LF eval_rational_series (uintC N, const cl_a_series& args, uintC len);

struct cl__series {
};
const cl_LF eval_rational_series (uintC N, const cl__series& args, uintC len);


// [Generalization.]
// Subroutine:
// Evaluates S = sum(N1 <= n < N2, (p(N1)...p(n))/(q(N1)...q(n)))
// and U = sum(N1 <= n < N2,
//             (c(N1)/d(N1)+...+c(n)/d(n))*(p(N1)...p(n))/(q(N1)...q(n)))
// and returns
//     P = p(N1)...p(N2-1),
//     Q = q(N1)...q(N2-1),
//     T = Q*S,
//     C/D = c(N1)/d(N1)+...+c(N2-1)/d(N2-1),
//     V = D*Q*U,
// all integers. On entry N1 < N2.
struct cl_pqcd_series_term {
	cl_I p;
	cl_I q;
	cl_I c;
	cl_I d;
};
template<class cl_T>
struct cl_pqcd_series_result {
	cl_T P;
	cl_T Q;
	cl_T T;
	cl_T C;
	cl_T D;
	cl_T V;
};
struct cl_pqcd_series_stream {
	cl_pqcd_series_term (*nextfn)(cl_pqcd_series_stream&);
	cl_pqcd_series_term next () { return nextfn(*this); }
	// Constructor.
	cl_pqcd_series_stream( cl_pqcd_series_term (*n)(cl_pqcd_series_stream&)) : nextfn (n) {}
};
void eval_pqcd_series_aux (uintC N, cl_pqcd_series_term* args, cl_pqcd_series_result<cl_I>& Z, bool rightmost = true);
void eval_pqcd_series_aux (uintC N, cl_pqcd_series_stream& args, cl_pqcd_series_result<cl_I>& Z, bool rightmost = true);
void eval_pqcd_series_aux (uintC N, cl_pqcd_series_stream& args, cl_pqcd_series_result<cl_R>& Z, uintC trunclen, bool rightmost = true);
// Ditto, but returns U/S.
const cl_LF eval_pqcd_series (uintC N, cl_pqcd_series_term* args, uintC len);
const cl_LF eval_pqcd_series (uintC N, cl_pqcd_series_stream& args, uintC len);
const cl_LF eval_pqcd_series (uintC N, cl_pqcd_series_stream& args, uintC len, uintC trunclen);

// [Special case c(n)=1.]
// Subroutine:
// Evaluates S = sum(N1 <= n < N2, (p(N1)...p(n))/(q(N1)...q(n)))
// and U = sum(N1 <= n < N2, (1/d(N1)+...+1/d(n))*(p(N1)...p(n))/(q(N1)...q(n)))
// and returns
//     P = p(N1)...p(N2-1),
//     Q = q(N1)...q(N2-1),
//     T = Q*S,
//     C/D = 1/d(N1)+...+1/d(N2-1),
//     V = D*Q*U,
// all integers. On entry N1 < N2.
struct cl_pqd_series_term {
	cl_I p;
	cl_I q;
	cl_I d;
};
template<class cl_T>
struct cl_pqd_series_result {
	cl_T P;
	cl_T Q;
	cl_T T;
	cl_T C;
	cl_T D;
	cl_T V;
};
struct cl_pqd_series_stream {
	cl_pqd_series_term (*nextfn)(cl_pqd_series_stream&);
	cl_pqd_series_term next () { return nextfn(*this); }
	// Constructor.
	cl_pqd_series_stream( cl_pqd_series_term (*n)(cl_pqd_series_stream&)) : nextfn (n) {}
};
void eval_pqd_series_aux (uintC N, cl_pqd_series_term* args, cl_pqd_series_result<cl_I>& Z, bool rightmost = true);
void eval_pqd_series_aux (uintC N, cl_pqd_series_stream& args, cl_pqd_series_result<cl_I>& Z, bool rightmost = true);
void eval_pqd_series_aux (uintC N, cl_pqd_series_stream& args, cl_pqd_series_result<cl_R>& Z, uintC trunclen, bool rightmost = true);
// Ditto, but returns U/S.
const cl_LF eval_pqd_series (uintC N, cl_pqd_series_term* args, uintC len);
const cl_LF eval_pqd_series (uintC N, cl_pqd_series_stream& args, uintC len);
const cl_LF eval_pqd_series (uintC N, cl_pqd_series_stream& args, uintC len, uintC trunclen);

// Helper function to divide q by the largest s such that 2^s divides q, returns s.
inline uintC
pullout_shiftcount(cl_I& q)
{
	var uintC qs = 0;
	if (!zerop(q)) {
		qs = ord2(q);
		if (qs > 0)
			q = q >> qs;
	}
	return qs;
}

// Helper function to convert integer of length > trunclen to long float of
// length = trunclen.
inline void
truncate_precision(cl_R& x, uintC trunclen)
{
	if (instanceof(x,cl_I_ring) &&
	    integer_length(the<cl_I>(x))>trunclen*intDsize) {
		x = cl_I_to_LF(the<cl_I>(x),trunclen);
	}
}

}  // namespace cln

#endif /* _CL_LF_TRAN_H */
