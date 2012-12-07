// Univariate Polynomials over some subring of the numbers.

#include "cln/SV_number.h"
#include "cln/number.h"
#include "cln/integer.h"
#include "cln/exception.h"

namespace cln {

// Assume a ring is a number ring.
  inline cl_heap_number_ring* TheNumberRing (const cl_ring& R)
  { return (cl_heap_number_ring*) R.heappointer; }

// Normalize a vector: remove leading zero coefficients.
// The result vector is known to have length len > 0.
static inline void num_normalize (cl_number_ring_ops<cl_number>& ops, cl_SV_number& result, uintL len)
{
	if (ops.zerop(result[len-1])) {
		len--;
		while (len > 0) {
			if (!ops.zerop(result[len-1]))
				break;
			len--;
		}
		var cl_SV_number newresult = cl_SV_number(cl_make_heap_SV_number_uninit(len));
		for (var sintL i = len-1; i >= 0; i--)
			init1(cl_number, newresult[i]) (result[i]);
		result = newresult;
	}
}

static void num_fprint (cl_heap_univpoly_ring* UPR, std::ostream& stream, const _cl_UP& x)
{{
	DeclarePoly(cl_SV_number,x);
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	var sintL xlen = x.size();
	if (xlen == 0)
		fprint(stream, "0");
	else {
		var const cl_ring& R = UPR->basering();
		var cl_string varname = get_varname(UPR);
		for (var sintL i = xlen-1; i >= 0; i--)
			if (!ops.zerop(x[i])) {
				if (i < xlen-1)
					fprint(stream, " + ");
				fprint(stream, cl_ring_element(R,x[i]));
				if (i > 0) {
					fprint(stream, "*");
					fprint(stream, varname);
					if (i != 1) {
						fprint(stream, "^");
						fprintdecimal(stream, i);
					}
				}
			}
	}
}}

static bool num_equal (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_SV_number,x);
	DeclarePoly(cl_SV_number,y);
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (!(xlen == ylen))
		return false;
	for (var sintL i = xlen-1; i >= 0; i--)
		if (!ops.equal(x[i],y[i]))
			return false;
	return true;
}}

static const _cl_UP num_zero (cl_heap_univpoly_ring* UPR)
{
	return _cl_UP(UPR, cl_null_SV_number);
}

static bool num_zerop (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{
	unused UPR;
 {	DeclarePoly(cl_SV_number,x);
	var sintL xlen = x.size();
	if (xlen == 0)
		return true;
	else
		return false;
}}

static const _cl_UP num_plus (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_SV_number,x);
	DeclarePoly(cl_SV_number,y);
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (xlen == 0)
		return _cl_UP(UPR, y);
	if (ylen == 0)
		return _cl_UP(UPR, x);
	// Now xlen > 0, ylen > 0.
	if (xlen > ylen) {
		var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(xlen));
		var sintL i;
		for (i = xlen-1; i >= ylen; i--)
			init1(cl_number, result[i]) (x[i]);
		for (i = ylen-1; i >= 0; i--)
			init1(cl_number, result[i]) (ops.plus(x[i],y[i]));
		return _cl_UP(UPR, result);
	}
	if (xlen < ylen) {
		var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(ylen));
		var sintL i;
		for (i = ylen-1; i >= xlen; i--)
			init1(cl_number, result[i]) (y[i]);
		for (i = xlen-1; i >= 0; i--)
			init1(cl_number, result[i]) (ops.plus(x[i],y[i]));
		return _cl_UP(UPR, result);
	}
	// Now xlen = ylen > 0. Add and normalize simultaneously.
	for (var sintL i = xlen-1; i >= 0; i--) {
		var cl_number hicoeff = ops.plus(x[i],y[i]);
		if (!ops.zerop(hicoeff)) {
			var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(i+1));
			init1(cl_number, result[i]) (hicoeff);
			for (i-- ; i >= 0; i--)
				init1(cl_number, result[i]) (ops.plus(x[i],y[i]));
			return _cl_UP(UPR, result);
		}
	}
	return _cl_UP(UPR, cl_null_SV_number);
}}

static const _cl_UP num_uminus (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{{
	DeclarePoly(cl_SV_number,x);
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	var sintL xlen = x.size();
	if (xlen == 0)
		return _cl_UP(UPR, x);
	// Now xlen > 0.
	// Negate. No normalization necessary, since the degree doesn't change.
	var sintL i = xlen-1;
	var cl_number hicoeff = ops.uminus(x[i]);
	if (ops.zerop(hicoeff)) throw runtime_exception();
	var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(xlen));
	init1(cl_number, result[i]) (hicoeff);
	for (i-- ; i >= 0; i--)
		init1(cl_number, result[i]) (ops.uminus(x[i]));
	return _cl_UP(UPR, result);
}}

static const _cl_UP num_minus (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_SV_number,x);
	DeclarePoly(cl_SV_number,y);
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (ylen == 0)
		return _cl_UP(UPR, x);
	if (xlen == 0)
		return num_uminus(UPR, _cl_UP(UPR, y));
	// Now xlen > 0, ylen > 0.
	if (xlen > ylen) {
		var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(xlen));
		var sintL i;
		for (i = xlen-1; i >= ylen; i--)
			init1(cl_number, result[i]) (x[i]);
		for (i = ylen-1; i >= 0; i--)
			init1(cl_number, result[i]) (ops.minus(x[i],y[i]));
		return _cl_UP(UPR, result);
	}
	if (xlen < ylen) {
		var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(ylen));
		var sintL i;
		for (i = ylen-1; i >= xlen; i--)
			init1(cl_number, result[i]) (ops.uminus(y[i]));
		for (i = xlen-1; i >= 0; i--)
			init1(cl_number, result[i]) (ops.minus(x[i],y[i]));
		return _cl_UP(UPR, result);
	}
	// Now xlen = ylen > 0. Add and normalize simultaneously.
	for (var sintL i = xlen-1; i >= 0; i--) {
		var cl_number hicoeff = ops.minus(x[i],y[i]);
		if (!ops.zerop(hicoeff)) {
			var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(i+1));
			init1(cl_number, result[i]) (hicoeff);
			for (i-- ; i >= 0; i--)
				init1(cl_number, result[i]) (ops.minus(x[i],y[i]));
			return _cl_UP(UPR, result);
		}
	}
	return _cl_UP(UPR, cl_null_SV_number);
}}

static const _cl_UP num_one (cl_heap_univpoly_ring* UPR)
{
	var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(1));
	init1(cl_number, result[0]) (1);
	return _cl_UP(UPR, result);
}

static const _cl_UP num_canonhom (cl_heap_univpoly_ring* UPR, const cl_I& x)
{
	var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(1));
	init1(cl_number, result[0]) (x);
	return _cl_UP(UPR, result);
}

static const _cl_UP num_mul (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_SV_number,x);
	DeclarePoly(cl_SV_number,y);
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (xlen == 0)
		return _cl_UP(UPR, x);
	if (ylen == 0)
		return _cl_UP(UPR, y);
	// Multiply.
	var sintL len = xlen + ylen - 1;
	var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(len));
	if (xlen < ylen) {
		{
			var sintL i = xlen-1;
			var cl_number xi = x[i];
			for (sintL j = ylen-1; j >= 0; j--)
				init1(cl_number, result[i+j]) (ops.mul(xi,y[j]));
		}
		for (sintL i = xlen-2; i >= 0; i--) {
			var cl_number xi = x[i];
			for (sintL j = ylen-1; j > 0; j--)
				result[i+j] = ops.plus(result[i+j],ops.mul(xi,y[j]));
			/* j=0 */ init1(cl_number, result[i]) (ops.mul(xi,y[0]));
		}
	} else {
		{
			var sintL j = ylen-1;
			var cl_number yj = y[j];
			for (sintL i = xlen-1; i >= 0; i--)
				init1(cl_number, result[i+j]) (ops.mul(x[i],yj));
		}
		for (sintL j = ylen-2; j >= 0; j--) {
			var cl_number yj = y[j];
			for (sintL i = xlen-1; i > 0; i--)
				result[i+j] = ops.plus(result[i+j],ops.mul(x[i],yj));
			/* i=0 */ init1(cl_number, result[j]) (ops.mul(x[0],yj));
		}
	}
	// Normalize (not necessary in integral domains).
	//num_normalize(ops,result,len);
	if (ops.zerop(result[len-1])) throw runtime_exception();
	return _cl_UP(UPR, result);
}}

static const _cl_UP num_square (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{{
	DeclarePoly(cl_SV_number,x);
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	var sintL xlen = x.size();
	if (xlen == 0)
		return cl_UP(UPR, x);
	var sintL len = 2*xlen-1;
	var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(len));
	if (xlen > 1) {
		// Loop through all 0 <= j < i <= xlen-1.
		{
			var sintL i = xlen-1;
			var cl_number xi = x[i];
			for (sintL j = i-1; j >= 0; j--)
				init1(cl_number, result[i+j]) (ops.mul(xi,x[j]));
		}
		{for (sintL i = xlen-2; i >= 1; i--) {
			var cl_number xi = x[i];
			for (sintL j = i-1; j >= 1; j--)
				result[i+j] = ops.plus(result[i+j],ops.mul(xi,x[j]));
			/* j=0 */ init1(cl_number, result[i]) (ops.mul(xi,x[0]));
		}}
		// Double.
		{for (sintL i = len-2; i >= 1; i--)
			result[i] = ops.plus(result[i],result[i]);
		}
		// Add squares.
		init1(cl_number, result[2*(xlen-1)]) (ops.square(x[xlen-1]));
		for (sintL i = xlen-2; i >= 1; i--)
			result[2*i] = ops.plus(result[2*i],ops.square(x[i]));
	}
	init1(cl_number, result[0]) (ops.square(x[0]));
	// Normalize (not necessary in integral domains).
	//num_normalize(ops,result,len);
	if (ops.zerop(result[len-1])) throw runtime_exception();
	return _cl_UP(UPR, result);
}}

static const _cl_UP num_exptpos (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const cl_I& y)
{
	var _cl_UP a = x;
	var cl_I b = y;
	while (!oddp(b)) { a = UPR->_square(a); b = b >> 1; }
	var _cl_UP c = a;
	until (b == 1)
	  { b = b >> 1;
	    a = UPR->_square(a);
	    if (oddp(b)) { c = UPR->_mul(a,c); }
	  }
	return c;
}

static const _cl_UP num_scalmul (cl_heap_univpoly_ring* UPR, const cl_ring_element& x, const _cl_UP& y)
{
	if (!(UPR->basering() == x.ring())) throw runtime_exception();
 {
	DeclarePoly(cl_number,x);
	DeclarePoly(cl_SV_number,y);
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	var sintL ylen = y.size();
	if (ylen == 0)
		return _cl_UP(UPR, y);
	if (ops.zerop(x))
		return _cl_UP(UPR, cl_null_SV_number);
	// Now ylen > 0.
	// No normalization necessary, since the degree doesn't change.
	var cl_SV_number result = cl_SV_number(cl_make_heap_SV_number_uninit(ylen));
	for (sintL i = ylen-1; i >= 0; i--)
		init1(cl_number, result[i]) (ops.mul(x,y[i]));
	return _cl_UP(UPR, result);
}}

static sintL num_degree (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{
	unused UPR;
 {	DeclarePoly(cl_SV_number,x);
	return (sintL) x.size() - 1;
}}

static sintL num_ldegree (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{{
	DeclarePoly(cl_SV_number,x);
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	var sintL xlen = x.size();
	for (sintL i = 0; i < xlen; i++) {
		if (!ops.zerop(x[i]))
			return i;
	}
	return -1;
}}

static const _cl_UP num_monomial (cl_heap_univpoly_ring* UPR, const cl_ring_element& x, uintL e)
{
	if (!(UPR->basering() == x.ring())) throw runtime_exception();
 {	DeclarePoly(cl_number,x);
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	if (ops.zerop(x))
		return _cl_UP(UPR, cl_null_SV_number);
	else {
		var sintL len = e+1;
		var cl_SV_number result = cl_SV_number(len);
		result[e] = x;
		return _cl_UP(UPR, result);
	}
}}

static const cl_ring_element num_coeff (cl_heap_univpoly_ring* UPR, const _cl_UP& x, uintL index)
{{
	DeclarePoly(cl_SV_number,x);
	var cl_heap_number_ring* R = TheNumberRing(UPR->basering());
	if (index < x.size())
		return cl_ring_element(R, x[index]);
	else
		return R->zero();
}}

static const _cl_UP num_create (cl_heap_univpoly_ring* UPR, sintL deg)
{
	if (deg < 0)
		return _cl_UP(UPR, cl_null_SV_number);
	else {
		var sintL len = deg+1;
		return _cl_UP(UPR, cl_SV_number(len));
	}
}

static void num_set_coeff (cl_heap_univpoly_ring* UPR, _cl_UP& x, uintL index, const cl_ring_element& y)
{{
	DeclareMutablePoly(cl_SV_number,x);
	if (!(UPR->basering() == y.ring())) throw runtime_exception();
  {	DeclarePoly(cl_number,y);
	if (!(index < x.size())) throw runtime_exception();
	x[index] = y;
}}}

static void num_finalize (cl_heap_univpoly_ring* UPR, _cl_UP& x)
{{
	DeclareMutablePoly(cl_SV_number,x); // NB: x is modified by reference!
	var cl_number_ring_ops<cl_number>& ops = *TheNumberRing(UPR->basering())->ops;
	var uintL len = x.size();
	if (len > 0)
		num_normalize(ops,x,len);
}}

static const cl_ring_element num_eval (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const cl_ring_element& y)
{{
	// Method:
	// If x = 0, return 0.
	// If y = 0, return x[0].
	// Else compute (...(x[len-1]*y+x[len-2])*y ...)*y + x[0].
	DeclarePoly(cl_SV_number,x);
	if (!(UPR->basering() == y.ring())) throw runtime_exception();
  {	DeclarePoly(cl_number,y);
	var cl_heap_number_ring* R = TheNumberRing(UPR->basering());
	var cl_number_ring_ops<cl_number>& ops = *R->ops;
	var uintL len = x.size();
	if (len==0)
		return R->zero();
	if (ops.zerop(y))
		return cl_ring_element(R, x[0]);
	var sintL i = len-1;
	var cl_number z = x[i];
	for ( ; --i >= 0; )
		z = ops.plus(ops.mul(z,y),x[i]);
	return cl_ring_element(R, z);
}}}

static cl_univpoly_setops num_setops = {
	num_fprint,
	num_equal
};

static cl_univpoly_addops num_addops = {
	num_zero,
	num_zerop,
	num_plus,
	num_minus,
	num_uminus
};

static cl_univpoly_mulops num_mulops = {
	num_one,
	num_canonhom,
	num_mul,
	num_square,
	num_exptpos
};

static cl_univpoly_modulops num_modulops = {
	num_scalmul
};

static cl_univpoly_polyops num_polyops = {
	num_degree,
	num_ldegree,
	num_monomial,
	num_coeff,
	num_create,
	num_set_coeff,
	num_finalize,
	num_eval
};

class cl_heap_num_univpoly_ring : public cl_heap_univpoly_ring {
	SUBCLASS_cl_heap_univpoly_ring()
public:
	// Constructor.
	cl_heap_num_univpoly_ring (const cl_ring& r);
	// Destructor.
	~cl_heap_num_univpoly_ring () {}
};

static void cl_heap_num_univpoly_ring_destructor (cl_heap* pointer)
{
	(*(cl_heap_num_univpoly_ring*)pointer).~cl_heap_num_univpoly_ring();
}

cl_class cl_class_num_univpoly_ring = {
	cl_heap_num_univpoly_ring_destructor,
	cl_class_flags_univpoly_ring
};

// Constructor.
inline cl_heap_num_univpoly_ring::cl_heap_num_univpoly_ring (const cl_ring& r)
	: cl_heap_univpoly_ring (r, &num_setops, &num_addops, &num_mulops, &num_modulops, &num_polyops)
{
	type = &cl_class_num_univpoly_ring;
}

}  // namespace cln
