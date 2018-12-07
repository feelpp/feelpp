// Univariate Polynomials over a general ring.

#include "cln/SV_ringelt.h"
#include "cln/integer.h"
#include "cln/exception.h"

namespace cln {

// Assume a ring is a ring.
  inline cl_heap_ring* TheRing (const cl_ring& R)
  { return (cl_heap_ring*) R.heappointer; }

// Normalize a vector: remove leading zero coefficients.
// The result vector is known to have length len > 0.
static inline void gen_normalize (cl_heap_ring* R, cl_SV_ringelt& result, uintL len)
{
	if (R->_zerop(result[len-1])) {
		len--;
		while (len > 0) {
			if (!R->_zerop(result[len-1]))
				break;
			len--;
		}
		var cl_SV_ringelt newresult = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(len));
		for (var sintL i = len-1; i >= 0; i--)
			init1(_cl_ring_element, newresult[i]) (result[i]);
		result = newresult;
	}
}

static void gen_fprint (cl_heap_univpoly_ring* UPR, std::ostream& stream, const _cl_UP& x)
{{
	DeclarePoly(cl_SV_ringelt,x);
	var cl_heap_ring* R = TheRing(UPR->basering());
	var sintL xlen = x.size();
	if (xlen == 0)
		fprint(stream, "0");
	else {
		var cl_string varname = get_varname(UPR);
		for (var sintL i = xlen-1; i >= 0; i--)
			if (!R->_zerop(x[i])) {
				if (i < xlen-1)
					fprint(stream, " + ");
				fprint(stream, "(");
				R->_fprint(stream, x[i]);
				fprint(stream, ")");
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

static bool gen_equal (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_SV_ringelt,x);
	DeclarePoly(cl_SV_ringelt,y);
	var cl_heap_ring* R = TheRing(UPR->basering());
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (!(xlen == ylen))
		return false;
	for (var sintL i = xlen-1; i >= 0; i--)
		if (!R->_equal(x[i],y[i]))
			return false;
	return true;
}}

static const _cl_UP gen_zero (cl_heap_univpoly_ring* UPR)
{
	return _cl_UP(UPR, cl_null_SV_ringelt);
}

static bool gen_zerop (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{
	unused UPR;
 {	DeclarePoly(cl_SV_ringelt,x);
	var sintL xlen = x.size();
	if (xlen == 0)
		return true;
	else
		return false;
}}

static const _cl_UP gen_plus (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_SV_ringelt,x);
	DeclarePoly(cl_SV_ringelt,y);
	var cl_heap_ring* R = TheRing(UPR->basering());
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (xlen == 0)
		return _cl_UP(UPR, y);
	if (ylen == 0)
		return _cl_UP(UPR, x);
	// Now xlen > 0, ylen > 0.
	if (xlen > ylen) {
		var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(xlen));
		var sintL i;
		for (i = xlen-1; i >= ylen; i--)
			init1(_cl_ring_element, result[i]) (x[i]);
		for (i = ylen-1; i >= 0; i--)
			init1(_cl_ring_element, result[i]) (R->_plus(x[i],y[i]));
		return _cl_UP(UPR, result);
	}
	if (xlen < ylen) {
		var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(ylen));
		var sintL i;
		for (i = ylen-1; i >= xlen; i--)
			init1(_cl_ring_element, result[i]) (y[i]);
		for (i = xlen-1; i >= 0; i--)
			init1(_cl_ring_element, result[i]) (R->_plus(x[i],y[i]));
		return _cl_UP(UPR, result);
	}
	// Now xlen = ylen > 0. Add and normalize simultaneously.
	for (var sintL i = xlen-1; i >= 0; i--) {
		var _cl_ring_element hicoeff = R->_plus(x[i],y[i]);
		if (!R->_zerop(hicoeff)) {
			var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(i+1));
			init1(_cl_ring_element, result[i]) (hicoeff);
			for (i-- ; i >= 0; i--)
				init1(_cl_ring_element, result[i]) (R->_plus(x[i],y[i]));
			return _cl_UP(UPR, result);
		}
	}
	return _cl_UP(UPR, cl_null_SV_ringelt);
}}

static const _cl_UP gen_uminus (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{{
	DeclarePoly(cl_SV_ringelt,x);
	var cl_heap_ring* R = TheRing(UPR->basering());
	var sintL xlen = x.size();
	if (xlen == 0)
		return _cl_UP(UPR, x);
	// Now xlen > 0.
	// Negate. No normalization necessary, since the degree doesn't change.
	var sintL i = xlen-1;
	var _cl_ring_element hicoeff = R->_uminus(x[i]);
	if (R->_zerop(hicoeff)) throw runtime_exception();
	var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(xlen));
	init1(_cl_ring_element, result[i]) (hicoeff);
	for (i-- ; i >= 0; i--)
		init1(_cl_ring_element, result[i]) (R->_uminus(x[i]));
	return _cl_UP(UPR, result);
}}

static const _cl_UP gen_minus (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_SV_ringelt,x);
	DeclarePoly(cl_SV_ringelt,y);
	var cl_heap_ring* R = TheRing(UPR->basering());
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (ylen == 0)
		return _cl_UP(UPR, x);
	if (xlen == 0)
		return gen_uminus(UPR,_cl_UP(UPR, y));
	// Now xlen > 0, ylen > 0.
	if (xlen > ylen) {
		var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(xlen));
		var sintL i;
		for (i = xlen-1; i >= ylen; i--)
			init1(_cl_ring_element, result[i]) (x[i]);
		for (i = ylen-1; i >= 0; i--)
			init1(_cl_ring_element, result[i]) (R->_minus(x[i],y[i]));
		return _cl_UP(UPR, result);
	}
	if (xlen < ylen) {
		var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(ylen));
		var sintL i;
		for (i = ylen-1; i >= xlen; i--)
			init1(_cl_ring_element, result[i]) (R->_uminus(y[i]));
		for (i = xlen-1; i >= 0; i--)
			init1(_cl_ring_element, result[i]) (R->_minus(x[i],y[i]));
		return _cl_UP(UPR, result);
	}
	// Now xlen = ylen > 0. Add and normalize simultaneously.
	for (var sintL i = xlen-1; i >= 0; i--) {
		var _cl_ring_element hicoeff = R->_minus(x[i],y[i]);
		if (!R->_zerop(hicoeff)) {
			var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(i+1));
			init1(_cl_ring_element, result[i]) (hicoeff);
			for (i-- ; i >= 0; i--)
				init1(_cl_ring_element, result[i]) (R->_minus(x[i],y[i]));
			return _cl_UP(UPR, result);
		}
	}
	return _cl_UP(UPR, cl_null_SV_ringelt);
}}

static const _cl_UP gen_one (cl_heap_univpoly_ring* UPR)
{
	var cl_heap_ring* R = TheRing(UPR->basering());
	var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(1));
	init1(_cl_ring_element, result[0]) (R->_one());
	return _cl_UP(UPR, result);
}

static const _cl_UP gen_canonhom (cl_heap_univpoly_ring* UPR, const cl_I& x)
{
	var cl_heap_ring* R = TheRing(UPR->basering());
	var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(1));
	init1(_cl_ring_element, result[0]) (R->_canonhom(x));
	return _cl_UP(UPR, result);
}

static const _cl_UP gen_mul (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_SV_ringelt,x);
	DeclarePoly(cl_SV_ringelt,y);
	var cl_heap_ring* R = TheRing(UPR->basering());
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (xlen == 0)
		return _cl_UP(UPR, x);
	if (ylen == 0)
		return _cl_UP(UPR, y);
	// Multiply.
	var sintL len = xlen + ylen - 1;
	var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(len));
	if (xlen < ylen) {
		{
			var sintL i = xlen-1;
			var _cl_ring_element xi = x[i];
			for (sintL j = ylen-1; j >= 0; j--)
				init1(_cl_ring_element, result[i+j]) (R->_mul(xi,y[j]));
		}
		for (sintL i = xlen-2; i >= 0; i--) {
			var _cl_ring_element xi = x[i];
			for (sintL j = ylen-1; j > 0; j--)
				result[i+j] = R->_plus(result[i+j],R->_mul(xi,y[j]));
			/* j=0 */ init1(_cl_ring_element, result[i]) (R->_mul(xi,y[0]));
		}
	} else {
		{
			var sintL j = ylen-1;
			var _cl_ring_element yj = y[j];
			for (sintL i = xlen-1; i >= 0; i--)
				init1(_cl_ring_element, result[i+j]) (R->_mul(x[i],yj));
		}
		for (sintL j = ylen-2; j >= 0; j--) {
			var _cl_ring_element yj = y[j];
			for (sintL i = xlen-1; i > 0; i--)
				result[i+j] = R->_plus(result[i+j],R->_mul(x[i],yj));
			/* i=0 */ init1(_cl_ring_element, result[j]) (R->_mul(x[0],yj));
		}
	}
	// Normalize (not necessary in integral domains).
	//gen_normalize(R,result,len);
	if (R->_zerop(result[len-1])) throw runtime_exception();
	return _cl_UP(UPR, result);
}}

static const _cl_UP gen_square (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{{
	DeclarePoly(cl_SV_ringelt,x);
	var cl_heap_ring* R = TheRing(UPR->basering());
	var sintL xlen = x.size();
	if (xlen == 0)
		return cl_UP(UPR, x);
	var sintL len = 2*xlen-1;
	var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(len));
	if (xlen > 1) {
		// Loop through all 0 <= j < i <= xlen-1.
		{
			var sintL i = xlen-1;
			var _cl_ring_element xi = x[i];
			for (sintL j = i-1; j >= 0; j--)
				init1(_cl_ring_element, result[i+j]) (R->_mul(xi,x[j]));
		}
		{for (sintL i = xlen-2; i >= 1; i--) {
			var _cl_ring_element xi = x[i];
			for (sintL j = i-1; j >= 1; j--)
				result[i+j] = R->_plus(result[i+j],R->_mul(xi,x[j]));
			/* j=0 */ init1(_cl_ring_element, result[i]) (R->_mul(xi,x[0]));
		}}
		// Double.
		{for (sintL i = len-2; i >= 1; i--)
			result[i] = R->_plus(result[i],result[i]);
		}
		// Add squares.
		init1(_cl_ring_element, result[2*(xlen-1)]) (R->_square(x[xlen-1]));
		for (sintL i = xlen-2; i >= 1; i--)
			result[2*i] = R->_plus(result[2*i],R->_square(x[i]));
	}
	init1(_cl_ring_element, result[0]) (R->_square(x[0]));
	// Normalize (not necessary in integral domains).
	//gen_normalize(R,result,len);
	if (R->_zerop(result[len-1])) throw runtime_exception();
	return _cl_UP(UPR, result);
}}

static const _cl_UP gen_exptpos (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const cl_I& y)
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

static const _cl_UP gen_scalmul (cl_heap_univpoly_ring* UPR, const cl_ring_element& x, const _cl_UP& y)
{
	if (!(UPR->basering() == x.ring())) throw runtime_exception();
 {
	DeclarePoly(cl_SV_ringelt,y);
	var cl_heap_ring* R = TheRing(UPR->basering());
	var sintL ylen = y.size();
	if (ylen == 0)
		return _cl_UP(UPR, y);
	if (R->zerop(x))
		return _cl_UP(UPR, cl_null_SV_ringelt);
	var cl_SV_ringelt result = cl_SV_ringelt(cl_make_heap_SV_ringelt_uninit(ylen));
	for (sintL i = ylen-1; i >= 0; i--)
		init1(_cl_ring_element, result[i]) (R->_mul(x,y[i]));
	// Normalize (not necessary in integral domains).
	//gen_normalize(R,result,ylen);
	if (R->_zerop(result[ylen-1])) throw runtime_exception();
	return _cl_UP(UPR, result);
}}

static sintL gen_degree (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{
	unused UPR;
 {	DeclarePoly(cl_SV_ringelt,x);
	return (sintL) x.size() - 1;
}}

static sintL gen_ldegree (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{{	DeclarePoly(cl_SV_ringelt,x);
	var cl_heap_ring* R = TheRing(UPR->basering());
	var sintL xlen = x.size();
	for (sintL i = 0; i < xlen; i++) {
		if (!R->_zerop(x[i]))
			return i;
	}
	return -1;
}}

static const _cl_UP gen_monomial (cl_heap_univpoly_ring* UPR, const cl_ring_element& x, uintL e)
{
	if (!(UPR->basering() == x.ring())) throw runtime_exception();
	var cl_heap_ring* R = TheRing(UPR->basering());
	if (R->_zerop(x))
		return _cl_UP(UPR, cl_null_SV_ringelt);
	else {
		var sintL len = e+1;
		var cl_SV_ringelt result = cl_SV_ringelt(len);
		result[e] = x;
		return _cl_UP(UPR, result);
	}
}

static const cl_ring_element gen_coeff (cl_heap_univpoly_ring* UPR, const _cl_UP& x, uintL index)
{{
	DeclarePoly(cl_SV_ringelt,x);
	var cl_heap_ring* R = TheRing(UPR->basering());
	if (index < x.size())
		return cl_ring_element(R, x[index]);
	else
		return R->zero();
}}

static const _cl_UP gen_create (cl_heap_univpoly_ring* UPR, sintL deg)
{
	if (deg < 0)
		return _cl_UP(UPR, cl_null_SV_ringelt);
	else {
		var sintL len = deg+1;
		return _cl_UP(UPR, cl_SV_ringelt(len));
	}
}

static void gen_set_coeff (cl_heap_univpoly_ring* UPR, _cl_UP& x, uintL index, const cl_ring_element& y)
{{
	DeclareMutablePoly(cl_SV_ringelt,x);
	if (!(UPR->basering() == y.ring())) throw runtime_exception();
	if (!(index < x.size())) throw runtime_exception();
	x[index] = y;
}}

static void gen_finalize (cl_heap_univpoly_ring* UPR, _cl_UP& x)
{{
	DeclareMutablePoly(cl_SV_ringelt,x); // NB: x is modified by reference!
	var cl_heap_ring* R = TheRing(UPR->basering());
	var uintL len = x.size();
	if (len > 0)
		gen_normalize(R,x,len);
}}

static const cl_ring_element gen_eval (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const cl_ring_element& y)
{{
	// Method:
	// If x = 0, return 0.
	// If y = 0, return x[0].
	// Else compute (...(x[len-1]*y+x[len-2])*y ...)*y + x[0].
	DeclarePoly(cl_SV_ringelt,x);
	var cl_heap_ring* R = TheRing(UPR->basering());
	if (!(y.ring() == R)) throw runtime_exception();
	var uintL len = x.size();
	if (len==0)
		return R->zero();
	if (R->_zerop(y))
		return cl_ring_element(R, x[0]);
	var sintL i = len-1;
	var _cl_ring_element z = x[i];
	for ( ; --i >= 0; )
		z = R->_plus(R->_mul(z,y),x[i]);
	return cl_ring_element(R, z);
}}

static cl_univpoly_setops gen_setops = {
	gen_fprint,
	gen_equal
};

static cl_univpoly_addops gen_addops = {
	gen_zero,
	gen_zerop,
	gen_plus,
	gen_minus,
	gen_uminus
};

static cl_univpoly_mulops gen_mulops = {
	gen_one,
	gen_canonhom,
	gen_mul,
	gen_square,
	gen_exptpos
};

static cl_univpoly_modulops gen_modulops = {
	gen_scalmul
};

static cl_univpoly_polyops gen_polyops = {
	gen_degree,
	gen_ldegree,
	gen_monomial,
	gen_coeff,
	gen_create,
	gen_set_coeff,
	gen_finalize,
	gen_eval
};

class cl_heap_gen_univpoly_ring : public cl_heap_univpoly_ring {
	SUBCLASS_cl_heap_univpoly_ring()
public:
	// Constructor.
	cl_heap_gen_univpoly_ring (const cl_ring& r);
	// Destructor
	~cl_heap_gen_univpoly_ring () {}
};

static void cl_heap_gen_univpoly_ring_destructor (cl_heap* pointer)
{
	(*(cl_heap_gen_univpoly_ring*)pointer).~cl_heap_gen_univpoly_ring();
}

cl_class cl_class_gen_univpoly_ring = {
	cl_heap_gen_univpoly_ring_destructor,
	cl_class_flags_univpoly_ring
};

// Constructor.
inline cl_heap_gen_univpoly_ring::cl_heap_gen_univpoly_ring (const cl_ring& r)
	: cl_heap_univpoly_ring (r, &gen_setops, &gen_addops, &gen_mulops, &gen_modulops, &gen_polyops)
{
	type = &cl_class_gen_univpoly_ring;
}

}  // namespace cln
