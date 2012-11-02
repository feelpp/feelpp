// Univariate Polynomials over a ring of modular integers.

#include "cln/GV_modinteger.h"
#include "cln/modinteger.h"
#include "cln/exception.h"

namespace cln {

// Assume a ring is a modint ring.
  inline cl_heap_modint_ring* TheModintRing (const cl_ring& R)
  { return (cl_heap_modint_ring*) R.heappointer; }

// Normalize a vector: remove leading zero coefficients.
// The result vector is known to have length len > 0.
static inline void modint_normalize (cl_heap_modint_ring* R, cl_GV_MI& result, uintL len)
{
	if (R->_zerop(result[len-1])) {
		len--;
		while (len > 0) {
			if (!R->_zerop(result[len-1]))
				break;
			len--;
		}
		var cl_GV_MI newresult = cl_GV_MI(len,R);
		#if 0
		for (var sintL i = len-1; i >= 0; i--)
			newresult[i] = result[i];
		#else
		cl_GV_MI::copy_elements(result,0,newresult,0,len);
		#endif
		result = newresult;
	}
}

static void modint_fprint (cl_heap_univpoly_ring* UPR, std::ostream& stream, const _cl_UP& x)
{{
	DeclarePoly(cl_GV_MI,x);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
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

static bool modint_equal (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_GV_MI,x);
	DeclarePoly(cl_GV_MI,y);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (!(xlen == ylen))
		return false;
	for (var sintL i = xlen-1; i >= 0; i--)
		if (!R->_equal(x[i],y[i]))
			return false;
	return true;
}}

static const _cl_UP modint_zero (cl_heap_univpoly_ring* UPR)
{
	return _cl_UP(UPR, cl_null_GV_I);
}

static bool modint_zerop (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{
	unused UPR;
 {	DeclarePoly(cl_GV_MI,x);
	var sintL xlen = x.size();
	if (xlen == 0)
		return true;
	else
		return false;
}}

static const _cl_UP modint_plus (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_GV_MI,x);
	DeclarePoly(cl_GV_MI,y);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (xlen == 0)
		return _cl_UP(UPR, y);
	if (ylen == 0)
		return _cl_UP(UPR, x);
	// Now xlen > 0, ylen > 0.
	if (xlen > ylen) {
		var cl_GV_MI result = cl_GV_MI(xlen,R);
		var sintL i;
		#if 0
		for (i = xlen-1; i >= ylen; i--)
			result[i] = x[i];
		#else
		cl_GV_MI::copy_elements(x,ylen,result,ylen,xlen-ylen);
		#endif
		for (i = ylen-1; i >= 0; i--)
			result[i] = R->_plus(x[i],y[i]);
		return _cl_UP(UPR, result);
	}
	if (xlen < ylen) {
		var cl_GV_MI result = cl_GV_MI(ylen,R);
		var sintL i;
		#if 0
		for (i = ylen-1; i >= xlen; i--)
			result[i] = y[i];
		#else
		cl_GV_MI::copy_elements(y,xlen,result,xlen,ylen-xlen);
		#endif
		for (i = xlen-1; i >= 0; i--)
			result[i] = R->_plus(x[i],y[i]);
		return _cl_UP(UPR, result);
	}
	// Now xlen = ylen > 0. Add and normalize simultaneously.
	for (var sintL i = xlen-1; i >= 0; i--) {
		var _cl_MI hicoeff = R->_plus(x[i],y[i]);
		if (!R->_zerop(hicoeff)) {
			var cl_GV_MI result = cl_GV_MI(i+1,R);
			result[i] = hicoeff;
			for (i-- ; i >= 0; i--)
				result[i] = R->_plus(x[i],y[i]);
			return _cl_UP(UPR, result);
		}
	}
	return _cl_UP(UPR, cl_null_GV_I);
}}

static const _cl_UP modint_uminus (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{{
	DeclarePoly(cl_GV_MI,x);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var sintL xlen = x.size();
	if (xlen == 0)
		return _cl_UP(UPR, x);
	// Now xlen > 0.
	// Negate. No normalization necessary, since the degree doesn't change.
	var sintL i = xlen-1;
	var _cl_MI hicoeff = R->_uminus(x[i]);
	if (R->_zerop(hicoeff)) throw runtime_exception();
	var cl_GV_MI result = cl_GV_MI(xlen,R);
	result[i] = hicoeff;
	for (i-- ; i >= 0; i--)
		result[i] = R->_uminus(x[i]);
	return _cl_UP(UPR, result);
}}

static const _cl_UP modint_minus (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_GV_MI,x);
	DeclarePoly(cl_GV_MI,y);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (ylen == 0)
		return _cl_UP(UPR, x);
	if (xlen == 0)
		return modint_uminus(UPR, _cl_UP(UPR, y));
	// Now xlen > 0, ylen > 0.
	if (xlen > ylen) {
		var cl_GV_MI result = cl_GV_MI(xlen,R);
		var sintL i;
		#if 0
		for (i = xlen-1; i >= ylen; i--)
			result[i] = x[i];
		#else
		cl_GV_MI::copy_elements(x,ylen,result,ylen,xlen-ylen);
		#endif
		for (i = ylen-1; i >= 0; i--)
			result[i] = R->_minus(x[i],y[i]);
		return _cl_UP(UPR, result);
	}
	if (xlen < ylen) {
		var cl_GV_MI result = cl_GV_MI(ylen,R);
		var sintL i;
		for (i = ylen-1; i >= xlen; i--)
			result[i] = R->_uminus(y[i]);
		for (i = xlen-1; i >= 0; i--)
			result[i] = R->_minus(x[i],y[i]);
		return _cl_UP(UPR, result);
	}
	// Now xlen = ylen > 0. Add and normalize simultaneously.
	for (var sintL i = xlen-1; i >= 0; i--) {
		var _cl_MI hicoeff = R->_minus(x[i],y[i]);
		if (!R->_zerop(hicoeff)) {
			var cl_GV_MI result = cl_GV_MI(i+1,R);
			result[i] = hicoeff;
			for (i-- ; i >= 0; i--)
				result[i] = R->_minus(x[i],y[i]);
			return _cl_UP(UPR, result);
		}
	}
	return _cl_UP(UPR, cl_null_GV_I);
}}

static const _cl_UP modint_one (cl_heap_univpoly_ring* UPR)
{
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var cl_GV_MI result = cl_GV_MI(1,R);
	result[0] = R->_one();
	return _cl_UP(UPR, result);
}

static const _cl_UP modint_canonhom (cl_heap_univpoly_ring* UPR, const cl_I& x)
{
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var cl_GV_MI result = cl_GV_MI(1,R);
	result[0] = R->_canonhom(x);
	return _cl_UP(UPR, result);
}

static const _cl_UP modint_mul (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const _cl_UP& y)
{{
	DeclarePoly(cl_GV_MI,x);
	DeclarePoly(cl_GV_MI,y);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var sintL xlen = x.size();
	var sintL ylen = y.size();
	if (xlen == 0)
		return _cl_UP(UPR, x);
	if (ylen == 0)
		return _cl_UP(UPR, y);
	// Multiply.
	var sintL len = xlen + ylen - 1;
	var cl_GV_MI result = cl_GV_MI(len,R);
	if (xlen < ylen) {
		{
			var sintL i = xlen-1;
			var _cl_MI xi = x[i];
			for (sintL j = ylen-1; j >= 0; j--)
				result[i+j] = R->_mul(xi,y[j]);
		}
		for (sintL i = xlen-2; i >= 0; i--) {
			var _cl_MI xi = x[i];
			for (sintL j = ylen-1; j > 0; j--)
				result[i+j] = R->_plus(result[i+j],R->_mul(xi,y[j]));
			/* j=0 */ result[i] = R->_mul(xi,y[0]);
		}
	} else {
		{
			var sintL j = ylen-1;
			var _cl_MI yj = y[j];
			for (sintL i = xlen-1; i >= 0; i--)
				result[i+j] = R->_mul(x[i],yj);
		}
		for (sintL j = ylen-2; j >= 0; j--) {
			var _cl_MI yj = y[j];
			for (sintL i = xlen-1; i > 0; i--)
				result[i+j] = R->_plus(result[i+j],R->_mul(x[i],yj));
			/* i=0 */ result[j] = R->_mul(x[0],yj);
		}
	}
	// Normalize (not necessary in integral domains).
	//modint_normalize(R,result,len);
	if (R->_zerop(result[len-1])) throw runtime_exception();
	return _cl_UP(UPR, result);
}}

static const _cl_UP modint_square (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{{
	DeclarePoly(cl_GV_MI,x);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var sintL xlen = x.size();
	if (xlen == 0)
		return cl_UP(UPR, x);
	var sintL len = 2*xlen-1;
	var cl_GV_MI result = cl_GV_MI(len,R);
	if (xlen > 1) {
		// Loop through all 0 <= j < i <= xlen-1.
		{
			var sintL i = xlen-1;
			var _cl_MI xi = x[i];
			for (sintL j = i-1; j >= 0; j--)
				result[i+j] = R->_mul(xi,x[j]);
		}
		{for (sintL i = xlen-2; i >= 1; i--) {
			var _cl_MI xi = x[i];
			for (sintL j = i-1; j >= 1; j--)
				result[i+j] = R->_plus(result[i+j],R->_mul(xi,x[j]));
			/* j=0 */ result[i] = R->_mul(xi,x[0]);
		}}
		// Double.
		{for (sintL i = len-2; i >= 1; i--)
			result[i] = R->_plus(result[i],result[i]);
		}
		// Add squares.
		result[2*(xlen-1)] = R->_square(x[xlen-1]);
		for (sintL i = xlen-2; i >= 1; i--)
			result[2*i] = R->_plus(result[2*i],R->_square(x[i]));
	}
	result[0] = R->_square(x[0]);
	// Normalize (not necessary in integral domains).
	//modint_normalize(R,result,len);
	if (R->_zerop(result[len-1])) throw runtime_exception();
	return _cl_UP(UPR, result);
}}

static const _cl_UP modint_exptpos (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const cl_I& y)
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

static const _cl_UP modint_scalmul (cl_heap_univpoly_ring* UPR, const cl_ring_element& x, const _cl_UP& y)
{
	if (!(UPR->basering() == x.ring())) throw runtime_exception();
 {
	DeclarePoly(_cl_MI,x);
	DeclarePoly(cl_GV_MI,y);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var sintL ylen = y.size();
	if (ylen == 0)
		return _cl_UP(UPR, y);
	if (R->_zerop(x))
		return _cl_UP(UPR, cl_null_GV_I);
	// Now ylen > 0.
	// No normalization necessary, since the degree doesn't change.
	var cl_GV_MI result = cl_GV_MI(ylen,R);
	for (sintL i = ylen-1; i >= 0; i--)
		result[i] = R->_mul(x,y[i]);
	return _cl_UP(UPR, result);
}}

static sintL modint_degree (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{
	unused UPR;
 {	DeclarePoly(cl_GV_MI,x);
	return (sintL) x.size() - 1;
}}

static sintL modint_ldegree (cl_heap_univpoly_ring* UPR, const _cl_UP& x)
{{
	DeclarePoly(cl_GV_MI,x);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var sintL xlen = x.size();
	for (sintL i = 0; i < xlen; i++) {
		if (!R->_zerop(x[i]))
			return i;
	}
	return -1;
}}

static const _cl_UP modint_monomial (cl_heap_univpoly_ring* UPR, const cl_ring_element& x, uintL e)
{
	if (!(UPR->basering() == x.ring())) throw runtime_exception();
 {	DeclarePoly(_cl_MI,x);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	if (R->_zerop(x))
		return _cl_UP(UPR, cl_null_GV_I);
	else {
		var sintL len = e+1;
		var cl_GV_MI result = cl_GV_MI(len,R);
		result[e] = x;
		return _cl_UP(UPR, result);
	}
}}

static const cl_ring_element modint_coeff (cl_heap_univpoly_ring* UPR, const _cl_UP& x, uintL index)
{{
	DeclarePoly(cl_GV_MI,x);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	if (index < x.size())
		return cl_MI(R, x[index]);
	else
		return R->zero();
}}

static const _cl_UP modint_create (cl_heap_univpoly_ring* UPR, sintL deg)
{
	if (deg < 0)
		return _cl_UP(UPR, cl_null_GV_I);
	else {
		var sintL len = deg+1;
		var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
		return _cl_UP(UPR, cl_GV_MI(len,R));
	}
}

static void modint_set_coeff (cl_heap_univpoly_ring* UPR, _cl_UP& x, uintL index, const cl_ring_element& y)
{{
	DeclareMutablePoly(cl_GV_MI,x);
	if (!(UPR->basering() == y.ring())) throw runtime_exception();
  {	DeclarePoly(_cl_MI,y);
	if (!(index < x.size())) throw runtime_exception();
	x[index] = y;
}}}

static void modint_finalize (cl_heap_univpoly_ring* UPR, _cl_UP& x)
{{
	DeclareMutablePoly(cl_GV_MI,x); // NB: x is modified by reference!
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var uintL len = x.size();
	if (len > 0)
		modint_normalize(R,x,len);
}}

static const cl_ring_element modint_eval (cl_heap_univpoly_ring* UPR, const _cl_UP& x, const cl_ring_element& y)
{{
	// Method:
	// If x = 0, return 0.
	// If y = 0, return x[0].
	// Else compute (...(x[len-1]*y+x[len-2])*y ...)*y + x[0].
	DeclarePoly(cl_GV_MI,x);
	if (!(UPR->basering() == y.ring())) throw runtime_exception();
  {	DeclarePoly(_cl_MI,y);
	var cl_heap_modint_ring* R = TheModintRing(UPR->basering());
	var uintL len = x.size();
	if (len==0)
		return R->zero();
	if (R->_zerop(y))
		return cl_MI(R, x[0]);
	var sintL i = len-1;
	var _cl_MI z = x[i];
	for ( ; --i >= 0; )
		z = R->_plus(R->_mul(z,y),x[i]);
	return cl_MI(R, z);
}}}

static cl_univpoly_setops modint_setops = {
	modint_fprint,
	modint_equal
};

static cl_univpoly_addops modint_addops = {
	modint_zero,
	modint_zerop,
	modint_plus,
	modint_minus,
	modint_uminus
};

static cl_univpoly_mulops modint_mulops = {
	modint_one,
	modint_canonhom,
	modint_mul,
	modint_square,
	modint_exptpos
};

static cl_univpoly_modulops modint_modulops = {
	modint_scalmul
};

static cl_univpoly_polyops modint_polyops = {
	modint_degree,
	modint_ldegree,
	modint_monomial,
	modint_coeff,
	modint_create,
	modint_set_coeff,
	modint_finalize,
	modint_eval
};

class cl_heap_modint_univpoly_ring : public cl_heap_univpoly_ring {
	SUBCLASS_cl_heap_univpoly_ring()
public:
	// Constructor.
	cl_heap_modint_univpoly_ring (const cl_ring& r);
	// Destructor.
	~cl_heap_modint_univpoly_ring () {}
};

static void cl_heap_modint_univpoly_ring_destructor (cl_heap* pointer)
{
	(*(cl_heap_modint_univpoly_ring*)pointer).~cl_heap_modint_univpoly_ring();
}

cl_class cl_class_modint_univpoly_ring = {
	cl_heap_modint_univpoly_ring_destructor,
	cl_class_flags_univpoly_ring
};

// Constructor.
inline cl_heap_modint_univpoly_ring::cl_heap_modint_univpoly_ring (const cl_ring& r)
	: cl_heap_univpoly_ring (r, &modint_setops, &modint_addops, &modint_mulops, &modint_modulops, &modint_polyops)
{
	type = &cl_class_modint_univpoly_ring;
}

}  // namespace cln
