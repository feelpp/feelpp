// Public short float operations.

#ifndef _CL_SFLOAT_H
#define _CL_SFLOAT_H

#include "cln/number.h"
#include "cln/sfloat_class.h"
#include "cln/integer_class.h"
#include "cln/float.h"

namespace cln {

CL_DEFINE_AS_CONVERSION(cl_SF)


// Liefert zu einem Short-Float x : (- x), ein SF.
extern const cl_SF operator- (const cl_SF& x);

// compare(x,y) vergleicht zwei Short-Floats x und y.
// Ergebnis: 0 falls x=y, +1 falls x>y, -1 falls x<y.
extern cl_signean compare (const cl_SF& x, const cl_SF& y);

// equal_hashcode(x) liefert einen equal-invarianten Hashcode für x.
extern uint32 equal_hashcode (const cl_SF& x);

inline bool operator== (const cl_SF& x, const cl_SF& y)
	{ return compare(x,y)==0; }
inline bool operator!= (const cl_SF& x, const cl_SF& y)
	{ return compare(x,y)!=0; }
inline bool operator<= (const cl_SF& x, const cl_SF& y)
	{ return compare(x,y)<=0; }
inline bool operator< (const cl_SF& x, const cl_SF& y)
	{ return compare(x,y)<0; }
inline bool operator>= (const cl_SF& x, const cl_SF& y)
	{ return compare(x,y)>=0; }
inline bool operator> (const cl_SF& x, const cl_SF& y)
	{ return compare(x,y)>0; }

// minusp(x) == (< x 0)
extern bool minusp (const cl_SF& x);

// zerop(x) stellt fest, ob ein Short-Float x = 0.0 ist.
extern bool zerop (const cl_SF& x);

// plusp(x) == (> x 0)
extern bool plusp (const cl_SF& x);

// Liefert zu zwei Short-Float x und y : (+ x y), ein SF.
extern const cl_SF operator+ (const cl_SF& x, const cl_SF& y);

// Liefert zu zwei Short-Float x und y : (- x y), ein SF.
extern const cl_SF operator- (const cl_SF& x, const cl_SF& y);

// Liefert zu zwei Short-Float x und y : (* x y), ein SF.
extern const cl_SF operator* (const cl_SF& x, const cl_SF& y);

// Liefert zu einem Short-Float x : (* x x), ein SF.
inline const cl_SF square (const cl_SF& x) { return x*x; }

// Liefert zu zwei Short-Float x und y : (/ x y), ein SF.
extern const cl_SF operator/ (const cl_SF& x, const cl_SF& y);

// Liefert zu einem Short-Float x>=0 : (sqrt x), ein SF.
extern const cl_SF sqrt (const cl_SF& x);

// recip(x) liefert (/ x), wo x ein Short-Float ist.
extern const cl_SF recip (const cl_SF& x);

// abs(x) liefert (abs x), wo x ein Short-Float ist.
extern const cl_SF abs (const cl_SF& x);


// (1+ x), wo x ein Short-Float ist.
inline const cl_SF plus1 (const cl_SF& x)
{
	extern const cl_SF cl_I_to_SF (const cl_I&);
	return x + cl_I_to_SF(cl_I(1));
}

// (1- x), wo x ein Short-Float ist.
inline const cl_SF minus1 (const cl_SF& x)
{
	extern const cl_SF cl_I_to_SF (const cl_I&);
	return x + cl_I_to_SF(cl_I(-1));
}


// ffloor(x) liefert (ffloor x), wo x ein SF ist.
extern const cl_SF ffloor (const cl_SF& x);

// fceiling(x) liefert (fceiling x), wo x ein SF ist.
extern const cl_SF fceiling (const cl_SF& x);

// ftruncate(x) liefert (ftruncate x), wo x ein SF ist.
extern const cl_SF ftruncate (const cl_SF& x);

// fround(x) liefert (fround x), wo x ein SF ist.
extern const cl_SF fround (const cl_SF& x);


// Return type for frounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_SF_fdiv_t {
	cl_SF quotient;
	cl_SF remainder;
// Constructor.
	cl_SF_fdiv_t () {}
	cl_SF_fdiv_t (const cl_SF& q, const cl_SF& r) : quotient(q), remainder(r) {}
};

// ffloor2(x) liefert (ffloor x), wo x ein SF ist.
inline const cl_SF_fdiv_t ffloor2 (const cl_SF& x)
	{ cl_SF q = ffloor(x); return cl_SF_fdiv_t(q,x-q); }

// fceiling2(x) liefert (fceiling x), wo x ein SF ist.
inline const cl_SF_fdiv_t fceiling2 (const cl_SF& x)
	{ cl_SF q = fceiling(x); return cl_SF_fdiv_t(q,x-q); }

// ftruncate2(x) liefert (ftruncate x), wo x ein SF ist.
inline const cl_SF_fdiv_t ftruncate2 (const cl_SF& x)
	{ cl_SF q = ftruncate(x); return cl_SF_fdiv_t(q,x-q); }

// fround2(x) liefert (fround x), wo x ein SF ist.
inline const cl_SF_fdiv_t fround2 (const cl_SF& x)
	{ cl_SF q = fround(x); return cl_SF_fdiv_t(q,x-q); }


// Return type for rounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_SF_div_t {
	cl_I quotient;
	cl_SF remainder;
// Constructor.
	cl_SF_div_t () {}
	cl_SF_div_t (const cl_I& q, const cl_SF& r) : quotient(q), remainder(r) {}
};

// floor2(x) liefert (floor x), wo x ein SF ist.
inline const cl_SF_div_t floor2 (const cl_SF& x)
{
	extern const cl_I cl_SF_to_I (const cl_SF& x);
	cl_SF q = ffloor(x);
	return cl_SF_div_t(cl_SF_to_I(q),x-q);
}
inline const cl_I floor1 (const cl_SF& x)
{
	extern const cl_I cl_SF_to_I (const cl_SF& x);
	return cl_SF_to_I(ffloor(x));
}

// ceiling2(x) liefert (ceiling x), wo x ein SF ist.
inline const cl_SF_div_t ceiling2 (const cl_SF& x)
{
	extern const cl_I cl_SF_to_I (const cl_SF& x);
	cl_SF q = fceiling(x);
	return cl_SF_div_t(cl_SF_to_I(q),x-q);
}
inline const cl_I ceiling1 (const cl_SF& x)
{
	extern const cl_I cl_SF_to_I (const cl_SF& x);
	return cl_SF_to_I(fceiling(x));
}

// truncate2(x) liefert (truncate x), wo x ein SF ist.
inline const cl_SF_div_t truncate2 (const cl_SF& x)
{
	extern const cl_I cl_SF_to_I (const cl_SF& x);
	cl_SF q = ftruncate(x);
	return cl_SF_div_t(cl_SF_to_I(q),x-q);
}
inline const cl_I truncate1 (const cl_SF& x)
{
	extern const cl_I cl_SF_to_I (const cl_SF& x);
	return cl_SF_to_I(ftruncate(x));
}

// round2(x) liefert (round x), wo x ein SF ist.
inline const cl_SF_div_t round2 (const cl_SF& x)
{
	extern const cl_I cl_SF_to_I (const cl_SF& x);
	cl_SF q = fround(x);
	return cl_SF_div_t(cl_SF_to_I(q),x-q);
}
inline const cl_I round1 (const cl_SF& x)
{
	extern const cl_I cl_SF_to_I (const cl_SF& x);
	return cl_SF_to_I(fround(x));
}

// floor2(x,y) liefert (floor x y).
extern const cl_SF_div_t floor2 (const cl_SF& x, const cl_SF& y);
inline const cl_I floor1 (const cl_SF& x, const cl_SF& y) { return floor1(x/y); }

// ceiling2(x,y) liefert (ceiling x y).
extern const cl_SF_div_t ceiling2 (const cl_SF& x, const cl_SF& y);
inline const cl_I ceiling1 (const cl_SF& x, const cl_SF& y) { return ceiling1(x/y); }

// truncate2(x,y) liefert (truncate x y).
extern const cl_SF_div_t truncate2 (const cl_SF& x, const cl_SF& y);
inline const cl_I truncate1 (const cl_SF& x, const cl_SF& y) { return truncate1(x/y); }

// round2(x,y) liefert (round x y).
extern const cl_SF_div_t round2 (const cl_SF& x, const cl_SF& y);
inline const cl_I round1 (const cl_SF& x, const cl_SF& y) { return round1(x/y); }


// Return type for decode_float:
struct decoded_sfloat {
	cl_SF mantissa;
	cl_I exponent;
	cl_SF sign;
// Constructor.
	decoded_sfloat () {}
	decoded_sfloat (const cl_SF& m, const cl_I& e, const cl_SF& s) : mantissa(m), exponent(e), sign(s) {}
};

// decode_float(x) liefert zu einem Float x: (decode-float x).
// x = 0.0 liefert (0.0, 0, 1.0).
// x = (-1)^s * 2^e * m liefert ((-1)^0 * 2^0 * m, e als Integer, (-1)^s).
extern const decoded_sfloat decode_float (const cl_SF& x);

// float_exponent(x) liefert zu einem Float x:
// den Exponenten von (decode-float x).
// x = 0.0 liefert 0.
// x = (-1)^s * 2^e * m liefert e.
extern sintE float_exponent (const cl_SF& x);

// float_radix(x) liefert (float-radix x), wo x ein Float ist.
inline sintL float_radix (const cl_SF& x)
{
	(void)x; // unused x
	return 2;
}

// float_sign(x) liefert (float-sign x), wo x ein Float ist.
extern const cl_SF float_sign (const cl_SF& x);

// float_digits(x) liefert (float-digits x), wo x ein Float ist.
// < ergebnis: ein uintC >0
extern uintC float_digits (const cl_SF& x);

// float_precision(x) liefert (float-precision x), wo x ein Float ist.
// < ergebnis: ein uintC >=0
extern uintC float_precision (const cl_SF& x);


// integer_decode_float(x) liefert zu einem Float x: (integer-decode-float x).
// x = 0.0 liefert (0, 0, 1).
// x = (-1)^s * 2^e * m bei Float-Precision p liefert
//   (Mantisse 2^p * m als Integer, e-p als Integer, (-1)^s als Fixnum).
extern const cl_idecoded_float integer_decode_float (const cl_SF& x);


// scale_float(x,delta) liefert x*2^delta, wo x ein SF ist.
extern const cl_SF scale_float (const cl_SF& x, sintC delta);
extern const cl_SF scale_float (const cl_SF& x, const cl_I& delta);


// max(x,y) liefert (max x y), wo x und y Floats sind.
extern const cl_SF max (const cl_SF& x, const cl_SF& y);

// min(x,y) liefert (min x y), wo x und y Floats sind.
extern const cl_SF min (const cl_SF& x, const cl_SF& y);

// signum(x) liefert (signum x), wo x ein Float ist.
extern const cl_SF signum (const cl_SF& x);


// Konversion zu einem C "float".
extern float float_approx (const cl_SF& x);

// Konversion zu einem C "double".
extern double double_approx (const cl_SF& x);


// This could be optimized to use in-place operations.
inline cl_SF& operator+= (cl_SF& x, const cl_SF& y) { return x = x + y; }
inline cl_SF& operator++ /* prefix */ (cl_SF& x) { return x = plus1(x); }
inline void operator++ /* postfix */ (cl_SF& x, int dummy) { (void)dummy; x = plus1(x); }
inline cl_SF& operator-= (cl_SF& x, const cl_SF& y) { return x = x - y; }
inline cl_SF& operator-- /* prefix */ (cl_SF& x) { return x = minus1(x); }
inline void operator-- /* postfix */ (cl_SF& x, int dummy) { (void)dummy; x = minus1(x); }
inline cl_SF& operator*= (cl_SF& x, const cl_SF& y) { return x = x * y; }
inline cl_SF& operator/= (cl_SF& x, const cl_SF& y) { return x = x / y; }


// Runtime typing support.
extern cl_class cl_class_sfloat;
CL_FORCE_LINK(cl_SF_classes_dummy, cl_class_sfloat)


// Debugging support.
#ifdef CL_DEBUG
extern int cl_SF_debug_module;
CL_FORCE_LINK(cl_SF_debug_dummy, cl_SF_debug_module)
#endif

}  // namespace cln

#endif /* _CL_SFLOAT_H */
