// Public double float operations.

#ifndef _CL_DFLOAT_H
#define _CL_DFLOAT_H

#include "cln/number.h"
#include "cln/dfloat_class.h"
#include "cln/integer_class.h"
#include "cln/float.h"

namespace cln {

CL_DEFINE_AS_CONVERSION(cl_DF)


// Liefert zu einem Double-Float x : (- x), ein DF.
extern const cl_DF operator- (const cl_DF& x);

// compare(x,y) vergleicht zwei Double-Floats x und y.
// Ergebnis: 0 falls x=y, +1 falls x>y, -1 falls x<y.
extern cl_signean compare (const cl_DF& x, const cl_DF& y);

// equal_hashcode(x) liefert einen equal-invarianten Hashcode für x.
extern uint32 equal_hashcode (const cl_DF& x);

inline bool operator== (const cl_DF& x, const cl_DF& y)
	{ return compare(x,y)==0; }
inline bool operator!= (const cl_DF& x, const cl_DF& y)
	{ return compare(x,y)!=0; }
inline bool operator<= (const cl_DF& x, const cl_DF& y)
	{ return compare(x,y)<=0; }
inline bool operator< (const cl_DF& x, const cl_DF& y)
	{ return compare(x,y)<0; }
inline bool operator>= (const cl_DF& x, const cl_DF& y)
	{ return compare(x,y)>=0; }
inline bool operator> (const cl_DF& x, const cl_DF& y)
	{ return compare(x,y)>0; }

// minusp(x) == (< x 0)
extern bool minusp (const cl_DF& x);

// zerop(x) stellt fest, ob ein Double-Float x = 0.0 ist.
extern bool zerop (const cl_DF& x);

// plusp(x) == (> x 0)
extern bool plusp (const cl_DF& x);

// Liefert zu zwei Double-Float x und y : (+ x y), ein DF.
extern const cl_DF operator+ (const cl_DF& x, const cl_DF& y);
// The C++ compiler may hesitate to do these conversions of its own:
inline const cl_DF operator+ (const cl_DF& x, const double y)
	{ return x + cl_DF(y); }
inline const cl_DF operator+ (const double x, const cl_DF& y)
	{ return cl_DF(x) + y; }

// Liefert zu zwei Double-Float x und y : (- x y), ein DF.
extern const cl_DF operator- (const cl_DF& x, const cl_DF& y);
// The C++ compiler may hesitate to do these conversions of its own:
inline const cl_DF operator- (const cl_DF& x, const double y)
	{ return x - cl_DF(y); }
inline const cl_DF operator- (const double x, const cl_DF& y)
	{ return cl_DF(x) - y; }

// Liefert zu zwei Double-Float x und y : (* x y), ein DF.
extern const cl_DF operator* (const cl_DF& x, const cl_DF& y);
// The C++ compiler may hesitate to do these conversions of its own:
inline const cl_DF operator* (const cl_DF& x, const double y)
	{ return x * cl_DF(y); }
inline const cl_DF operator* (const double x, const cl_DF& y)
	{ return cl_DF(x) * y; }

// Liefert zu einem Double-Float x : (* x x), ein DF.
inline const cl_DF square (const cl_DF& x) { return x*x; }

// Liefert zu zwei Double-Float x und y : (/ x y), ein DF.
extern const cl_DF operator/ (const cl_DF& x, const cl_DF& y);
// The C++ compiler may hesitate to do these conversions of its own:
inline const cl_DF operator/ (const cl_DF& x, const double y)
	{ return x / cl_DF(y); }
inline const cl_DF operator/ (const double x, const cl_DF& y)
	{ return cl_DF(x) / y; }

// Liefert zu einem Double-Float x>=0 : (sqrt x), ein DF.
extern const cl_DF sqrt (const cl_DF& x);

// recip(x) liefert (/ x), wo x ein Double-Float ist.
extern const cl_DF recip (const cl_DF& x);

// abs(x) liefert (abs x), wo x ein Double-Float ist.
extern const cl_DF abs (const cl_DF& x);


// (1+ x), wo x ein Double-Float ist.
inline const cl_DF plus1 (const cl_DF& x)
{
	extern const cl_DF cl_I_to_DF (const cl_I&);
	return x + cl_I_to_DF(cl_I(1));
}

// (1- x), wo x ein Double-Float ist.
inline const cl_DF minus1 (const cl_DF& x)
{
	extern const cl_DF cl_I_to_DF (const cl_I&);
	return x + cl_I_to_DF(cl_I(-1));
}


// ffloor(x) liefert (ffloor x), wo x ein DF ist.
extern const cl_DF ffloor (const cl_DF& x);

// fceiling(x) liefert (fceiling x), wo x ein DF ist.
extern const cl_DF fceiling (const cl_DF& x);

// ftruncate(x) liefert (ftruncate x), wo x ein DF ist.
extern const cl_DF ftruncate (const cl_DF& x);

// fround(x) liefert (fround x), wo x ein DF ist.
extern const cl_DF fround (const cl_DF& x);


// Return type for frounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_DF_fdiv_t {
	cl_DF quotient;
	cl_DF remainder;
// Constructor.
	cl_DF_fdiv_t () {}
	cl_DF_fdiv_t (const cl_DF& q, const cl_DF& r) : quotient(q), remainder(r) {}
};

// ffloor2(x) liefert (ffloor x), wo x ein DF ist.
inline const cl_DF_fdiv_t ffloor2 (const cl_DF& x)
	{ cl_DF q = ffloor(x); return cl_DF_fdiv_t(q,x-q); }

// fceiling2(x) liefert (fceiling x), wo x ein DF ist.
inline const cl_DF_fdiv_t fceiling2 (const cl_DF& x)
	{ cl_DF q = fceiling(x); return cl_DF_fdiv_t(q,x-q); }

// ftruncate2(x) liefert (ftruncate x), wo x ein DF ist.
inline const cl_DF_fdiv_t ftruncate2 (const cl_DF& x)
	{ cl_DF q = ftruncate(x); return cl_DF_fdiv_t(q,x-q); }

// fround2(x) liefert (fround x), wo x ein DF ist.
inline const cl_DF_fdiv_t fround2 (const cl_DF& x)
	{ cl_DF q = fround(x); return cl_DF_fdiv_t(q,x-q); }


// Return type for rounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_DF_div_t {
	cl_I quotient;
	cl_DF remainder;
// Constructor.
	cl_DF_div_t () {}
	cl_DF_div_t (const cl_I& q, const cl_DF& r) : quotient(q), remainder(r) {}
};

// floor2(x) liefert (floor x), wo x ein DF ist.
inline const cl_DF_div_t floor2 (const cl_DF& x)
{
	extern const cl_I cl_DF_to_I (const cl_DF& x);
	cl_DF q = ffloor(x);
	return cl_DF_div_t(cl_DF_to_I(q),x-q);
}
inline const cl_I floor1 (const cl_DF& x)
{
	extern const cl_I cl_DF_to_I (const cl_DF& x);
	return cl_DF_to_I(ffloor(x));
}

// ceiling2(x) liefert (ceiling x), wo x ein DF ist.
inline const cl_DF_div_t ceiling2 (const cl_DF& x)
{
	extern const cl_I cl_DF_to_I (const cl_DF& x);
	cl_DF q = fceiling(x);
	return cl_DF_div_t(cl_DF_to_I(q),x-q);
}
inline const cl_I ceiling1 (const cl_DF& x)
{
	extern const cl_I cl_DF_to_I (const cl_DF& x);
	return cl_DF_to_I(fceiling(x));
}

// truncate2(x) liefert (truncate x), wo x ein DF ist.
inline const cl_DF_div_t truncate2 (const cl_DF& x)
{
	extern const cl_I cl_DF_to_I (const cl_DF& x);
	cl_DF q = ftruncate(x);
	return cl_DF_div_t(cl_DF_to_I(q),x-q);
}
inline const cl_I truncate1 (const cl_DF& x)
{
	extern const cl_I cl_DF_to_I (const cl_DF& x);
	return cl_DF_to_I(ftruncate(x));
}

// round2(x) liefert (round x), wo x ein DF ist.
inline const cl_DF_div_t round2 (const cl_DF& x)
{
	extern const cl_I cl_DF_to_I (const cl_DF& x);
	cl_DF q = fround(x);
	return cl_DF_div_t(cl_DF_to_I(q),x-q);
}
inline const cl_I round1 (const cl_DF& x)
{
	extern const cl_I cl_DF_to_I (const cl_DF& x);
	return cl_DF_to_I(fround(x));
}

// floor2(x,y) liefert (floor x y).
extern const cl_DF_div_t floor2 (const cl_DF& x, const cl_DF& y);
inline const cl_I floor1 (const cl_DF& x, const cl_DF& y) { return floor1(x/y); }

// ceiling2(x,y) liefert (ceiling x y).
extern const cl_DF_div_t ceiling2 (const cl_DF& x, const cl_DF& y);
inline const cl_I ceiling1 (const cl_DF& x, const cl_DF& y) { return ceiling1(x/y); }

// truncate2(x,y) liefert (truncate x y).
extern const cl_DF_div_t truncate2 (const cl_DF& x, const cl_DF& y);
inline const cl_I truncate1 (const cl_DF& x, const cl_DF& y) { return truncate1(x/y); }

// round2(x,y) liefert (round x y).
extern const cl_DF_div_t round2 (const cl_DF& x, const cl_DF& y);
inline const cl_I round1 (const cl_DF& x, const cl_DF& y) { return round1(x/y); }


// Return type for decode_float:
struct decoded_dfloat {
	cl_DF mantissa;
	cl_I exponent;
	cl_DF sign;
// Constructor.
	decoded_dfloat () {}
	decoded_dfloat (const cl_DF& m, const cl_I& e, const cl_DF& s) : mantissa(m), exponent(e), sign(s) {}
};

// decode_float(x) liefert zu einem Float x: (decode-float x).
// x = 0.0 liefert (0.0, 0, 1.0).
// x = (-1)^s * 2^e * m liefert ((-1)^0 * 2^0 * m, e als Integer, (-1)^s).
extern const decoded_dfloat decode_float (const cl_DF& x);

// float_exponent(x) liefert zu einem Float x:
// den Exponenten von (decode-float x).
// x = 0.0 liefert 0.
// x = (-1)^s * 2^e * m liefert e.
extern sintE float_exponent (const cl_DF& x);

// float_radix(x) liefert (float-radix x), wo x ein Float ist.
inline sintL float_radix (const cl_DF& x)
{
	(void)x; // unused x
	return 2;
}

// float_sign(x) liefert (float-sign x), wo x ein Float ist.
extern const cl_DF float_sign (const cl_DF& x);

// float_digits(x) liefert (float-digits x), wo x ein Float ist.
// < ergebnis: ein uintC >0
extern uintC float_digits (const cl_DF& x);

// float_precision(x) liefert (float-precision x), wo x ein Float ist.
// < ergebnis: ein uintC >=0
extern uintC float_precision (const cl_DF& x);


// integer_decode_float(x) liefert zu einem Float x: (integer-decode-float x).
// x = 0.0 liefert (0, 0, 1).
// x = (-1)^s * 2^e * m bei Float-Precision p liefert
//   (Mantisse 2^p * m als Integer, e-p als Integer, (-1)^s als Fixnum).
extern const cl_idecoded_float integer_decode_float (const cl_DF& x);


// scale_float(x,delta) liefert x*2^delta, wo x ein DF ist.
extern const cl_DF scale_float (const cl_DF& x, sintC delta);
extern const cl_DF scale_float (const cl_DF& x, const cl_I& delta);


// max(x,y) liefert (max x y), wo x und y Floats sind.
extern const cl_DF max (const cl_DF& x, const cl_DF& y);

// min(x,y) liefert (min x y), wo x und y Floats sind.
extern const cl_DF min (const cl_DF& x, const cl_DF& y);

// signum(x) liefert (signum x), wo x ein Float ist.
extern const cl_DF signum (const cl_DF& x);


// Konversion zu einem C "float".
extern float float_approx (const cl_DF& x);

// Konversion zu einem C "double".
extern double double_approx (const cl_DF& x);


// This could be optimized to use in-place operations.
inline cl_DF& operator+= (cl_DF& x, const cl_DF& y) { return x = x + y; }
inline cl_DF& operator+= (cl_DF& x, const double y) { return x = x + y; }
inline cl_DF& operator++ /* prefix */ (cl_DF& x) { return x = plus1(x); }
inline void operator++ /* postfix */ (cl_DF& x, int dummy) { (void)dummy; x = plus1(x); }
inline cl_DF& operator-= (cl_DF& x, const cl_DF& y) { return x = x - y; }
inline cl_DF& operator-= (cl_DF& x, const double y) { return x = x - y; }
inline cl_DF& operator-- /* prefix */ (cl_DF& x) { return x = minus1(x); }
inline void operator-- /* postfix */ (cl_DF& x, int dummy) { (void)dummy; x = minus1(x); }
inline cl_DF& operator*= (cl_DF& x, const cl_DF& y) { return x = x * y; }
inline cl_DF& operator*= (cl_DF& x, const double y) { return x = x * y; }
inline cl_DF& operator/= (cl_DF& x, const cl_DF& y) { return x = x / y; }
inline cl_DF& operator/= (cl_DF& x, const double y) { return x = x / y; }


/* */


// Runtime typing support.
extern cl_class cl_class_dfloat;


// Debugging support.
#ifdef CL_DEBUG
extern int cl_DF_debug_module;
CL_FORCE_LINK(cl_DF_debug_dummy, cl_DF_debug_module)
#endif

}  // namespace cln

#endif /* _CL_DFLOAT_H */
