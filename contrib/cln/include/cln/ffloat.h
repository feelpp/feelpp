// Public single float operations.

#ifndef _CL_FFLOAT_H
#define _CL_FFLOAT_H

#include "cln/number.h"
#include "cln/ffloat_class.h"
#include "cln/integer_class.h"
#include "cln/float.h"

namespace cln {

CL_DEFINE_AS_CONVERSION(cl_FF)


// Liefert zu einem Single-Float x : (- x), ein FF.
extern const cl_FF operator- (const cl_FF& x);

// compare(x,y) vergleicht zwei Single-Floats x und y.
// Ergebnis: 0 falls x=y, +1 falls x>y, -1 falls x<y.
extern cl_signean compare (const cl_FF& x, const cl_FF& y);

// equal_hashcode(x) liefert einen equal-invarianten Hashcode für x.
extern uint32 equal_hashcode (const cl_FF& x);

inline bool operator== (const cl_FF& x, const cl_FF& y)
	{ return compare(x,y)==0; }
inline bool operator!= (const cl_FF& x, const cl_FF& y)
	{ return compare(x,y)!=0; }
inline bool operator<= (const cl_FF& x, const cl_FF& y)
	{ return compare(x,y)<=0; }
inline bool operator< (const cl_FF& x, const cl_FF& y)
	{ return compare(x,y)<0; }
inline bool operator>= (const cl_FF& x, const cl_FF& y)
	{ return compare(x,y)>=0; }
inline bool operator> (const cl_FF& x, const cl_FF& y)
	{ return compare(x,y)>0; }

// minusp(x) == (< x 0)
extern bool minusp (const cl_FF& x);

// zerop(x) stellt fest, ob ein Single-Float x = 0.0 ist.
extern bool zerop (const cl_FF& x);

// plusp(x) == (> x 0)
extern bool plusp (const cl_FF& x);

// Liefert zu zwei Single-Float x und y : (+ x y), ein FF.
extern const cl_FF operator+ (const cl_FF& x, const cl_FF& y);
// The C++ compiler may hesitate to do these conversions of its own:
inline const cl_FF operator+ (const cl_FF& x, const float y)
	{ return x + cl_FF(y); }
inline const cl_FF operator+ (const float x, const cl_FF& y)
	{ return cl_FF(x) + y; }

// Liefert zu zwei Single-Float x und y : (- x y), ein FF.
extern const cl_FF operator- (const cl_FF& x, const cl_FF& y);
// The C++ compiler may hesitate to do these conversions of its own:
inline const cl_FF operator- (const cl_FF& x, const float y)
	{ return x - cl_FF(y); }
inline const cl_FF operator- (const float x, const cl_FF& y)
	{ return cl_FF(x) - y; }

// Liefert zu zwei Single-Float x und y : (* x y), ein FF.
extern const cl_FF operator* (const cl_FF& x, const cl_FF& y);
// The C++ compiler may hesitate to do these conversions of its own:
inline const cl_FF operator* (const cl_FF& x, const float y)
	{ return x * cl_FF(y); }
inline const cl_FF operator* (const float x, const cl_FF& y)
	{ return cl_FF(x) * y; }

// Liefert zu einem Single-Float x : (* x x), ein FF.
inline const cl_FF square (const cl_FF& x) { return x*x; }

// Liefert zu zwei Single-Float x und y : (/ x y), ein FF.
extern const cl_FF operator/ (const cl_FF& x, const cl_FF& y);
// The C++ compiler may hesitate to do these conversions of its own:
inline const cl_FF operator/ (const cl_FF& x, const float y)
	{ return x / cl_FF(y); }
inline const cl_FF operator/ (const float x, const cl_FF& y)
	{ return cl_FF(x) / y; }

// Liefert zu einem Single-Float x>=0 : (sqrt x), ein FF.
extern const cl_FF sqrt (const cl_FF& x);

// recip(x) liefert (/ x), wo x ein Single-Float ist.
extern const cl_FF recip (const cl_FF& x);

// abs(x) liefert (abs x), wo x ein Single-Float ist.
extern const cl_FF abs (const cl_FF& x);


// (1+ x), wo x ein Single-Float ist.
inline const cl_FF plus1 (const cl_FF& x)
{
	extern const cl_FF cl_I_to_FF (const cl_I&);
	return x + cl_I_to_FF(cl_I(1));
}

// (1- x), wo x ein Single-Float ist.
inline const cl_FF minus1 (const cl_FF& x)
{
	extern const cl_FF cl_I_to_FF (const cl_I&);
	return x + cl_I_to_FF(cl_I(-1));
}


// ffloor(x) liefert (ffloor x), wo x ein FF ist.
extern const cl_FF ffloor (const cl_FF& x);

// fceiling(x) liefert (fceiling x), wo x ein FF ist.
extern const cl_FF fceiling (const cl_FF& x);

// ftruncate(x) liefert (ftruncate x), wo x ein FF ist.
extern const cl_FF ftruncate (const cl_FF& x);

// fround(x) liefert (fround x), wo x ein FF ist.
extern const cl_FF fround (const cl_FF& x);


// Return type for frounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_FF_fdiv_t {
	cl_FF quotient;
	cl_FF remainder;
// Constructor.
	cl_FF_fdiv_t () {}
	cl_FF_fdiv_t (const cl_FF& q, const cl_FF& r) : quotient(q), remainder(r) {}
};

// ffloor2(x) liefert (ffloor x), wo x ein FF ist.
inline const cl_FF_fdiv_t ffloor2 (const cl_FF& x)
	{ cl_FF q = ffloor(x); return cl_FF_fdiv_t(q,x-q); }

// fceiling2(x) liefert (fceiling x), wo x ein FF ist.
inline const cl_FF_fdiv_t fceiling2 (const cl_FF& x)
	{ cl_FF q = fceiling(x); return cl_FF_fdiv_t(q,x-q); }

// ftruncate2(x) liefert (ftruncate x), wo x ein FF ist.
inline const cl_FF_fdiv_t ftruncate2 (const cl_FF& x)
	{ cl_FF q = ftruncate(x); return cl_FF_fdiv_t(q,x-q); }

// fround2(x) liefert (fround x), wo x ein FF ist.
inline const cl_FF_fdiv_t fround2 (const cl_FF& x)
	{ cl_FF q = fround(x); return cl_FF_fdiv_t(q,x-q); }


// Return type for rounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_FF_div_t {
	cl_I quotient;
	cl_FF remainder;
// Constructor.
	cl_FF_div_t () {}
	cl_FF_div_t (const cl_I& q, const cl_FF& r) : quotient(q), remainder(r) {}
};

// floor2(x) liefert (floor x), wo x ein FF ist.
inline const cl_FF_div_t floor2 (const cl_FF& x)
{
	extern const cl_I cl_FF_to_I (const cl_FF& x);
	cl_FF q = ffloor(x);
	return cl_FF_div_t(cl_FF_to_I(q),x-q);
}
inline const cl_I floor1 (const cl_FF& x)
{
	extern const cl_I cl_FF_to_I (const cl_FF& x);
	return cl_FF_to_I(ffloor(x));
}

// ceiling2(x) liefert (ceiling x), wo x ein FF ist.
inline const cl_FF_div_t ceiling2 (const cl_FF& x)
{
	extern const cl_I cl_FF_to_I (const cl_FF& x);
	cl_FF q = fceiling(x);
	return cl_FF_div_t(cl_FF_to_I(q),x-q);
}
inline const cl_I ceiling1 (const cl_FF& x)
{
	extern const cl_I cl_FF_to_I (const cl_FF& x);
	return cl_FF_to_I(fceiling(x));
}

// truncate2(x) liefert (truncate x), wo x ein FF ist.
inline const cl_FF_div_t truncate2 (const cl_FF& x)
{
	extern const cl_I cl_FF_to_I (const cl_FF& x);
	cl_FF q = ftruncate(x);
	return cl_FF_div_t(cl_FF_to_I(q),x-q);
}
inline const cl_I truncate1 (const cl_FF& x)
{
	extern const cl_I cl_FF_to_I (const cl_FF& x);
	return cl_FF_to_I(ftruncate(x));
}

// round2(x) liefert (round x), wo x ein FF ist.
inline const cl_FF_div_t round2 (const cl_FF& x)
{
	extern const cl_I cl_FF_to_I (const cl_FF& x);
	cl_FF q = fround(x);
	return cl_FF_div_t(cl_FF_to_I(q),x-q);
}
inline const cl_I round1 (const cl_FF& x)
{
	extern const cl_I cl_FF_to_I (const cl_FF& x);
	return cl_FF_to_I(fround(x));
}

// floor2(x,y) liefert (floor x y).
extern const cl_FF_div_t floor2 (const cl_FF& x, const cl_FF& y);
inline const cl_I floor1 (const cl_FF& x, const cl_FF& y) { return floor1(x/y); }

// ceiling2(x,y) liefert (ceiling x y).
extern const cl_FF_div_t ceiling2 (const cl_FF& x, const cl_FF& y);
inline const cl_I ceiling1 (const cl_FF& x, const cl_FF& y) { return ceiling1(x/y); }

// truncate2(x,y) liefert (truncate x y).
extern const cl_FF_div_t truncate2 (const cl_FF& x, const cl_FF& y);
inline const cl_I truncate1 (const cl_FF& x, const cl_FF& y) { return truncate1(x/y); }

// round2(x,y) liefert (round x y).
extern const cl_FF_div_t round2 (const cl_FF& x, const cl_FF& y);
inline const cl_I round1 (const cl_FF& x, const cl_FF& y) { return round1(x/y); }


// Return type for decode_float:
struct decoded_ffloat {
	cl_FF mantissa;
	cl_I exponent;
	cl_FF sign;
// Constructor.
	decoded_ffloat () {}
	decoded_ffloat (const cl_FF& m, const cl_I& e, const cl_FF& s) : mantissa(m), exponent(e), sign(s) {}
};

// decode_float(x) liefert zu einem Float x: (decode-float x).
// x = 0.0 liefert (0.0, 0, 1.0).
// x = (-1)^s * 2^e * m liefert ((-1)^0 * 2^0 * m, e als Integer, (-1)^s).
extern const decoded_ffloat decode_float (const cl_FF& x);

// float_exponent(x) liefert zu einem Float x:
// den Exponenten von (decode-float x).
// x = 0.0 liefert 0.
// x = (-1)^s * 2^e * m liefert e.
extern sintE float_exponent (const cl_FF& x);

// float_radix(x) liefert (float-radix x), wo x ein Float ist.
inline sintL float_radix (const cl_FF& x)
{
	(void)x; // unused x
	return 2;
}

// float_sign(x) liefert (float-sign x), wo x ein Float ist.
extern const cl_FF float_sign (const cl_FF& x);

// float_digits(x) liefert (float-digits x), wo x ein Float ist.
// < ergebnis: ein uintC >0
extern uintC float_digits (const cl_FF& x);

// float_precision(x) liefert (float-precision x), wo x ein Float ist.
// < ergebnis: ein uintC >=0
extern uintC float_precision (const cl_FF& x);


// integer_decode_float(x) liefert zu einem Float x: (integer-decode-float x).
// x = 0.0 liefert (0, 0, 1).
// x = (-1)^s * 2^e * m bei Float-Precision p liefert
//   (Mantisse 2^p * m als Integer, e-p als Integer, (-1)^s als Fixnum).
extern const cl_idecoded_float integer_decode_float (const cl_FF& x);


// scale_float(x,delta) liefert x*2^delta, wo x ein FF ist.
extern const cl_FF scale_float (const cl_FF& x, sintC delta);
extern const cl_FF scale_float (const cl_FF& x, const cl_I& delta);


// max(x,y) liefert (max x y), wo x und y Floats sind.
extern const cl_FF max (const cl_FF& x, const cl_FF& y);

// min(x,y) liefert (min x y), wo x und y Floats sind.
extern const cl_FF min (const cl_FF& x, const cl_FF& y);

// signum(x) liefert (signum x), wo x ein Float ist.
extern const cl_FF signum (const cl_FF& x);


// Konversion zu einem C "float".
extern float float_approx (const cl_FF& x);

// Konversion zu einem C "double".
extern double double_approx (const cl_FF& x);


// This could be optimized to use in-place operations.
inline cl_FF& operator+= (cl_FF& x, const cl_FF& y) { return x = x + y; }
inline cl_FF& operator+= (cl_FF& x, const float y) { return x = x + y; }
inline cl_FF& operator++ /* prefix */ (cl_FF& x) { return x = plus1(x); }
inline void operator++ /* postfix */ (cl_FF& x, int dummy) { (void)dummy; x = plus1(x); }
inline cl_FF& operator-= (cl_FF& x, const cl_FF& y) { return x = x - y; }
inline cl_FF& operator-= (cl_FF& x, const float y) { return x = x - y; }
inline cl_FF& operator-- /* prefix */ (cl_FF& x) { return x = minus1(x); }
inline void operator-- /* postfix */ (cl_FF& x, int dummy) { (void)dummy; x = minus1(x); }
inline cl_FF& operator*= (cl_FF& x, const cl_FF& y) { return x = x * y; }
inline cl_FF& operator*= (cl_FF& x, const float y) { return x = x * y; }
inline cl_FF& operator/= (cl_FF& x, const cl_FF& y) { return x = x / y; }
inline cl_FF& operator/= (cl_FF& x, const float y) { return x = x / y; }


/* */


// Runtime typing support.
extern cl_class cl_class_ffloat;
#ifdef CL_WIDE_POINTERS
CL_FORCE_LINK(cl_FF_classes_dummy, cl_class_ffloat)
#endif


// Debugging support.
#ifdef CL_DEBUG
extern int cl_FF_debug_module;
CL_FORCE_LINK(cl_FF_debug_dummy, cl_FF_debug_module)
#endif

}  // namespace cln

#endif /* _CL_FFLOAT_H */
