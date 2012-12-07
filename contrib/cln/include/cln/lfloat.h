// Public long float operations.

#ifndef _CL_LFLOAT_H
#define _CL_LFLOAT_H

#include "cln/number.h"
#include "cln/lfloat_class.h"
#include "cln/integer_class.h"
#include "cln/float.h"

namespace cln {

CL_DEFINE_AS_CONVERSION(cl_LF)


// Liefert zu einem Long-Float x : (- x), ein LF.
extern const cl_LF operator- (const cl_LF& x);

// compare(x,y) vergleicht zwei Long-Floats x und y.
// Ergebnis: 0 falls x=y, +1 falls x>y, -1 falls x<y.
extern cl_signean compare (const cl_LF& x, const cl_LF& y);

// equal_hashcode(x) liefert einen equal-invarianten Hashcode für x.
extern uint32 equal_hashcode (const cl_LF& x);

inline bool operator== (const cl_LF& x, const cl_LF& y)
	{ return compare(x,y)==0; }
inline bool operator!= (const cl_LF& x, const cl_LF& y)
	{ return compare(x,y)!=0; }
inline bool operator<= (const cl_LF& x, const cl_LF& y)
	{ return compare(x,y)<=0; }
inline bool operator< (const cl_LF& x, const cl_LF& y)
	{ return compare(x,y)<0; }
inline bool operator>= (const cl_LF& x, const cl_LF& y)
	{ return compare(x,y)>=0; }
inline bool operator> (const cl_LF& x, const cl_LF& y)
	{ return compare(x,y)>0; }

// minusp(x) == (< x 0)
extern bool minusp (const cl_LF& x);

// zerop(x) stellt fest, ob ein Long-Float x = 0.0 ist.
extern bool zerop (const cl_LF& x);

// plusp(x) == (> x 0)
extern bool plusp (const cl_LF& x);

// Liefert zu zwei Long-Float x und y : (+ x y), ein LF.
extern const cl_LF operator+ (const cl_LF& x, const cl_LF& y);

// Liefert zu zwei Long-Float x und y : (- x y), ein LF.
extern const cl_LF operator- (const cl_LF& x, const cl_LF& y);

// Liefert zu zwei Long-Float x und y : (* x y), ein LF.
extern const cl_LF operator* (const cl_LF& x, const cl_LF& y);
// Spezialfall x oder y Integer oder rationale Zahl.
inline const cl_R operator* (const cl_LF& x, const cl_I& y)
{
	extern const cl_R cl_LF_I_mul (const cl_LF&, const cl_I&);
	return cl_LF_I_mul(x,y);
}
inline const cl_R operator* (const cl_I& x, const cl_LF& y)
{
	extern const cl_R cl_LF_I_mul (const cl_LF&, const cl_I&);
	return cl_LF_I_mul(y,x);
}
inline const cl_R operator* (const cl_LF& x, const cl_RA& y)
{
	extern const cl_R cl_LF_RA_mul (const cl_LF&, const cl_RA&);
	return cl_LF_RA_mul(x,y);
}
inline const cl_R operator* (const cl_RA& x, const cl_LF& y)
{
	extern const cl_R cl_LF_RA_mul (const cl_LF&, const cl_RA&);
	return cl_LF_RA_mul(y,x);
}
// Dem C++-Compiler muß man auch das Folgende sagen (wg. `int * cl_LF' u.ä.):
inline const cl_R operator* (const int x, const cl_LF& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const unsigned int x, const cl_LF& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const long x, const cl_LF& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const unsigned long x, const cl_LF& y)
	{ return cl_I(x) * y; }
#ifdef HAVE_LONGLONG
inline const cl_R operator* (const long long x, const cl_LF& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const unsigned long long x, const cl_LF& y)
	{ return cl_I(x) * y; }
#endif
inline const cl_R operator* (const cl_LF& x, const int y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_LF& x, const unsigned int y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_LF& x, const long y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_LF& x, const unsigned long y)
	{ return x * cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_R operator* (const cl_LF& x, const long long y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_LF& x, const unsigned long long y)
	{ return x * cl_I(y); }
#endif
// Spezialfall x = y.
// Liefert zu einem Long-Float x : (* x x), ein LF.
extern const cl_LF square (const cl_LF& x);

// Liefert zu zwei Long-Float x und y : (/ x y), ein LF.
extern const cl_LF operator/ (const cl_LF& x, const cl_LF& y);
// Spezialfall x oder y Integer oder rationale Zahl.
inline const cl_LF operator/ (const cl_LF& x, const cl_I& y)
{
	extern const cl_LF cl_LF_I_div (const cl_LF& x, const cl_I& y);
	return cl_LF_I_div(x,y);
}
inline const cl_R operator/ (const cl_I& x, const cl_LF& y)
{
	extern const cl_R cl_I_LF_div (const cl_I& x, const cl_LF& y);
	return cl_I_LF_div(x,y);
}
inline const cl_LF operator/ (const cl_LF& x, const cl_RA& y)
{
	extern const cl_LF cl_LF_RA_div (const cl_LF& x, const cl_RA& y);
	return cl_LF_RA_div(x,y);
}
inline const cl_R operator/ (const cl_RA& x, const cl_LF& y)
{
	extern const cl_R cl_RA_LF_div (const cl_RA& x, const cl_LF& y);
	return cl_RA_LF_div(x,y);
}
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_LF operator/ (const cl_LF& x, const int y)
	{ return x / cl_I(y); }
inline const cl_LF operator/ (const cl_LF& x, const unsigned int y)
	{ return x / cl_I(y); }
inline const cl_LF operator/ (const cl_LF& x, const long y)
	{ return x / cl_I(y); }
inline const cl_LF operator/ (const cl_LF& x, const unsigned long y)
	{ return x / cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_LF operator/ (const cl_LF& x, const long long y)
	{ return x / cl_I(y); }
inline const cl_LF operator/ (const cl_LF& x, const unsigned long long y)
	{ return x / cl_I(y); }
#endif
inline const cl_R operator/ (const int x, const cl_LF& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const unsigned int x, const cl_LF& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const long x, const cl_LF& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const unsigned long x, const cl_LF& y)
	{ return cl_I(x) / y; }
#ifdef HAVE_LONGLONG
inline const cl_R operator/ (const long long x, const cl_LF& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const unsigned long long x, const cl_LF& y)
	{ return cl_I(x) / y; }
#endif

// Liefert zu einem Long-Float x>=0 : (sqrt x), ein LF.
extern const cl_LF sqrt (const cl_LF& x);

// recip(x) liefert (/ x), wo x ein Long-Float ist.
extern const cl_LF recip (const cl_LF& x);

// abs(x) liefert (abs x), wo x ein Long-Float ist.
extern const cl_LF abs (const cl_LF& x);


// (1+ x), wo x ein Long-Float ist.
extern const cl_LF plus1 (const cl_LF& x);

// (1- x), wo x ein Long-Float ist.
extern const cl_LF minus1 (const cl_LF& x);


// ffloor(x) liefert (ffloor x), wo x ein LF ist.
extern const cl_LF ffloor (const cl_LF& x);

// fceiling(x) liefert (fceiling x), wo x ein LF ist.
extern const cl_LF fceiling (const cl_LF& x);

// ftruncate(x) liefert (ftruncate x), wo x ein LF ist.
extern const cl_LF ftruncate (const cl_LF& x);

// fround(x) liefert (fround x), wo x ein LF ist.
extern const cl_LF fround (const cl_LF& x);


// Return type for frounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_LF_fdiv_t {
	cl_LF quotient;
	cl_LF remainder;
// Constructor.
	cl_LF_fdiv_t () {}
	cl_LF_fdiv_t (const cl_LF& q, const cl_LF& r) : quotient(q), remainder(r) {}
};

// ffloor2(x) liefert (ffloor x), wo x ein LF ist.
inline const cl_LF_fdiv_t ffloor2 (const cl_LF& x)
{
	extern const cl_LF LF_LF_minus_LF (const cl_LF&, const cl_LF&);
	cl_LF q = ffloor(x);
	return cl_LF_fdiv_t(q,LF_LF_minus_LF(x,q));
}

// fceiling2(x) liefert (fceiling x), wo x ein LF ist.
inline const cl_LF_fdiv_t fceiling2 (const cl_LF& x)
{
	extern const cl_LF LF_LF_minus_LF (const cl_LF&, const cl_LF&);
	cl_LF q = fceiling(x);
	return cl_LF_fdiv_t(q,LF_LF_minus_LF(x,q));
}

// ftruncate2(x) liefert (ftruncate x), wo x ein LF ist.
inline const cl_LF_fdiv_t ftruncate2 (const cl_LF& x)
{
	extern const cl_LF LF_LF_minus_LF (const cl_LF&, const cl_LF&);
	cl_LF q = ftruncate(x);
	return cl_LF_fdiv_t(q,LF_LF_minus_LF(x,q));
}

// fround2(x) liefert (fround x), wo x ein LF ist.
inline const cl_LF_fdiv_t fround2 (const cl_LF& x)
{
	extern const cl_LF LF_LF_minus_LF (const cl_LF&, const cl_LF&);
	cl_LF q = fround(x);
	return cl_LF_fdiv_t(q,LF_LF_minus_LF(x,q));
}


// Return type for rounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_LF_div_t {
	cl_I quotient;
	cl_LF remainder;
// Constructor.
	cl_LF_div_t () {}
	cl_LF_div_t (const cl_I& q, const cl_LF& r) : quotient(q), remainder(r) {}
};

// floor2(x) liefert (floor x), wo x ein LF ist.
inline const cl_LF_div_t floor2 (const cl_LF& x)
{
	extern const cl_LF LF_LF_minus_LF (const cl_LF&, const cl_LF&);
	extern const cl_I cl_LF_to_I (const cl_LF& x);
	cl_LF q = ffloor(x);
	return cl_LF_div_t(cl_LF_to_I(q),LF_LF_minus_LF(x,q));
}
inline const cl_I floor1 (const cl_LF& x)
{
	extern const cl_I cl_LF_to_I (const cl_LF& x);
	return cl_LF_to_I(ffloor(x));
}

// ceiling2(x) liefert (ceiling x), wo x ein LF ist.
inline const cl_LF_div_t ceiling2 (const cl_LF& x)
{
	extern const cl_LF LF_LF_minus_LF (const cl_LF&, const cl_LF&);
	extern const cl_I cl_LF_to_I (const cl_LF& x);
	cl_LF q = fceiling(x);
	return cl_LF_div_t(cl_LF_to_I(q),LF_LF_minus_LF(x,q));
}
inline const cl_I ceiling1 (const cl_LF& x)
{
	extern const cl_I cl_LF_to_I (const cl_LF& x);
	return cl_LF_to_I(fceiling(x));
}

// truncate2(x) liefert (truncate x), wo x ein LF ist.
inline const cl_LF_div_t truncate2 (const cl_LF& x)
{
	extern const cl_LF LF_LF_minus_LF (const cl_LF&, const cl_LF&);
	extern const cl_I cl_LF_to_I (const cl_LF& x);
	cl_LF q = ftruncate(x);
	return cl_LF_div_t(cl_LF_to_I(q),LF_LF_minus_LF(x,q));
}
inline const cl_I truncate1 (const cl_LF& x)
{
	extern const cl_I cl_LF_to_I (const cl_LF& x);
	return cl_LF_to_I(ftruncate(x));
}

// round2(x) liefert (round x), wo x ein LF ist.
inline const cl_LF_div_t round2 (const cl_LF& x)
{
	extern const cl_LF LF_LF_minus_LF (const cl_LF&, const cl_LF&);
	extern const cl_I cl_LF_to_I (const cl_LF& x);
	cl_LF q = fround(x);
	return cl_LF_div_t(cl_LF_to_I(q),LF_LF_minus_LF(x,q));
}
inline const cl_I round1 (const cl_LF& x)
{
	extern const cl_I cl_LF_to_I (const cl_LF& x);
	return cl_LF_to_I(fround(x));
}

// floor2(x,y) liefert (floor x y).
extern const cl_LF_div_t floor2 (const cl_LF& x, const cl_LF& y);
inline const cl_I floor1 (const cl_LF& x, const cl_LF& y) { return floor1(x/y); }

// ceiling2(x,y) liefert (ceiling x y).
extern const cl_LF_div_t ceiling2 (const cl_LF& x, const cl_LF& y);
inline const cl_I ceiling1 (const cl_LF& x, const cl_LF& y) { return ceiling1(x/y); }

// truncate2(x,y) liefert (truncate x y).
extern const cl_LF_div_t truncate2 (const cl_LF& x, const cl_LF& y);
inline const cl_I truncate1 (const cl_LF& x, const cl_LF& y) { return truncate1(x/y); }

// round2(x,y) liefert (round x y).
extern const cl_LF_div_t round2 (const cl_LF& x, const cl_LF& y);
inline const cl_I round1 (const cl_LF& x, const cl_LF& y) { return round1(x/y); }


// cl_float(x,y) returns a long float if y is a long float.
inline const cl_LF cl_float (const cl_F& x, const cl_LF& y)
{
	extern const cl_F cl_float (const cl_F& x, const cl_F& y);
	return The(cl_LF)(cl_float(x,(const cl_F&)y));
}
inline const cl_LF cl_float (const cl_I& x, const cl_LF& y)
{
	extern const cl_F cl_float (const cl_I& x, const cl_F& y);
	return The(cl_LF)(cl_float(x,(const cl_F&)y));
}
inline const cl_LF cl_float (const cl_RA& x, const cl_LF& y)
{
	extern const cl_F cl_float (const cl_RA& x, const cl_F& y);
	return The(cl_LF)(cl_float(x,(const cl_F&)y));
}
inline const cl_LF cl_float (int x, const cl_LF& y)
	{ return cl_float(cl_I(x),y); }
inline const cl_LF cl_float (unsigned int x, const cl_LF& y)
	{ return cl_float(cl_I(x),y); }


// Return type for decode_float:
struct decoded_lfloat {
	cl_LF mantissa;
	cl_I exponent;
	cl_LF sign;
// Constructor.
	decoded_lfloat () {}
	decoded_lfloat (const cl_LF& m, const cl_I& e, const cl_LF& s) : mantissa(m), exponent(e), sign(s) {}
};

// decode_float(x) liefert zu einem Float x: (decode-float x).
// x = 0.0 liefert (0.0, 0, 1.0).
// x = (-1)^s * 2^e * m liefert ((-1)^0 * 2^0 * m, e als Integer, (-1)^s).
extern const decoded_lfloat decode_float (const cl_LF& x);

// float_exponent(x) liefert zu einem Float x:
// den Exponenten von (decode-float x).
// x = 0.0 liefert 0.
// x = (-1)^s * 2^e * m liefert e.
extern sintE float_exponent (const cl_LF& x);

// float_radix(x) liefert (float-radix x), wo x ein Float ist.
inline sintL float_radix (const cl_LF& x)
{
	(void)x; // unused x
	return 2;
}

// float_sign(x) liefert (float-sign x), wo x ein Float ist.
extern const cl_LF float_sign (const cl_LF& x);

// float_digits(x) liefert (float-digits x), wo x ein Float ist.
// < ergebnis: ein uintC >0
extern uintC float_digits (const cl_LF& x);

// float_precision(x) liefert (float-precision x), wo x ein Float ist.
// < ergebnis: ein uintC >=0
extern uintC float_precision (const cl_LF& x);


// integer_decode_float(x) liefert zu einem Float x: (integer-decode-float x).
// x = 0.0 liefert (0, 0, 1).
// x = (-1)^s * 2^e * m bei Float-Precision p liefert
//   (Mantisse 2^p * m als Integer, e-p als Integer, (-1)^s als Fixnum).
extern const cl_idecoded_float integer_decode_float (const cl_LF& x);


// scale_float(x,delta) liefert x*2^delta, wo x ein LF ist.
extern const cl_LF scale_float (const cl_LF& x, sintC delta);
extern const cl_LF scale_float (const cl_LF& x, const cl_I& delta);


// max(x,y) liefert (max x y), wo x und y Floats sind.
extern const cl_LF max (const cl_LF& x, const cl_LF& y);

// min(x,y) liefert (min x y), wo x und y Floats sind.
extern const cl_LF min (const cl_LF& x, const cl_LF& y);

// signum(x) liefert (signum x), wo x ein Float ist.
extern const cl_LF signum (const cl_LF& x);


// Konversion zu einem C "float".
extern float float_approx (const cl_LF& x);

// Konversion zu einem C "double".
extern double double_approx (const cl_LF& x);


// This could be optimized to use in-place operations.
inline cl_LF& operator+= (cl_LF& x, const cl_LF& y) { return x = x + y; }
inline cl_LF& operator++ /* prefix */ (cl_LF& x) { return x = plus1(x); }
inline void operator++ /* postfix */ (cl_LF& x, int dummy) { (void)dummy; x = plus1(x); }
inline cl_LF& operator-= (cl_LF& x, const cl_LF& y) { return x = x - y; }
inline cl_LF& operator-- /* prefix */ (cl_LF& x) { return x = minus1(x); }
inline void operator-- /* postfix */ (cl_LF& x, int dummy) { (void)dummy; x = minus1(x); }
inline cl_LF& operator*= (cl_LF& x, const cl_LF& y) { return x = x * y; }
inline cl_LF& operator/= (cl_LF& x, const cl_LF& y) { return x = x / y; }


// Runtime typing support.
extern cl_class cl_class_lfloat;


// Debugging support.
#ifdef CL_DEBUG
extern int cl_LF_debug_module;
CL_FORCE_LINK(cl_LF_debug_dummy, cl_LF_debug_module)
#endif

}  // namespace cln

#endif /* _CL_LFLOAT_H */
