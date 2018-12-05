// Public float operations.

#ifndef _CL_FLOAT_H
#define _CL_FLOAT_H

#include "cln/number.h"
#include "cln/float_class.h"
#include "cln/floatformat.h"
#include "cln/random.h"
#include "cln/integer_class.h"
#include "cln/sfloat_class.h"
#include "cln/ffloat_class.h"
#include "cln/dfloat_class.h"
#include "cln/lfloat_class.h"
#include "cln/exception.h"

namespace cln {

CL_DEFINE_AS_CONVERSION(cl_F)


// Return type for integer_decode_float:
struct cl_idecoded_float {
	cl_I mantissa;
	cl_I exponent;
	cl_I sign;
// Constructor.
	cl_idecoded_float () {}
	cl_idecoded_float (const cl_I& m, const cl_I& e, const cl_I& s) : mantissa(m), exponent(e), sign(s) {}
};


// zerop(x) testet, ob (= x 0).
extern bool zerop (const cl_F& x);

// minusp(x) testet, ob (< x 0).
extern bool minusp (const cl_F& x);

// plusp(x) testet, ob (> x 0).
extern bool plusp (const cl_F& x);


// cl_F_to_SF(x) wandelt ein Float x in ein Short-Float um und rundet dabei.
extern const cl_SF cl_F_to_SF (const cl_F& x);

// cl_F_to_FF(x) wandelt ein Float x in ein Single-Float um und rundet dabei.
extern const cl_FF cl_F_to_FF (const cl_F& x);

// cl_F_to_DF(x) wandelt ein Float x in ein Double-Float um und rundet dabei.
extern const cl_DF cl_F_to_DF (const cl_F& x);

// cl_F_to_LF(x,len) wandelt ein Float x in ein Long-Float mit len Digits um
// und rundet dabei.
// > uintC len: gewünschte Anzahl Digits, >=LF_minlen
extern const cl_LF cl_F_to_LF (const cl_F& x, uintC len);


// The default float format used when converting rational numbers to floats.
extern float_format_t default_float_format;

// Returns the smallest float format which guarantees at least n decimal digits
// in the mantissa (after the decimal point).
extern float_format_t float_format (uintE n);

// cl_float(x,y) wandelt ein Float x in das Float-Format des Floats y um
// und rundet dabei nötigenfalls.
// > x,y: Floats
// < ergebnis: (float x y)
extern const cl_F cl_float (const cl_F& x, const cl_F& y);

// cl_float(x,f) wandelt ein Float x in das Float-Format f um
// und rundet dabei nötigenfalls.
// > x: ein Float
// > f: eine Float-Format-Spezifikation
// < ergebnis: (float x f)
extern const cl_F cl_float (const cl_F& x, float_format_t f);

// cl_float(x) wandelt eine reelle Zahl x in ein Float um
// und rundet dabei nötigenfalls.
// > x: eine reelle Zahl
// < ergebnis: (float x)
// Abhängig von default_float_format.
inline const cl_F cl_float (const cl_F& x) { return x; }

// cl_float(x,y) wandelt ein Integer x in das Float-Format des Floats y um
// und rundet dabei nötigenfalls.
// > x: ein Integer
// > y: ein Float
// < ergebnis: (float x y)
extern const cl_F cl_float (const cl_I& x, const cl_F& y);

// cl_float(x,y) wandelt ein Integer x in das Float-Format f um
// und rundet dabei nötigenfalls.
// > x: ein Integer
// > f: eine Float-Format-Spezifikation
// < ergebnis: (float x f)
extern const cl_F cl_float (const cl_I& x, float_format_t f);

// cl_float(x) wandelt ein Integer x in ein Float um und rundet dabei.
// > x: ein Integer
// < ergebnis: (float x)
// Abhängig von default_float_format.
extern const cl_F cl_float (const cl_I& x);

// cl_float(x,y) wandelt eine rationale Zahl x in das Float-Format des
// Floats y um und rundet dabei nötigenfalls.
// > x: eine rationale Zahl
// > y: ein Float
// < ergebnis: (float x y)
extern const cl_F cl_float (const cl_RA& x, const cl_F& y);

// cl_float(x,y) wandelt eine rationale Zahl x in das Float-Format f um
// und rundet dabei nötigenfalls.
// > x: eine rationale Zahl
// > f: eine Float-Format-Spezifikation
// < ergebnis: (float x f)
extern const cl_F cl_float (const cl_RA& x, float_format_t f);

// cl_float(x) wandelt eine rationale Zahl x in ein Float um und rundet dabei.
// > x: eine rationale Zahl
// < ergebnis: (float x)
// Abhängig von default_float_format.
extern const cl_F cl_float (const cl_RA& x);

// The C++ compilers are not clever enough to guess this:
inline const cl_F cl_float (int x, const cl_F& y)
	{ return cl_float(cl_I(x),y); }
inline const cl_F cl_float (unsigned int x, const cl_F& y)
	{ return cl_float(cl_I(x),y); }
inline const cl_F cl_float (int x, float_format_t y)
	{ return cl_float(cl_I(x),y); }
inline const cl_F cl_float (unsigned int x, float_format_t y)
	{ return cl_float(cl_I(x),y); }
inline const cl_F cl_float (int x)
	{ return cl_float(cl_I(x)); }
inline const cl_F cl_float (unsigned int x)
	{ return cl_float(cl_I(x)); }
// The C++ compilers could hardly guess the following:
inline const cl_F cl_float (float x, const cl_F& y)
	{ return cl_float(cl_FF(x),y); }
inline const cl_F cl_float (double x, const cl_F& y)
	{ return cl_float(cl_DF(x),y); }
inline const cl_F cl_float (float x, float_format_t y)
	{ return cl_float(cl_FF(x),y); }
inline const cl_F cl_float (double x, float_format_t y)
	{ return cl_float(cl_DF(x),y); }
inline const cl_F cl_float (float x)
	{ return cl_float(cl_FF(x)); }
inline const cl_F cl_float (double x)
	{ return cl_float(cl_DF(x)); }


// Liefert (- x), wo x ein Float ist.
extern const cl_F operator- (const cl_F& x);

// Liefert (+ x y), wo x und y Floats sind.
extern const cl_F operator+ (const cl_F& x, const cl_F& y);
// The C++ compilers could hardly guess the following:
inline const cl_F operator+ (const cl_RA& x, const cl_F& y)
	{ return cl_float(x,y) + y; }
inline const cl_F operator+ (const cl_I& x, const cl_F& y)
	{ return cl_float(x,y) + y; }
inline const cl_F operator+ (const cl_F& x, const cl_RA& y)
	{ return x + cl_float(y,x); }
inline const cl_F operator+ (const cl_F& x, const cl_I& y)
	{ return x + cl_float(y,x); }
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_F operator+ (const int x, const cl_F& y)
	{ return cl_I(x) + y; }
inline const cl_F operator+ (const unsigned int x, const cl_F& y)
	{ return cl_I(x) + y; }
inline const cl_F operator+ (const long x, const cl_F& y)
	{ return cl_I(x) + y; }
inline const cl_F operator+ (const unsigned long x, const cl_F& y)
	{ return cl_I(x) + y; }
#ifdef HAVE_LONGLONG
inline const cl_F operator+ (const long long x, const cl_F& y)
	{ return cl_I(x) + y; }
inline const cl_F operator+ (const unsigned long long x, const cl_F& y)
	{ return cl_I(x) + y; }
#endif
inline const cl_F operator+ (const float x, const cl_F& y)
	{ return cl_F(x) + y; }
inline const cl_F operator+ (const double x, const cl_F& y)
	{ return cl_F(x) + y; }
inline const cl_F operator+ (const cl_F& x, const int y)
	{ return x + cl_I(y); }
inline const cl_F operator+ (const cl_F& x, const unsigned int y)
	{ return x + cl_I(y); }
inline const cl_F operator+ (const cl_F& x, const long y)
	{ return x + cl_I(y); }
inline const cl_F operator+ (const cl_F& x, const unsigned long y)
	{ return x + cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_F operator+ (const cl_F& x, const long long y)
	{ return x + cl_I(y); }
inline const cl_F operator+ (const cl_F& x, const unsigned long long y)
	{ return x + cl_I(y); }
#endif
inline const cl_F operator+ (const cl_F& x, const float y)
	{ return x + cl_F(y); }
inline const cl_F operator+ (const cl_F& x, const double y)
	{ return x + cl_F(y); }

// Liefert (- x y), wo x und y Floats sind.
extern const cl_F operator- (const cl_F& x, const cl_F& y);
// The C++ compilers could hardly guess the following:
inline const cl_F operator- (const cl_RA& x, const cl_F& y)
	{ return cl_float(x,y) - y; }
inline const cl_F operator- (const cl_I& x, const cl_F& y)
	{ return cl_float(x,y) - y; }
inline const cl_F operator- (const cl_F& x, const cl_RA& y)
	{ return x - cl_float(y,x); }
inline const cl_F operator- (const cl_F& x, const cl_I& y)
	{ return x - cl_float(y,x); }
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_F operator- (const int x, const cl_F& y)
	{ return cl_I(x) - y; }
inline const cl_F operator- (const unsigned int x, const cl_F& y)
	{ return cl_I(x) - y; }
inline const cl_F operator- (const long x, const cl_F& y)
	{ return cl_I(x) - y; }
inline const cl_F operator- (const unsigned long x, const cl_F& y)
	{ return cl_I(x) - y; }
#ifdef HAVE_LONGLONG
inline const cl_F operator- (const long long x, const cl_F& y)
	{ return cl_I(x) - y; }
inline const cl_F operator- (const unsigned long long x, const cl_F& y)
	{ return cl_I(x) - y; }
#endif
inline const cl_F operator- (const float x, const cl_F& y)
	{ return cl_F(x) - y; }
inline const cl_F operator- (const double x, const cl_F& y)
	{ return cl_F(x) - y; }
inline const cl_F operator- (const cl_F& x, const int y)
	{ return x - cl_I(y); }
inline const cl_F operator- (const cl_F& x, const unsigned int y)
	{ return x - cl_I(y); }
inline const cl_F operator- (const cl_F& x, const long y)
	{ return x - cl_I(y); }
inline const cl_F operator- (const cl_F& x, const unsigned long y)
	{ return x - cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_F operator- (const cl_F& x, const long long y)
	{ return x - cl_I(y); }
inline const cl_F operator- (const cl_F& x, const unsigned long long y)
	{ return x - cl_I(y); }
#endif
inline const cl_F operator- (const cl_F& x, const float y)
	{ return x - cl_F(y); }
inline const cl_F operator- (const cl_F& x, const double y)
	{ return x - cl_F(y); }

// Liefert (* x y), wo x und y Floats sind.
extern const cl_F operator* (const cl_F& x, const cl_F& y);
// Spezialfall x oder y Integer oder rationale Zahl.
inline const cl_R operator* (const cl_F& x, const cl_I& y)
{
	extern const cl_R cl_F_I_mul (const cl_F&, const cl_I&);
	return cl_F_I_mul(x,y);
}
inline const cl_R operator* (const cl_I& x, const cl_F& y)
{
	extern const cl_R cl_F_I_mul (const cl_F&, const cl_I&);
	return cl_F_I_mul(y,x);
}
inline const cl_R operator* (const cl_F& x, const cl_RA& y)
{
	extern const cl_R cl_F_RA_mul (const cl_F&, const cl_RA&);
	return cl_F_RA_mul(x,y);
}
inline const cl_R operator* (const cl_RA& x, const cl_F& y)
{
	extern const cl_R cl_F_RA_mul (const cl_F&, const cl_RA&);
	return cl_F_RA_mul(y,x);
}
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R operator* (const int x, const cl_F& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const unsigned int x, const cl_F& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const long x, const cl_F& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const unsigned long x, const cl_F& y)
	{ return cl_I(x) * y; }
#ifdef HAVE_LONGLONG
inline const cl_R operator* (const long long x, const cl_F& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const unsigned long long x, const cl_F& y)
	{ return cl_I(x) * y; }
#endif
inline const cl_F operator* (const float x, const cl_F& y)
	{ return cl_F(x) * y; }
inline const cl_F operator* (const double x, const cl_F& y)
	{ return cl_F(x) * y; }
inline const cl_R operator* (const cl_F& x, const int y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_F& x, const unsigned int y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_F& x, const long y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_F& x, const unsigned long y)
	{ return x * cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_R operator* (const cl_F& x, const long long y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_F& x, const unsigned long long y)
	{ return x * cl_I(y); }
#endif
inline const cl_F operator* (const cl_F& x, const float y)
	{ return x * cl_F(y); }
inline const cl_F operator* (const cl_F& x, const double y)
	{ return x * cl_F(y); }

// Liefert (* x x), wo x ein Float ist.
extern const cl_F square (const cl_F& x);

// Liefert (/ x y), wo x und y Floats sind.
extern const cl_F operator/ (const cl_F& x, const cl_F& y);
// Liefert (/ x y), wo x und y ein Float und eine rationale Zahl sind.
extern const cl_F operator/ (const cl_F& x, const cl_RA& y);
extern const cl_F operator/ (const cl_F& x, const cl_I& y);
extern const cl_R operator/ (const cl_RA& x, const cl_F& y);
extern const cl_R operator/ (const cl_I& x, const cl_F& y);
// The C++ compilers could hardly guess the following:
inline const cl_F operator/ (const cl_F& x, const int y)
	{ return x / cl_I(y); }
inline const cl_F operator/ (const cl_F& x, const unsigned int y)
	{ return x / cl_I(y); }
inline const cl_F operator/ (const cl_F& x, const long y)
	{ return x / cl_I(y); }
inline const cl_F operator/ (const cl_F& x, const unsigned long y)
	{ return x / cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_F operator/ (const cl_F& x, const long long y)
	{ return x / cl_I(y); }
inline const cl_F operator/ (const cl_F& x, const unsigned long long y)
	{ return x / cl_I(y); }
#endif
inline const cl_F operator/ (const cl_F& x, const float y)
	{ return x / cl_F(y); }
inline const cl_F operator/ (const cl_F& x, const double y)
	{ return x / cl_F(y); }
inline const cl_R operator/ (const int x, const cl_F& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const unsigned int x, const cl_F& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const long x, const cl_F& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const unsigned long x, const cl_F& y)
	{ return cl_I(x) / y; }
#ifdef HAVE_LONGLONG
inline const cl_R operator/ (const long long x, const cl_F& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const unsigned long long x, const cl_F& y)
	{ return cl_I(x) / y; }
#endif
inline const cl_F operator/ (const float x, const cl_F& y)
	{ return cl_F(x) / y; }
inline const cl_F operator/ (const double x, const cl_F& y)
	{ return cl_F(x) / y; }

// Liefert (abs x), wo x ein Float ist.
extern const cl_F abs (const cl_F& x);

// Liefert zu einem Float x>=0 : (sqrt x), ein Float.
extern const cl_F sqrt (const cl_F& x);

// recip(x) liefert (/ x), wo x ein Float ist.
extern const cl_F recip (const cl_F& x);

// (1+ x), wo x ein Float ist.
inline const cl_F plus1 (const cl_F& x) // { return x + cl_I(1); }
{
	return x + cl_float(1,x);
}

// (1- x), wo x ein Float ist.
inline const cl_F minus1 (const cl_F& x) // { return x + cl_I(-1); }
{
	return x + cl_float(-1,x);
}

// compare(x,y) vergleicht zwei Floats x und y.
// Ergebnis: 0 falls x=y, +1 falls x>y, -1 falls x<y.
extern cl_signean compare (const cl_F& x, const cl_F& y);

// equal_hashcode(x) liefert einen equal-invarianten Hashcode für x.
extern uint32 equal_hashcode (const cl_F& x);

inline bool operator== (const cl_F& x, const cl_F& y)
	{ return compare(x,y)==0; }
inline bool operator!= (const cl_F& x, const cl_F& y)
	{ return compare(x,y)!=0; }
inline bool operator<= (const cl_F& x, const cl_F& y)
	{ return compare(x,y)<=0; }
inline bool operator< (const cl_F& x, const cl_F& y)
	{ return compare(x,y)<0; }
inline bool operator>= (const cl_F& x, const cl_F& y)
	{ return compare(x,y)>=0; }
inline bool operator> (const cl_F& x, const cl_F& y)
	{ return compare(x,y)>0; }


// ffloor(x) liefert (ffloor x), wo x ein Float ist.
extern const cl_F ffloor (const cl_F& x);

// fceiling(x) liefert (fceiling x), wo x ein Float ist.
extern const cl_F fceiling (const cl_F& x);

// ftruncate(x) liefert (ftruncate x), wo x ein Float ist.
extern const cl_F ftruncate (const cl_F& x);

// fround(x) liefert (fround x), wo x ein Float ist.
extern const cl_F fround (const cl_F& x);


// Return type for frounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_F_fdiv_t {
	cl_F quotient;
	cl_F remainder;
// Constructor.
	cl_F_fdiv_t () {}
	cl_F_fdiv_t (const cl_F& q, const cl_F& r) : quotient(q), remainder(r) {}
};

// ffloor2(x) liefert (ffloor x), wo x ein F ist.
extern const cl_F_fdiv_t ffloor2 (const cl_F& x);

// fceiling2(x) liefert (fceiling x), wo x ein F ist.
extern const cl_F_fdiv_t fceiling2 (const cl_F& x);

// ftruncate2(x) liefert (ftruncate x), wo x ein F ist.
extern const cl_F_fdiv_t ftruncate2 (const cl_F& x);

// fround2(x) liefert (fround x), wo x ein F ist.
extern const cl_F_fdiv_t fround2 (const cl_F& x);


// Return type for rounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_F_div_t {
	cl_I quotient;
	cl_F remainder;
// Constructor.
	cl_F_div_t () {}
	cl_F_div_t (const cl_I& q, const cl_F& r) : quotient(q), remainder(r) {}
};

// floor2(x) liefert (floor x), wo x ein F ist.
extern const cl_F_div_t floor2 (const cl_F& x);
extern const cl_I floor1 (const cl_F& x);

// ceiling2(x) liefert (ceiling x), wo x ein F ist.
extern const cl_F_div_t ceiling2 (const cl_F& x);
extern const cl_I ceiling1 (const cl_F& x);

// truncate2(x) liefert (truncate x), wo x ein F ist.
extern const cl_F_div_t truncate2 (const cl_F& x);
extern const cl_I truncate1 (const cl_F& x);

// round2(x) liefert (round x), wo x ein F ist.
extern const cl_F_div_t round2 (const cl_F& x);
extern const cl_I round1 (const cl_F& x);

// floor2(x,y) liefert (floor x y), wo x und y Floats sind.
extern const cl_F_div_t floor2 (const cl_F& x, const cl_F& y);
inline const cl_I floor1 (const cl_F& x, const cl_F& y) { return floor1(x/y); }

// ceiling2(x,y) liefert (ceiling x y), wo x und y Floats sind.
extern const cl_F_div_t ceiling2 (const cl_F& x, const cl_F& y);
inline const cl_I ceiling1 (const cl_F& x, const cl_F& y) { return ceiling1(x/y); }

// truncate2(x,y) liefert (truncate x y), wo x und y Floats sind.
extern const cl_F_div_t truncate2 (const cl_F& x, const cl_F& y);
inline const cl_I truncate1 (const cl_F& x, const cl_F& y) { return truncate1(x/y); }

// round2(x,y) liefert (round x y), wo x und y Floats sind.
extern const cl_F_div_t round2 (const cl_F& x, const cl_F& y);
inline const cl_I round1 (const cl_F& x, const cl_F& y) { return round1(x/y); }


// Return type for decode_float:
struct decoded_float {
	cl_F mantissa;
	cl_I exponent;
	cl_F sign;
// Constructor.
	decoded_float () {}
	decoded_float (const cl_F& m, const cl_I& e, const cl_F& s) : mantissa(m), exponent(e), sign(s) {}
};

// decode_float(x) liefert zu einem Float x: (decode-float x).
// x = 0.0 liefert (0.0, 0, 1.0).
// x = (-1)^s * 2^e * m liefert ((-1)^0 * 2^0 * m, e als Integer, (-1)^s).
extern const decoded_float decode_float (const cl_F& x);

// float_exponent(x) liefert zu einem Float x:
// den Exponenten von (decode-float x).
// x = 0.0 liefert 0.
// x = (-1)^s * 2^e * m liefert e.
extern sintE float_exponent (const cl_F& x);

// float_radix(x) liefert (float-radix x), wo x ein Float ist.
inline sintL float_radix (const cl_F& x)
{
	(void)x; // unused x
	return 2;
}

// float_sign(x) liefert (float-sign x), wo x ein Float ist.
extern const cl_F float_sign (const cl_F& x);

// float_sign(x,y) liefert (float-sign x y), wo x und y Floats sind.
extern const cl_F float_sign (const cl_F& x, const cl_F& y);

// float_digits(x) liefert (float-digits x), wo x ein Float ist.
// < ergebnis: ein uintC >0
extern uintC float_digits (const cl_F& x);

// float_precision(x) liefert (float-precision x), wo x ein Float ist.
// < ergebnis: ein uintC >=0
extern uintC float_precision (const cl_F& x);

// Returns the floating point format of a float.
inline float_format_t float_format (const cl_F& x)
	{ return (float_format_t) float_digits(x); }


// integer_decode_float(x) liefert zu einem Float x: (integer-decode-float x).
// x = 0.0 liefert (0, 0, 1).
// x = (-1)^s * 2^e * m bei Float-Precision p liefert
//   (Mantisse 2^p * m als Integer, e-p als Integer, (-1)^s als Fixnum).
extern const cl_idecoded_float integer_decode_float (const cl_F& x);


// rational(x) liefert (rational x), wo x ein Float ist.
extern const cl_RA rational (const cl_F& x);


// scale_float(x,delta) liefert x*2^delta, wo x ein Float ist.
extern const cl_F scale_float (const cl_F& x, sintC delta);
extern const cl_F scale_float (const cl_F& x, const cl_I& delta);


// max(x,y) liefert (max x y), wo x und y Floats sind.
extern const cl_F max (const cl_F& x, const cl_F& y);

// min(x,y) liefert (min x y), wo x und y Floats sind.
extern const cl_F min (const cl_F& x, const cl_F& y);

// signum(x) liefert (signum x), wo x ein Float ist.
extern const cl_F signum (const cl_F& x);


// Returns the largest (most positive) floating point number in float format f.
extern const cl_F most_positive_float (float_format_t f);

// Returns the smallest (most negative) floating point number in float format f.
extern const cl_F most_negative_float (float_format_t f);

// Returns the least positive floating point number (i.e. > 0 but closest to 0)
// in float format f.
extern const cl_F least_positive_float (float_format_t f);

// Returns the least negative floating point number (i.e. < 0 but closest to 0)
// in float format f.
extern const cl_F least_negative_float (float_format_t f);

// Returns the smallest floating point number e > 0 such that 1+e != 1.
extern const cl_F float_epsilon (float_format_t f);

// Returns the smallest floating point number e > 0 such that 1-e != 1.
extern const cl_F float_negative_epsilon (float_format_t f);


// Konversion zu einem C "float".
extern float float_approx (const cl_F& x);

// Konversion zu einem C "double".
extern double double_approx (const cl_F& x);


// Transcendental functions


// pi(y) liefert die Zahl pi im selben Float-Format wie y.
// > y: ein Float
extern const cl_F pi (const cl_F& y);

// pi(y) liefert die Zahl pi im Float-Format f.
// > f: eine Float-Format-Spezifikation
extern const cl_F pi (float_format_t f);

// pi() liefert die Zahl pi im Default-Float-Format.
extern const cl_F pi (void);


// sin(x) liefert den Sinus (sin x) eines Float x.
extern const cl_F sin (const cl_F& x);

// cos(x) liefert den Cosinus (cos x) eines Float x.
extern const cl_F cos (const cl_F& x);

// Return type for cos_sin():
struct cos_sin_t {
	cl_R cos;
	cl_R sin;
// Constructor:
	cos_sin_t () {}
	cos_sin_t (const cl_R& u, const cl_R& v) : cos (u), sin (v) {}
};

// cos_sin(x) liefert ((cos x),(sin x)), beide Werte.
extern const cos_sin_t cos_sin (const cl_F& x);

// tan(x) liefert den Tangens (tan x) eines Float x.
extern const cl_F tan (const cl_F& x);


// exp1(y) liefert die Zahl e = exp(1) im selben Float-Format wie y.
// > y: ein Float
extern const cl_F exp1 (const cl_F& y);

// exp1(y) liefert die Zahl e = exp(1) im Float-Format f.
// > f: eine Float-Format-Spezifikation
extern const cl_F exp1 (float_format_t f);

// exp1() liefert die Zahl e = exp(1) im Default-Float-Format.
extern const cl_F exp1 (void);


// ln(x) liefert zu einem Float x>0 die Zahl ln(x).
extern const cl_F ln (const cl_F& x);
// Spezialfall: x Long-Float -> Ergebnis Long-Float
inline const cl_LF ln (const cl_LF& x) { return The(cl_LF)(ln(The(cl_F)(x))); }

// exp(x) liefert zu einem Float x die Zahl exp(x).
extern const cl_F exp (const cl_F& x);

// sinh(x) liefert zu einem Float x die Zahl sinh(x).
extern const cl_F sinh (const cl_F& x);

// cosh(x) liefert zu einem Float x die Zahl cosh(x).
extern const cl_F cosh (const cl_F& x);

// Return type for cosh_sinh():
struct cosh_sinh_t {
	cl_R cosh;
	cl_R sinh;
// Constructor:
	cosh_sinh_t () {}
	cosh_sinh_t (const cl_R& u, const cl_R& v) : cosh (u), sinh (v) {}
};

// cosh_sinh(x) liefert ((cosh x),(sinh x)), beide Werte.
extern const cosh_sinh_t cosh_sinh (const cl_F& x);

// tanh(x) liefert zu einem Float x die Zahl tanh(x).
extern const cl_F tanh (const cl_F& x);


// eulerconst(y) liefert die Eulersche Konstante
// im selben Float-Format wie y.
// > y: ein Float
extern const cl_F eulerconst (const cl_F& y);

// eulerconst(y) liefert die Eulersche Konstante im Float-Format f.
// > f: eine Float-Format-Spezifikation
extern const cl_F eulerconst (float_format_t f);

// eulerconst() liefert die Eulersche Konstante im Default-Float-Format.
extern const cl_F eulerconst (void);

// catalanconst(y) liefert die Catalansche Konstante
// im selben Float-Format wie y.
// > y: ein Float
extern const cl_F catalanconst (const cl_F& y);

// catalanconst(y) liefert die Catalansche Konstante im Float-Format f.
// > f: eine Float-Format-Spezifikation
extern const cl_F catalanconst (float_format_t f);

// catalanconst() liefert die Catalansche Konstante im Default-Float-Format.
extern const cl_F catalanconst (void);

// zeta(s) returns the Riemann zeta function at s>1.
extern const cl_F zeta (int s, const cl_F& y);
extern const cl_F zeta (int s, float_format_t f);
extern const cl_F zeta (int s);


// random_F(randomstate,n) liefert zu einem Float n>0 ein zufälliges
// Float x mit 0 <= x < n.
// > randomstate: ein Random-State, wird verändert
extern const cl_F random_F (random_state& randomstate, const cl_F& n);

inline const cl_F random_F (const cl_F& n)
	{ return random_F(default_random_state,n); }


// This could be optimized to use in-place operations.
inline cl_F& operator+= (cl_F& x, const cl_F& y) { return x = x + y; }
inline cl_F& operator+= (cl_F& x, const float y) { return x = x + y; }
inline cl_F& operator+= (cl_F& x, const double y) { return x = x + y; }
inline cl_F& operator++ /* prefix */ (cl_F& x) { return x = plus1(x); }
inline void operator++ /* postfix */ (cl_F& x, int dummy) { (void)dummy; x = plus1(x); }
inline cl_F& operator-= (cl_F& x, const cl_F& y) { return x = x - y; }
inline cl_F& operator-= (cl_F& x, const float y) { return x = x - y; }
inline cl_F& operator-= (cl_F& x, const double y) { return x = x - y; }
inline cl_F& operator-- /* prefix */ (cl_F& x) { return x = minus1(x); }
inline void operator-- /* postfix */ (cl_F& x, int dummy) { (void)dummy; x = minus1(x); }
inline cl_F& operator*= (cl_F& x, const cl_F& y) { return x = x * y; }
inline cl_F& operator*= (cl_F& x, const float y) { return x = x * y; }
inline cl_F& operator*= (cl_F& x, const double y) { return x = x * y; }
inline cl_F& operator/= (cl_F& x, const cl_F& y) { return x = x / y; }
inline cl_F& operator/= (cl_F& x, const float y) { return x = x / y; }
inline cl_F& operator/= (cl_F& x, const double y) { return x = x / y; }

// Thrown when a floating-point exception occurs.
class floating_point_exception : public runtime_exception {
public:
	explicit floating_point_exception(const std::string & what)
		: runtime_exception(what) {}
};

// Thrown when NaN occurs.
class floating_point_nan_exception : public floating_point_exception {
public:
	floating_point_nan_exception();
};

// Thrown when overflow occurs.
class floating_point_overflow_exception : public floating_point_exception {
public:
	floating_point_overflow_exception();
};

// Thrown when underflow occurs.
class floating_point_underflow_exception : public floating_point_exception {
public:
	floating_point_underflow_exception();
};




// If this is true, floating point underflow returns zero instead of throwing an exception.
extern bool cl_inhibit_floating_point_underflow;

}  // namespace cln

#endif /* _CL_FLOAT_H */
