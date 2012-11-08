// Public real number operations.

#ifndef _CL_REAL_H
#define _CL_REAL_H

#include "cln/number.h"
#include "cln/real_class.h"
#include "cln/rational_class.h"
#include "cln/integer_class.h"
#include "cln/float.h"
#include "cln/floatformat.h"
#include "cln/random.h"

namespace cln {

CL_DEFINE_AS_CONVERSION(cl_R)


// zerop(x) testet, ob (= x 0).
extern bool zerop (const cl_R& x);

// minusp(x) testet, ob (< x 0).
extern bool minusp (const cl_R& x);

// plusp(x) testet, ob (> x 0).
extern bool plusp (const cl_R& x);


// R_to_SF(x) wandelt eine reelle Zahl x in ein Short-Float um.
// < ergebnis: (coerce x 'short-float)
extern const cl_SF cl_R_to_SF (const cl_R& x);

// R_to_FF(x) wandelt eine reelle Zahl x in ein Single-Float um.
// < ergebnis: (coerce x 'single-float)
extern const cl_FF cl_R_to_FF (const cl_R& x);

// R_to_DF(x) wandelt eine reelle Zahl x in ein Double-Float um.
// < ergebnis: (coerce x 'double-float)
extern const cl_DF cl_R_to_DF (const cl_R& x);

// R_to_LF(x,len) wandelt eine reelle Zahl x in ein Long-Float mit len Digits um.
// > uintC len: gewünschte Anzahl Digits, >=LF_minlen
// < ergebnis: (coerce x `(long-float ,len))
extern const cl_LF cl_R_to_LF (const cl_R& x, uintC len);

// cl_float(x,y) wandelt eine reelle Zahl x in das Float-Format des
// Floats y um und rundet dabei nötigenfalls.
// > x: eine reelle Zahl
// > y: ein Float
// < ergebnis: (float x y)
extern const cl_F cl_float (const cl_R& x, const cl_F& y);

// cl_float(x,f) wandelt eine reelle Zahl x in das Float-Format f um
// und rundet dabei nötigenfalls.
// > x: eine reelle Zahl
// > f: eine Float-Format-Spezifikation
// < ergebnis: (float x f)
extern const cl_F cl_float (const cl_R& x, float_format_t f);

// cl_float(x) wandelt eine reelle Zahl x in ein Float um
// und rundet dabei nötigenfalls.
// > x: eine reelle Zahl
// < ergebnis: (float x)
// Abhängig von default_float_format.
extern const cl_F cl_float (const cl_R& x);


// Liefert (- x), wo x eine reelle Zahl ist.
extern const cl_R operator- (const cl_R& x);

// Liefert (+ x y), wo x und y reelle Zahlen sind.
extern const cl_R operator+ (const cl_R& x, const cl_R& y);
// Spezialfall: x oder y Float -> Ergebnis Float
inline const cl_F operator+ (const cl_R& x, const cl_F& y)
	{ return The(cl_F)(x + The(cl_R)(y)); }
inline const cl_F operator+ (const cl_F& x, const cl_R& y)
	{ return The(cl_F)(The(cl_R)(x) + y); }
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R operator+ (const int x, const cl_R& y)
	{ return cl_I(x) + y; }
inline const cl_R operator+ (const unsigned int x, const cl_R& y)
	{ return cl_I(x) + y; }
inline const cl_R operator+ (const long x, const cl_R& y)
	{ return cl_I(x) + y; }
inline const cl_R operator+ (const unsigned long x, const cl_R& y)
	{ return cl_I(x) + y; }
#ifdef HAVE_LONGLONG
inline const cl_R operator+ (const long long x, const cl_R& y)
	{ return cl_I(x) + y; }
inline const cl_R operator+ (const unsigned long long x, const cl_R& y)
	{ return cl_I(x) + y; }
#endif
inline const cl_F operator+ (const float x, const cl_R& y)
	{ return The(cl_F)(cl_R(x) + y); }
inline const cl_F operator+ (const double x, const cl_R& y)
	{ return The(cl_F)(cl_R(x) + y); }
inline const cl_R operator+ (const cl_R& x, const int y)
	{ return x + cl_I(y); }
inline const cl_R operator+ (const cl_R& x, const unsigned int y)
	{ return x + cl_I(y); }
inline const cl_R operator+ (const cl_R& x, const long y)
	{ return x + cl_I(y); }
inline const cl_R operator+ (const cl_R& x, const unsigned long y)
	{ return x + cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_R operator+ (const cl_R& x, const long long y)
	{ return x + cl_I(y); }
inline const cl_R operator+ (const cl_R& x, const unsigned long long y)
	{ return x + cl_I(y); }
#endif
inline const cl_F operator+ (const cl_R& x, const float y)
	{ return The(cl_F)(x + cl_R(y)); }
inline const cl_F operator+ (const cl_R& x, const double y)
	{ return The(cl_F)(x + cl_R(y)); }

// Liefert (- x y), wo x und y reelle Zahlen sind.
extern const cl_R operator- (const cl_R& x, const cl_R& y);
// Spezialfall: x oder y Float -> Ergebnis Float
inline const cl_F operator- (const cl_R& x, const cl_F& y)
	{ return The(cl_F)(x - The(cl_R)(y)); }
inline const cl_F operator- (const cl_F& x, const cl_R& y)
	{ return The(cl_F)(The(cl_R)(x) - y); }
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R operator- (const int x, const cl_R& y)
	{ return cl_I(x) - y; }
inline const cl_R operator- (const unsigned int x, const cl_R& y)
	{ return cl_I(x) - y; }
inline const cl_R operator- (const long x, const cl_R& y)
	{ return cl_I(x) - y; }
inline const cl_R operator- (const unsigned long x, const cl_R& y)
	{ return cl_I(x) - y; }
#ifdef HAVE_LONGLONG
inline const cl_R operator- (const long long x, const cl_R& y)
	{ return cl_I(x) - y; }
inline const cl_R operator- (const unsigned long long x, const cl_R& y)
	{ return cl_I(x) - y; }
#endif
inline const cl_F operator- (const float x, const cl_R& y)
	{ return The(cl_F)(cl_R(x) - y); }
inline const cl_F operator- (const double x, const cl_R& y)
	{ return The(cl_F)(cl_R(x) - y); }
inline const cl_R operator- (const cl_R& x, const int y)
	{ return x - cl_I(y); }
inline const cl_R operator- (const cl_R& x, const unsigned int y)
	{ return x - cl_I(y); }
inline const cl_R operator- (const cl_R& x, const long y)
	{ return x - cl_I(y); }
inline const cl_R operator- (const cl_R& x, const unsigned long y)
	{ return x - cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_R operator- (const cl_R& x, const long long y)
	{ return x - cl_I(y); }
inline const cl_R operator- (const cl_R& x, const unsigned long long y)
	{ return x - cl_I(y); }
#endif
inline const cl_F operator- (const cl_R& x, const float y)
	{ return The(cl_F)(x - cl_R(y)); }
inline const cl_F operator- (const cl_R& x, const double y)
	{ return The(cl_F)(x - cl_R(y)); }

// Liefert (* x y), wo x und y reelle Zahlen sind.
extern const cl_R operator* (const cl_R& x, const cl_R& y);
// Dem C++-Compiler muß man auch das Folgende sagen (wg. `int * cl_F' u.ä.):
inline const cl_R operator* (const int x, const cl_R& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const unsigned int x, const cl_R& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const long x, const cl_R& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const unsigned long x, const cl_R& y)
	{ return cl_I(x) * y; }
#ifdef HAVE_LONGLONG
inline const cl_R operator* (const long long x, const cl_R& y)
	{ return cl_I(x) * y; }
inline const cl_R operator* (const unsigned long long x, const cl_R& y)
	{ return cl_I(x) * y; }
#endif
inline const cl_R operator* (const float x, const cl_R& y)
	{ return cl_R(x) * y; }
inline const cl_R operator* (const double x, const cl_R& y)
	{ return cl_R(x) * y; }
inline const cl_R operator* (const cl_R& x, const int y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_R& x, const unsigned int y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_R& x, const long y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_R& x, const unsigned long y)
	{ return x * cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_R operator* (const cl_R& x, const long long y)
	{ return x * cl_I(y); }
inline const cl_R operator* (const cl_R& x, const unsigned long long y)
	{ return x * cl_I(y); }
#endif
inline const cl_R operator* (const cl_R& x, const float y)
	{ return x * cl_R(y); }
inline const cl_R operator* (const cl_R& x, const double y)
	{ return x * cl_R(y); }

// Liefert (* x x), wo x eine reelle Zahl ist.
extern const cl_R square (const cl_R& x);

// Liefert (/ x y), wo x und y reelle Zahlen sind.
extern const cl_R operator/ (const cl_R& x, const cl_R& y);
// Spezialfall: x Float -> Ergebnis Float
inline const cl_F operator/ (const cl_F& x, const cl_R& y)
	{ return The(cl_F)(The(cl_R)(x) / y); }
// Dem C++-Compiler muß man auch das Folgende sagen (wg. `int / cl_F' u.ä.):
inline const cl_R operator/ (const int x, const cl_R& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const unsigned int x, const cl_R& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const long x, const cl_R& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const unsigned long x, const cl_R& y)
	{ return cl_I(x) / y; }
#ifdef HAVE_LONGLONG
inline const cl_R operator/ (const long long x, const cl_R& y)
	{ return cl_I(x) / y; }
inline const cl_R operator/ (const unsigned long long x, const cl_R& y)
	{ return cl_I(x) / y; }
#endif
inline const cl_F operator/ (const float x, const cl_R& y)
	{ return The(cl_F)(cl_R(x) / y); }
inline const cl_F operator/ (const double x, const cl_R& y)
	{ return The(cl_F)(cl_R(x) / y); }
inline const cl_R operator/ (const cl_R& x, const int y)
	{ return x / cl_I(y); }
inline const cl_R operator/ (const cl_R& x, const unsigned int y)
	{ return x / cl_I(y); }
inline const cl_R operator/ (const cl_R& x, const long y)
	{ return x / cl_I(y); }
inline const cl_R operator/ (const cl_R& x, const unsigned long y)
	{ return x / cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_R operator/ (const cl_R& x, const long long y)
	{ return x / cl_I(y); }
inline const cl_R operator/ (const cl_R& x, const unsigned long long y)
	{ return x / cl_I(y); }
#endif
inline const cl_R operator/ (const cl_R& x, const float y)
	{ return x / cl_R(y); }
inline const cl_R operator/ (const cl_R& x, const double y)
	{ return x / cl_R(y); }

// Liefert (abs x), wo x eine reelle Zahl ist.
extern const cl_R abs (const cl_R& x);

// recip(x) liefert (/ x), wo x eine reelle Zahl ist.
extern const cl_R recip (const cl_R& x);

// (1+ x), wo x eine reelle Zahl ist.
extern const cl_R plus1 (const cl_R& x);

// (1- x), wo x eine reelle Zahl ist.
extern const cl_R minus1 (const cl_R& x);


// Return type for rounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_R_div_t {
	cl_I quotient;
	cl_R remainder;
// Constructor.
	cl_R_div_t () {}
	cl_R_div_t (const cl_I& q, const cl_R& r) : quotient(q), remainder(r) {}
	cl_R_div_t (const struct cl_I_div_t &);
	cl_R_div_t (const struct cl_RA_div_t &);
	cl_R_div_t (const struct cl_F_div_t &);
};

// floor2(x) liefert (floor x), wo x eine reelle Zahl ist.
extern const cl_R_div_t floor2 (const cl_R& x);
extern const cl_I floor1 (const cl_R& x);

// ceiling2(x) liefert (ceiling x), wo x eine reelle Zahl ist.
extern const cl_R_div_t ceiling2 (const cl_R& x);
extern const cl_I ceiling1 (const cl_R& x);

// truncate2(x) liefert (truncate x), wo x eine reelle Zahl ist.
extern const cl_R_div_t truncate2 (const cl_R& x);
extern const cl_I truncate1 (const cl_R& x);

// round2(x) liefert (round x), wo x eine reelle Zahl ist.
extern const cl_R_div_t round2 (const cl_R& x);
extern const cl_I round1 (const cl_R& x);

// floor2(x,y) liefert (floor x y), wo x und y reelle Zahlen sind.
extern const cl_R_div_t floor2 (const cl_R& x, const cl_R& y);
extern const cl_I floor1 (const cl_R& x, const cl_R& y);

// ceiling2(x,y) liefert (ceiling x y), wo x und y reelle Zahlen sind.
extern const cl_R_div_t ceiling2 (const cl_R& x, const cl_R& y);
extern const cl_I ceiling1 (const cl_R& x, const cl_R& y);

// truncate2(x,y) liefert (truncate x y), wo x und y reelle Zahlen sind.
extern const cl_R_div_t truncate2 (const cl_R& x, const cl_R& y);
extern const cl_I truncate1 (const cl_R& x, const cl_R& y);

// round2(x,y) liefert (round x y), wo x und y reelle Zahlen sind.
extern const cl_R_div_t round2 (const cl_R& x, const cl_R& y);
extern const cl_I round1 (const cl_R& x, const cl_R& y);


// Return type for frounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_R_fdiv_t {
	cl_F quotient;
	cl_R remainder;
// Constructor.
	cl_R_fdiv_t () {}
	cl_R_fdiv_t (const cl_F& q, const cl_R& r) : quotient(q), remainder(r) {}
	cl_R_fdiv_t (const struct cl_F_fdiv_t &);
};

// ffloor2(x) liefert (ffloor x), wo x eine reelle Zahl ist.
extern const cl_R_fdiv_t ffloor2 (const cl_R& x);
extern const cl_F ffloor (const cl_R& x);

// fceiling2(x) liefert (fceiling x), wo x eine reelle Zahl ist.
extern const cl_R_fdiv_t fceiling2 (const cl_R& x);
extern const cl_F fceiling (const cl_R& x);

// ftruncate2(x) liefert (ftruncate x), wo x eine reelle Zahl ist.
extern const cl_R_fdiv_t ftruncate2 (const cl_R& x);
extern const cl_F ftruncate (const cl_R& x);

// fround2(x) liefert (fround x), wo x eine reelle Zahl ist.
extern const cl_R_fdiv_t fround2 (const cl_R& x);
extern const cl_F fround (const cl_R& x);

// ffloor2(x,y) liefert (ffloor x y), wo x und y reelle Zahlen sind.
extern const cl_R_fdiv_t ffloor2 (const cl_R& x, const cl_R& y);
extern const cl_F ffloor (const cl_R& x, const cl_R& y);

// fceiling2(x,y) liefert (fceiling x y), wo x und y reelle Zahlen sind.
extern const cl_R_fdiv_t fceiling2 (const cl_R& x, const cl_R& y);
extern const cl_F fceiling (const cl_R& x, const cl_R& y);

// ftruncate2(x,y) liefert (ftruncate x y), wo x und y reelle Zahlen sind.
extern const cl_R_fdiv_t ftruncate2 (const cl_R& x, const cl_R& y);
extern const cl_F ftruncate (const cl_R& x, const cl_R& y);

// fround2(x,y) liefert (fround x y), wo x und y reelle Zahlen sind.
extern const cl_R_fdiv_t fround2 (const cl_R& x, const cl_R& y);
extern const cl_F fround (const cl_R& x, const cl_R& y);


// mod(x,y) = (mod x y), wo x und y reelle Zahlen sind.
extern const cl_R mod (const cl_R& x, const cl_R& y);

// rem(x,y) = (rem x y), wo x und y reelle Zahlen sind.
extern const cl_R rem (const cl_R& x, const cl_R& y);


// rational(x) liefert (rational x), wo x eine reelle Zahl ist.
extern const cl_RA rational (const cl_R& x);
// Spezialfall:
inline const cl_RA rational (const cl_RA& x) { return x; }


// equal(x,y) vergleicht zwei reelle Zahlen x und y auf Gleichheit.
extern bool equal (const cl_R& x, const cl_R& y);
// equal_hashcode(x) liefert einen equal-invarianten Hashcode für x.
extern uint32 equal_hashcode (const cl_R& x);

// compare(x,y) vergleicht zwei reelle Zahlen x und y.
// Ergebnis: 0 falls x=y, +1 falls x>y, -1 falls x<y.
extern cl_signean compare (const cl_R& x, const cl_R& y);

inline bool operator== (const cl_R& x, const cl_R& y)
	{ return equal(x,y); }
inline bool operator!= (const cl_R& x, const cl_R& y)
	{ return !equal(x,y); }
inline bool operator<= (const cl_R& x, const cl_R& y)
	{ return compare(x,y)<=0; }
inline bool operator< (const cl_R& x, const cl_R& y)
	{ return compare(x,y)<0; }
inline bool operator>= (const cl_R& x, const cl_R& y)
	{ return compare(x,y)>=0; }
inline bool operator> (const cl_R& x, const cl_R& y)
	{ return compare(x,y)>0; }

// max(x,y) liefert (max x y), wo x und y reelle Zahlen sind.
extern const cl_R max (const cl_R& x, const cl_R& y);

// min(x,y) liefert (min x y), wo x und y reelle Zahlen sind.
extern const cl_R min (const cl_R& x, const cl_R& y);

// signum(x) liefert (signum x), wo x eine reelle Zahl ist.
extern const cl_R signum (const cl_R& x);

// sqrt(x) = (sqrt x) zieht die Wurzel aus einer reellen Zahl x >=0.
extern const cl_R sqrt (const cl_R& x);
// sqrt(x) = (sqrt x) zieht die Wurzel aus einer rationalen Zahl x >=0.
extern const cl_R sqrt (const cl_RA& x);

// (expt x y), wo x eine reelle Zahl und y ein Integer ist.
extern const cl_R expt (const cl_R& x, sintL y);
extern const cl_R expt (const cl_R& x, const cl_I& y);

// rationalize(x) liefert (rationalize x), wo x eine reelle Zahl ist.
extern const cl_RA rationalize (const cl_R& x);


// Konversion zu einem C "float".
extern float float_approx (const cl_R& x);

// Konversion zu einem C "double".
extern double double_approx (const cl_R& x);


// Transcendental functions


// atan(x,y) liefert zu zwei reellen Zahlen x, y den Winkel von (x,y)
// in Polarkoordinaten. Ergebnis rational nur, wenn x>0 und y=0.
extern const cl_R atan (const cl_R& x, const cl_R& y);
// Spezialfall: y Float -> Ergebnis Float
inline const cl_F atan (const cl_R& x, const cl_F& y)
	{ return The(cl_F)(atan(x,The(cl_R)(y))); }
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R atan (const cl_R& x, const int y)
	{ return atan(x,cl_I(y)); }
inline const cl_R atan (const cl_R& x, const unsigned int y)
	{ return atan(x,cl_I(y)); }
inline const cl_R atan (const cl_R& x, const long y)
	{ return atan(x,cl_I(y)); }
inline const cl_R atan (const cl_R& x, const unsigned long y)
	{ return atan(x,cl_I(y)); }

// atan(x) liefert den Arctan einer reellen Zahl x.
// Ergebnis rational nur, wenn x=0.
extern const cl_R atan (const cl_R& x);
// Spezialfall: x Float -> Ergebnis Float
inline const cl_F atan (const cl_F& x) { return The(cl_F)(atan(The(cl_R)(x))); }
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R atan (const int x) { return atan(cl_I(x)); }
inline const cl_R atan (const unsigned int x) { return atan(cl_I(x)); }
inline const cl_R atan (const long x) { return atan(cl_I(x)); }
inline const cl_R atan (const unsigned long x) { return atan(cl_I(x)); }

// sin(x) liefert den Sinus (sin x) einer reellen Zahl x.
extern const cl_R sin (const cl_R& x);
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R sin (const int x) { return sin(cl_I(x)); }
inline const cl_R sin (const unsigned int x) { return sin(cl_I(x)); }
inline const cl_R sin (const long x) { return sin(cl_I(x)); }
inline const cl_R sin (const unsigned long x) { return sin(cl_I(x)); }

// cos(x) liefert den Cosinus (cos x) einer reellen Zahl x.
extern const cl_R cos (const cl_R& x);
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R cos (const int x) { return cos(cl_I(x)); }
inline const cl_R cos (const unsigned int x) { return cos(cl_I(x)); }
inline const cl_R cos (const long x) { return cos(cl_I(x)); }
inline const cl_R cos (const unsigned long x) { return cos(cl_I(x)); }

// cos_sin(x) liefert ((cos x),(sin x)), beide Werte.
extern const cos_sin_t cos_sin (const cl_R& x);

// tan(x) liefert den Tangens (tan x) einer reellen Zahl x.
extern const cl_R tan (const cl_R& x);
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R tan (const int x) { return tan(cl_I(x)); }
inline const cl_R tan (const unsigned int x) { return tan(cl_I(x)); }
inline const cl_R tan (const long x) { return tan(cl_I(x)); }
inline const cl_R tan (const unsigned long x) { return tan(cl_I(x)); }

// ln(x) liefert zu einer reellen Zahl x>0 die Zahl ln(x).
extern const cl_R ln (const cl_R& x);
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R ln (const int x) { return ln(cl_I(x)); }
inline const cl_R ln (const unsigned int x) { return ln(cl_I(x)); }
inline const cl_R ln (const long x) { return ln(cl_I(x)); }
inline const cl_R ln (const unsigned long x) { return ln(cl_I(x)); }

// log(a,b) liefert zu reellen Zahlen a>0, b>0 die Zahl
// log(a,b)=ln(a)/ln(b).
// Ergebnis rational nur, wenn a=1 oder a und b rational.
extern const cl_R log (const cl_R& a, const cl_R& b);

// exp(x) liefert zu einer reellen Zahl x die Zahl exp(x).
extern const cl_R exp (const cl_R& x);
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R exp (const int x) { return exp(cl_I(x)); }
inline const cl_R exp (const unsigned int x) { return exp(cl_I(x)); }
inline const cl_R exp (const long x) { return exp(cl_I(x)); }
inline const cl_R exp (const unsigned long x) { return exp(cl_I(x)); }

// sinh(x) liefert zu einer reellen Zahl x die Zahl sinh(x).
extern const cl_R sinh (const cl_R& x);
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R sinh (const int x) { return sinh(cl_I(x)); }
inline const cl_R sinh (const unsigned int x) { return sinh(cl_I(x)); }
inline const cl_R sinh (const long x) { return sinh(cl_I(x)); }
inline const cl_R sinh (const unsigned long x) { return sinh(cl_I(x)); }

// cosh(x) liefert zu einer reellen Zahl x die Zahl cosh(x).
extern const cl_R cosh (const cl_R& x);
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R cosh (const int x) { return cosh(cl_I(x)); }
inline const cl_R cosh (const unsigned int x) { return cosh(cl_I(x)); }
inline const cl_R cosh (const long x) { return cosh(cl_I(x)); }
inline const cl_R cosh (const unsigned long x) { return cosh(cl_I(x)); }

// cosh_sinh(x) liefert ((cosh x),(sinh x)), beide Werte.
extern const cosh_sinh_t cosh_sinh (const cl_R& x);

// tanh(x) liefert zu einer reellen Zahl x die Zahl tanh(x).
extern const cl_R tanh (const cl_R& x);
// Dem C++-Compiler muß man nun auch das Folgende sagen:
inline const cl_R tanh (const int x) { return tanh(cl_I(x)); }
inline const cl_R tanh (const unsigned int x) { return tanh(cl_I(x)); }
inline const cl_R tanh (const long x) { return tanh(cl_I(x)); }
inline const cl_R tanh (const unsigned long x) { return tanh(cl_I(x)); }


// random_R(randomstate,n) liefert zu einer reellen Zahl n>0 eine Zufallszahl
// x mit 0 <= x < n.
extern const cl_R random_R (random_state& randomstate, const cl_R& n);

inline const cl_R random_R (const cl_R& n)
	{ return random_R(default_random_state,n); }


// This could be optimized to use in-place operations.
inline cl_R& operator+= (cl_R& x, const cl_R& y) { return x = x + y; }
inline cl_F& operator+= (cl_F& x, const cl_R& y) { return x = x + y; }
inline cl_F& operator+= (cl_F& x, const cl_RA& y) { return x = x + y; }
inline cl_F& operator+= (cl_F& x, const cl_I& y) { return x = x + y; }
inline cl_R& operator+= (cl_R& x, const int y) { return x = x + y; }
inline cl_R& operator+= (cl_R& x, const unsigned int y) { return x = x + y; }
inline cl_R& operator+= (cl_R& x, const long y) { return x = x + y; }
inline cl_R& operator+= (cl_R& x, const unsigned long y) { return x = x + y; }
#ifdef HAVE_LONGLONG
inline cl_R& operator+= (cl_R& x, const long long y) { return x = x + y; }
inline cl_R& operator+= (cl_R& x, const unsigned long long y) { return x = x + y; }
#endif
inline cl_F& operator+= (cl_R& x, const float y) { return static_cast<cl_F&>(x = x + y); }
inline cl_F& operator+= (cl_R& x, const double y) { return static_cast<cl_F&>(x = x + y); }
inline cl_F& operator+= (cl_F& x, const int y) { return x = x + y; }
inline cl_F& operator+= (cl_F& x, const unsigned int y) { return x = x + y; }
inline cl_F& operator+= (cl_F& x, const long y) { return x = x + y; }
inline cl_F& operator+= (cl_F& x, const unsigned long y) { return x = x + y; }
#ifdef HAVE_LONGLONG
inline cl_F& operator+= (cl_F& x, const long long y) { return x = x + y; }
inline cl_F& operator+= (cl_F& x, const unsigned long long y) { return x = x + y; }
#endif
inline cl_R& operator++ /* prefix */ (cl_R& x) { return x = plus1(x); }
inline void operator++ /* postfix */ (cl_R& x, int dummy) { (void)dummy; x = plus1(x); }
inline cl_R& operator-= (cl_R& x, const cl_R& y) { return x = x - y; }
inline cl_F& operator-= (cl_F& x, const cl_R& y) { return x = x - y; }
inline cl_F& operator-= (cl_F& x, const cl_RA& y) { return x = x - y; }
inline cl_F& operator-= (cl_F& x, const cl_I& y) { return x = x - y; }
inline cl_R& operator-= (cl_R& x, const int y) { return x = x - y; }
inline cl_R& operator-= (cl_R& x, const unsigned int y) { return x = x - y; }
inline cl_R& operator-= (cl_R& x, const long y) { return x = x - y; }
inline cl_R& operator-= (cl_R& x, const unsigned long y) { return x = x - y; }
#ifdef HAVE_LONGLONG
inline cl_R& operator-= (cl_R& x, const long long y) { return x = x - y; }
inline cl_R& operator-= (cl_R& x, const unsigned long long y) { return x = x - y; }
#endif
inline cl_F& operator-= (cl_R& x, const float y) { return static_cast<cl_F&>(x = x - y); }
inline cl_F& operator-= (cl_R& x, const double y) { return static_cast<cl_F&>(x = x - y); }
inline cl_F& operator-= (cl_F& x, const int y) { return x = x - y; }
inline cl_F& operator-= (cl_F& x, const unsigned int y) { return x = x - y; }
inline cl_F& operator-= (cl_F& x, const long y) { return x = x - y; }
inline cl_F& operator-= (cl_F& x, const unsigned long y) { return x = x - y; }
#ifdef HAVE_LONGLONG
inline cl_F& operator-= (cl_F& x, const long long y) { return x = x - y; }
inline cl_F& operator-= (cl_F& x, const unsigned long long y) { return x = x - y; }
#endif
inline cl_R& operator-- /* prefix */ (cl_R& x) { return x = minus1(x); }
inline void operator-- /* postfix */ (cl_R& x, int dummy) { (void)dummy; x = minus1(x); }
inline cl_R& operator*= (cl_R& x, const cl_R& y) { return x = x * y; }
inline cl_R& operator*= (cl_R& x, const int y) { return x = x * y; }
inline cl_R& operator*= (cl_R& x, const unsigned int y) { return x = x * y; }
inline cl_R& operator*= (cl_R& x, const long y) { return x = x * y; }
inline cl_R& operator*= (cl_R& x, const unsigned long y) { return x = x * y; }
#ifdef HAVE_LONGLONG
inline cl_R& operator*= (cl_R& x, const long long y) { return x = x * y; }
inline cl_R& operator*= (cl_R& x, const unsigned long long y) { return x = x * y; }
#endif
inline cl_R& operator*= (cl_R& x, const float y) { return x = x * y; }
inline cl_R& operator*= (cl_R& x, const double y) { return x = x * y; }
inline cl_R& operator/= (cl_R& x, const cl_R& y) { return x = x / y; }
inline cl_F& operator/= (cl_F& x, const cl_R& y) { return x = x / y; }
inline cl_F& operator/= (cl_F& x, const cl_RA& y) { return x = x / y; }
inline cl_F& operator/= (cl_F& x, const cl_I& y) { return x = x / y; }
inline cl_R& operator/= (cl_R& x, const int y) { return x = x / y; }
inline cl_R& operator/= (cl_R& x, const unsigned int y) { return x = x / y; }
inline cl_R& operator/= (cl_R& x, const long y) { return x = x / y; }
inline cl_R& operator/= (cl_R& x, const unsigned long y) { return x = x / y; }
#ifdef HAVE_LONGLONG
inline cl_R& operator/= (cl_R& x, const long long y) { return x = x / y; }
inline cl_R& operator/= (cl_R& x, const unsigned long long y) { return x = x / y; }
#endif
inline cl_R& operator/= (cl_R& x, const float y) { return x = x / y; }
inline cl_R& operator/= (cl_R& x, const double y) { return x = x / y; }
inline cl_F& operator/= (cl_F& x, const int y) { return x = x / y; }
inline cl_F& operator/= (cl_F& x, const unsigned int y) { return x = x / y; }
inline cl_F& operator/= (cl_F& x, const long y) { return x = x / y; }
inline cl_F& operator/= (cl_F& x, const unsigned long y) { return x = x / y; }
#ifdef HAVE_LONGLONG
inline cl_F& operator/= (cl_F& x, const long long y) { return x = x / y; }
inline cl_F& operator/= (cl_F& x, const unsigned long long y) { return x = x / y; }
#endif


// Complex operations, trivial for reals

inline const cl_R realpart (const cl_R& x)
{
	return x;
}
inline const cl_R imagpart (const cl_R& x)
{
	(void)x; // unused x
	return 0;
}
inline const cl_R conjugate (const cl_R& x)
{
	return x;
}


// Debugging support.
#ifdef CL_DEBUG
extern int cl_R_debug_module;
CL_FORCE_LINK(cl_R_debug_dummy, cl_R_debug_module)
#endif

}  // namespace cln

#endif /* _CL_REAL_H */
