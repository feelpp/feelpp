// Public rational number operations.

#ifndef _CL_RATIONAL_H
#define _CL_RATIONAL_H

#include "cln/number.h"
#include "cln/rational_class.h"
#include "cln/integer_class.h"
#include "cln/exception.h"

namespace cln {

CL_DEFINE_AS_CONVERSION(cl_RA)


// numerator(r) liefert den Zähler der rationalen Zahl r.
extern const cl_I numerator (const cl_RA& r);

// denominator(r) liefert den Nenner (> 0) der rationalen Zahl r.
extern const cl_I denominator (const cl_RA& r);


// Liefert (- r), wo r eine rationale Zahl ist.
extern const cl_RA operator- (const cl_RA& r);

// (+ r s), wo r und s rationale Zahlen sind.
extern const cl_RA operator+ (const cl_RA& r, const cl_RA& s);
// Dem C++-Compiler muß man auch das Folgende sagen:
inline const cl_RA operator+ (const int x, const cl_RA& y)
	{ return cl_I(x) + y; }
inline const cl_RA operator+ (const unsigned int x, const cl_RA& y)
	{ return cl_I(x) + y; }
inline const cl_RA operator+ (const long x, const cl_RA& y)
	{ return cl_I(x) + y; }
inline const cl_RA operator+ (const unsigned long x, const cl_RA& y)
	{ return cl_I(x) + y; }
#ifdef HAVE_LONGLONG
inline const cl_RA operator+ (const long long x, const cl_RA& y)
	{ return cl_I(x) + y; }
inline const cl_RA operator+ (const unsigned long long x, const cl_RA& y)
	{ return cl_I(x) + y; }
#endif
inline const cl_RA operator+ (const cl_RA& x, const int y)
	{ return x + cl_I(y); }
inline const cl_RA operator+ (const cl_RA& x, const unsigned int y)
	{ return x + cl_I(y); }
inline const cl_RA operator+ (const cl_RA& x, const long y)
	{ return x + cl_I(y); }
inline const cl_RA operator+ (const cl_RA& x, const unsigned long y)
	{ return x + cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_RA operator+ (const cl_RA& x, const long long y)
	{ return x + cl_I(y); }
inline const cl_RA operator+ (const cl_RA& x, const unsigned long long y)
	{ return x + cl_I(y); }
#endif

// (- r s), wo r und s rationale Zahlen sind.
extern const cl_RA operator- (const cl_RA& r, const cl_RA& s);
// Dem C++-Compiler muß man auch das Folgende sagen:
inline const cl_RA operator- (const int x, const cl_RA& y)
	{ return cl_I(x) - y; }
inline const cl_RA operator- (const unsigned int x, const cl_RA& y)
	{ return cl_I(x) - y; }
inline const cl_RA operator- (const long x, const cl_RA& y)
	{ return cl_I(x) - y; }
inline const cl_RA operator- (const unsigned long x, const cl_RA& y)
	{ return cl_I(x) - y; }
#ifdef HAVE_LONGLONG
inline const cl_RA operator- (const long long x, const cl_RA& y)
	{ return cl_I(x) - y; }
inline const cl_RA operator- (const unsigned long long x, const cl_RA& y)
	{ return cl_I(x) - y; }
#endif
inline const cl_RA operator- (const cl_RA& x, const int y)
	{ return x - cl_I(y); }
inline const cl_RA operator- (const cl_RA& x, const unsigned int y)
	{ return x - cl_I(y); }
inline const cl_RA operator- (const cl_RA& x, const long y)
	{ return x - cl_I(y); }
inline const cl_RA operator- (const cl_RA& x, const unsigned long y)
	{ return x - cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_RA operator- (const cl_RA& x, const long long y)
	{ return x - cl_I(y); }
inline const cl_RA operator- (const cl_RA& x, const unsigned long long y)
	{ return x - cl_I(y); }
#endif

// (1+ r), wo r eine rationale Zahl ist.
extern const cl_RA plus1 (const cl_RA& r);

// (1- r), wo r eine rationale Zahl ist.
extern const cl_RA minus1 (const cl_RA& r);

// (abs r), wo r eine rationale Zahl ist.
extern const cl_RA abs (const cl_RA& r);

// equal(r,s) vergleicht zwei rationale Zahlen r und s auf Gleichheit.
extern bool equal (const cl_RA& r, const cl_RA& s);
// equal_hashcode(r) liefert einen equal-invarianten Hashcode für r.
extern uint32 equal_hashcode (const cl_RA& r);

// compare(r,s) vergleicht zwei rationale Zahlen r und s.
// Ergebnis: 0 falls r=s, +1 falls r>s, -1 falls r<s.
extern cl_signean compare (const cl_RA& r, const cl_RA& s);

inline bool operator== (const cl_RA& x, const cl_RA& y)
	{ return equal(x,y); }
inline bool operator!= (const cl_RA& x, const cl_RA& y)
	{ return !equal(x,y); }
inline bool operator<= (const cl_RA& x, const cl_RA& y)
	{ return compare(x,y)<=0; }
inline bool operator< (const cl_RA& x, const cl_RA& y)
	{ return compare(x,y)<0; }
inline bool operator>= (const cl_RA& x, const cl_RA& y)
	{ return compare(x,y)>=0; }
inline bool operator> (const cl_RA& x, const cl_RA& y)
	{ return compare(x,y)>0; }

// minusp(x) == (< x 0)
extern bool minusp (const cl_RA& x);

// zerop(x) stellt fest, ob eine rationale Zahl = 0 ist.
extern bool zerop (const cl_RA& x);

// plusp(x) == (> x 0)
extern bool plusp (const cl_RA& x);

// Kehrwert (/ r), wo r eine rationale Zahl ist.
extern const cl_RA recip (const cl_RA& r);

// Liefert (* r s), wo r und s rationale Zahlen sind.
extern const cl_RA operator* (const cl_RA& r, const cl_RA& s);
// Dem C++-Compiler muß man auch das Folgende sagen:
inline const cl_RA operator* (const int x, const cl_RA& y)
	{ return cl_I(x) * y; }
inline const cl_RA operator* (const unsigned int x, const cl_RA& y)
	{ return cl_I(x) * y; }
inline const cl_RA operator* (const long x, const cl_RA& y)
	{ return cl_I(x) * y; }
inline const cl_RA operator* (const unsigned long x, const cl_RA& y)
	{ return cl_I(x) * y; }
#ifdef HAVE_LONGLONG
inline const cl_RA operator* (const long long x, const cl_RA& y)
	{ return cl_I(x) * y; }
inline const cl_RA operator* (const unsigned long long x, const cl_RA& y)
	{ return cl_I(x) * y; }
#endif
inline const cl_RA operator* (const cl_RA& x, const int y)
	{ return x * cl_I(y); }
inline const cl_RA operator* (const cl_RA& x, const unsigned int y)
	{ return x * cl_I(y); }
inline const cl_RA operator* (const cl_RA& x, const long y)
	{ return x * cl_I(y); }
inline const cl_RA operator* (const cl_RA& x, const unsigned long y)
	{ return x * cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_RA operator* (const cl_RA& x, const long long y)
	{ return x * cl_I(y); }
inline const cl_RA operator* (const cl_RA& x, const unsigned long long y)
	{ return x * cl_I(y); }
#endif

// Quadrat (* r r), wo r eine rationale Zahl ist.
extern const cl_RA square (const cl_RA& r);

// Liefert (/ r s), wo r und s rationale Zahlen sind.
extern const cl_RA operator/ (const cl_RA& r, const cl_RA& s);
// Dem C++-Compiler muß man auch das Folgende sagen:
inline const cl_RA operator/ (const int x, const cl_RA& y)
	{ return cl_I(x) / y; }
inline const cl_RA operator/ (const unsigned int x, const cl_RA& y)
	{ return cl_I(x) / y; }
inline const cl_RA operator/ (const long x, const cl_RA& y)
	{ return cl_I(x) / y; }
inline const cl_RA operator/ (const unsigned long x, const cl_RA& y)
	{ return cl_I(x) / y; }
#ifdef HAVE_LONGLONG
inline const cl_RA operator/ (const long long x, const cl_RA& y)
	{ return cl_I(x) / y; }
inline const cl_RA operator/ (const unsigned long long x, const cl_RA& y)
	{ return cl_I(x) / y; }
#endif
inline const cl_RA operator/ (const cl_RA& x, const int y)
	{ return x / cl_I(y); }
inline const cl_RA operator/ (const cl_RA& x, const unsigned int y)
	{ return x / cl_I(y); }
inline const cl_RA operator/ (const cl_RA& x, const long y)
	{ return x / cl_I(y); }
inline const cl_RA operator/ (const cl_RA& x, const unsigned long y)
	{ return x / cl_I(y); }
#ifdef HAVE_LONGLONG
inline const cl_RA operator/ (const cl_RA& x, const long long y)
	{ return x / cl_I(y); }
inline const cl_RA operator/ (const cl_RA& x, const unsigned long long y)
	{ return x / cl_I(y); }
#endif

// Return type for rounding operators.
// x / y  --> (q,r) with x = y*q+r.
struct cl_RA_div_t {
	cl_I quotient;
	cl_RA remainder;
// Constructor.
	cl_RA_div_t () {}
	cl_RA_div_t (const cl_I& q, const cl_RA& r) : quotient(q), remainder(r) {}
};

// Liefert ganzzahligen und gebrochenen Anteil einer rationalen Zahl.
// (q,r) := (floor x)
// floor2(x)
// > x: rationale Zahl
// < q,r: Quotient q, ein Integer, Rest r, eine rationale Zahl
  extern const cl_RA_div_t floor2 (const cl_RA& x);
  extern const cl_I floor1 (const cl_RA& x);

// Liefert ganzzahligen und gebrochenen Anteil einer rationalen Zahl.
// (q,r) := (ceiling x)
// ceiling2(x)
// > x: rationale Zahl
// < q,r: Quotient q, ein Integer, Rest r, eine rationale Zahl
  extern const cl_RA_div_t ceiling2 (const cl_RA& x);
  extern const cl_I ceiling1 (const cl_RA& x);

// Liefert ganzzahligen und gebrochenen Anteil einer rationalen Zahl.
// (q,r) := (truncate x)
// truncate2(x)
// > x: rationale Zahl
// < q,r: Quotient q, ein Integer, Rest r, eine rationale Zahl
  extern const cl_RA_div_t truncate2 (const cl_RA& x);
  extern const cl_I truncate1 (const cl_RA& x);

// Liefert ganzzahligen und gebrochenen Anteil einer rationalen Zahl.
// (q,r) := (round x)
// round2(x)
// > x: rationale Zahl
// < q,r: Quotient q, ein Integer, Rest r, eine rationale Zahl
  extern const cl_RA_div_t round2 (const cl_RA& x);
  extern const cl_I round1 (const cl_RA& x);

// floor2(x,y) liefert (floor x y).
extern const cl_RA_div_t floor2 (const cl_RA& x, const cl_RA& y);
extern const cl_I floor1 (const cl_RA& x, const cl_RA& y);

// ceiling2(x,y) liefert (ceiling x y).
extern const cl_RA_div_t ceiling2 (const cl_RA& x, const cl_RA& y);
extern const cl_I ceiling1 (const cl_RA& x, const cl_RA& y);

// truncate2(x,y) liefert (truncate x y).
extern const cl_RA_div_t truncate2 (const cl_RA& x, const cl_RA& y);
extern const cl_I truncate1 (const cl_RA& x, const cl_RA& y);

// round2(x,y) liefert (round x y).
extern const cl_RA_div_t round2 (const cl_RA& x, const cl_RA& y);
extern const cl_I round1 (const cl_RA& x, const cl_RA& y);

// max(x,y) liefert (max x y), wo x und y rationale Zahlen sind.
extern const cl_RA max (const cl_RA& x, const cl_RA& y);

// min(x,y) liefert (min x y), wo x und y rationale Zahlen sind.
extern const cl_RA min (const cl_RA& x, const cl_RA& y);

// signum(x) liefert (signum x), wo x eine rationale Zahl ist.
extern const cl_RA signum (const cl_RA& x);

// (expt x y), wo x eine rationale Zahl und y ein Integer >0 ist.
extern const cl_RA expt_pos (const cl_RA& x, uintL y);
extern const cl_RA expt_pos (const cl_RA& x, const cl_I& y);

// (expt x y), wo x eine rationale Zahl und y ein Integer ist.
extern const cl_RA expt (const cl_RA& x, sintL y);
extern const cl_RA expt (const cl_RA& x, const cl_I& y);

// Stellt fest, ob eine rationale Zahl >=0 das Quadrat einer rationalen Zahl
// ist.
// sqrtp(x,&w)
// > x: eine rationale Zahl >=0
// < w: rationale Zahl (sqrt x) falls x Quadratzahl
// < ergebnis: true      ..................., false sonst
  extern bool sqrtp (const cl_RA& x, cl_RA* w);

// Stellt fest, ob eine rationale Zahl >=0 die n-te Potenz einer rationalen Zahl
// ist.
// rootp(x,n,&w)
// > x: eine rationale Zahl >=0
// > n: ein Integer >0
// < w: exakte n-te Wurzel (expt x (/ n)) falls x eine n-te Potenz
// < ergebnis: true                       ........................, false sonst
  extern bool rootp (const cl_RA& x, uintL n, cl_RA* w);
  extern bool rootp (const cl_RA& x, const cl_I& n, cl_RA* w);

// Liefert zu Integers a>0, b>1 den Logarithmus log(a,b),
// falls er eine rationale Zahl ist.
// logp(a,b,&l)
// > a: ein Integer >0
// > b: ein Integer >1
// < l: log(a,b)       falls er eine exakte rationale Zahl ist
// < ergebnis: true    ......................................., false sonst
  extern bool logp (const cl_I& a, const cl_I& b, cl_RA* l);

// Liefert zu rationalen Zahlen a>0, b>0 den Logarithmus log(a,b),
// falls er eine rationale Zahl ist.
// logp(a,b,&l)
// > a: eine rationale Zahl >0
// > b: eine rationale Zahl >0, /=1
// < l: log(a,b)       falls er eine exakte rationale Zahl ist
// < ergebnis: true    ......................................., false sonst
  extern bool logp (const cl_RA& a, const cl_RA& b, cl_RA* l);

// Konversion zu einem C "float".
extern float float_approx (const cl_RA& x);

// Konversion zu einem C "double".
extern double double_approx (const cl_RA& x);


// This could be optimized to use in-place operations.
inline cl_RA& operator+= (cl_RA& x, const cl_RA& y) { return x = x + y; }
inline cl_RA& operator+= (cl_RA& x, const int y) { return x = x + y; }
inline cl_RA& operator+= (cl_RA& x, const unsigned int y) { return x = x + y; }
inline cl_RA& operator+= (cl_RA& x, const long y) { return x = x + y; }
inline cl_RA& operator+= (cl_RA& x, const unsigned long y) { return x = x + y; }
#ifdef HAVE_LONGLONG
inline cl_RA& operator+= (cl_RA& x, const long long y) { return x = x + y; }
inline cl_RA& operator+= (cl_RA& x, const unsigned long long y) { return x = x + y; }
#endif
inline cl_RA& operator++ /* prefix */ (cl_RA& x) { return x = plus1(x); }
inline void operator++ /* postfix */ (cl_RA& x, int dummy) { (void)dummy; x = plus1(x); }
inline cl_RA& operator-= (cl_RA& x, const cl_RA& y) { return x = x - y; }
inline cl_RA& operator-= (cl_RA& x, const int y) { return x = x - y; }
inline cl_RA& operator-= (cl_RA& x, const unsigned int y) { return x = x - y; }
inline cl_RA& operator-= (cl_RA& x, const long y) { return x = x - y; }
inline cl_RA& operator-= (cl_RA& x, const unsigned long y) { return x = x - y; }
#ifdef HAVE_LONGLONG
inline cl_RA& operator-= (cl_RA& x, const long long y) { return x = x - y; }
inline cl_RA& operator-= (cl_RA& x, const unsigned long long y) { return x = x - y; }
#endif
inline cl_RA& operator-- /* prefix */ (cl_RA& x) { return x = minus1(x); }
inline void operator-- /* postfix */ (cl_RA& x, int dummy) { (void)dummy; x = minus1(x); }
inline cl_RA& operator*= (cl_RA& x, const cl_RA& y) { return x = x * y; }
inline cl_RA& operator/= (cl_RA& x, const cl_RA& y) { return x = x / y; }


// Runtime typing support.
extern cl_class cl_class_ratio;


// Debugging support.
#ifdef CL_DEBUG
extern int cl_RA_debug_module;
CL_FORCE_LINK(cl_RA_debug_dummy, cl_RA_debug_module)
#endif

}  // namespace cln

#endif /* _CL_RATIONAL_H */
