// cl_LF internals

#ifndef _CL_LF_H
#define _CL_LF_H

#include "cln/number.h"
#include "cln/lfloat_class.h"
#include "cln/integer_class.h"

namespace cln {

struct cl_heap_lfloat : cl_heap {
	uintC len;		// length of mantissa (in digits)
	int sign;		// sign (0 or -1)
	uintE expo;		// exponent
	uintD data[1];		// mantissa
};

// Minimum number of mantissa digits,
// so that a LF has not fewer mantissa bits than a DF.
  #define LF_minlen  ceiling(53,intDsize)
// Exponent.
#if (intEsize==64)
  #define LF_exp_low  1
  #define LF_exp_mid  0x8000000000000000ULL
  #define LF_exp_high 0xFFFFFFFFFFFFFFFFULL
#else
  #define LF_exp_low  1
  #define LF_exp_mid  0x80000000U
  #define LF_exp_high 0xFFFFFFFFU
#endif

inline cl_heap_lfloat* TheLfloat (cl_heap_lfloat* p)
	{ return p; }
inline cl_heap_lfloat* TheLfloat (const cl_number& obj)
	{ return (cl_heap_lfloat*)(obj.pointer); }


// Liefert zu einem Long-Float x : (futruncate x), ein LF.
// x wird von der 0 weg zur nächsten ganzen Zahl gerundet.
extern const cl_LF futruncate (const cl_LF& x);

// shorten(x,len) verkürzt ein Long-Float x auf gegebene Länge len
// und rundet dabei.
// > cl_LF x: ein Long-FLoat
// > uintC len: gewünschte Länge (>= LF_minlen, < TheLfloat(x)->len)
// < cl_LF ergebnis: verkürztes Long-Float
extern const cl_LF shorten (const cl_LF& x, uintC len);

// extend(x,len) verlängert ein Long-Float x auf gegebene Länge len.
// > cl_LF x: ein Long-FLoat
// > uintC len: gewünschte Länge (> TheLfloat(x)->len)
// < cl_LF ergebnis: verlängertes Long-Float
extern const cl_LF extend (const cl_LF& x, uintC len);

// LF_to_LF(x,len) wandelt ein Long-Float x in ein Long-Float gegebener Länge
// len um und rundet dabei nötigenfalls.
// > cl_LF x: ein Long-FLoat
// > uintC len: gewünschte Länge (>= LF_minlen)
// < cl_LF ergebnis: Long-Float gegebener Länge
extern const cl_LF LF_to_LF (const cl_LF& x, uintC len);

// GEN_LF_OP2(arg1,arg2,LF_OP,ergebnis_zuweisung)
// generates the body of a LF operation with two arguments.
// LF_OP is only executed once both arguments have been converted to the same
// float format (the longer one of arg1 and arg2). The result is then
// converted the shorter of the two float formats.
#define GEN_LF_OP2(arg1,arg2,LF_OP,ergebnis_zuweisung)  \
{									\
	var uintC len1 = TheLfloat(arg1)->len;				\
	var uintC len2 = TheLfloat(arg2)->len;				\
	if (len1==len2) /* gleich -> direkt ausführen */		\
		return LF_OP(arg1,arg2);				\
	elif (len1>len2) /* -> arg2 auf die Länge von arg1 bringen */	\
		return shorten(LF_OP(arg1,extend(arg2,len1)),len2);	\
	else /* (len1<len2) -> arg1 auf die Länge von arg2 bringen */	\
		return shorten(LF_OP(extend(arg1,len2),arg2),len1);	\
}

// Liefert zu zwei gleichlangen Long-Float x und y : (+ x y), ein LF.
// LF_LF_plus_LF(x)
extern const cl_LF LF_LF_plus_LF (const cl_LF& x, const cl_LF& y);

// Liefert zu zwei gleichlangen Long-Float x und y : (- x y), ein LF.
// LF_LF_minus_LF(x)
extern const cl_LF LF_LF_minus_LF (const cl_LF& x, const cl_LF& y);

// Use this macro if ALL of your cl_LF operations (+, -, *, /) in the
// rest of your file ALWAYS get two operands of the same precision.
#define ALL_cl_LF_OPERATIONS_SAME_PRECISION()  \
									\
inline const cl_LF operator+ (const cl_LF& x, const cl_LF& y)		\
{									\
	return LF_LF_plus_LF(x,y);					\
}									\
									\
inline const cl_LF operator- (const cl_LF& x, const cl_LF& y)		\
{									\
	return LF_LF_minus_LF(x,y);					\
}

// LF_to_I(x) wandelt ein Long-Float x, das eine ganze Zahl darstellt,
// in ein Integer um.
extern const cl_I cl_LF_to_I (const cl_LF& x);

// cl_I_to_LF(x,len) wandelt ein Integer x in ein Long-Float um und rundet dabei.
extern const cl_LF cl_I_to_LF (const cl_I& x, uintC len);

// cl_RA_to_LF(x,len) wandelt eine rationale Zahl x in ein Long-Float um
// und rundet dabei.
extern const cl_LF cl_RA_to_LF (const cl_RA& x, uintC len);

// cl_LF_I_mul(x,y) multipliziert ein Long-Float x und ein Integer y.
extern const cl_R cl_LF_I_mul (const cl_LF& x, const cl_I& y);

// cl_LF_I_div(x,y) dividiert ein Long-Float x durch ein Integer y.
extern const cl_LF cl_LF_I_div (const cl_LF& x, const cl_I& y);

// cl_I_LF_div(x,y) dividiert ein Integer x durch ein Long-Float y.
extern const cl_R cl_I_LF_div (const cl_I& x, const cl_LF& y);

// cl_LF_RA_mul(x,y) multipliziert ein Long-Float x und eine rationale Zahl y.
extern const cl_R cl_LF_RA_mul (const cl_LF& x, const cl_RA& y);

// cl_LF_RA_div(x,y) dividiert ein Long-Float x durch eine rationale Zahl y.
extern const cl_LF cl_LF_RA_div (const cl_LF& x, const cl_RA& y);

// cl_RA_LF_div(x,y) dividiert eine rationale Zahl x durch ein Long-Float y.
extern const cl_R cl_RA_LF_div (const cl_RA& x, const cl_LF& y);

// Vergrößert eine Long-Float-Länge n, so daß aus d = intDsize*n
// mindestens d+sqrt(d)+2 wird.
extern uintC cl_LF_len_incsqrt (uintC len);

// Vergrößert eine Long-Float-Länge n, so daß aus d = intDsize*n
// mindestens d+sqrt(d)+2+(LF_exp_len-1) wird.
extern uintC cl_LF_len_incsqrtx (uintC len);

// cl_LF_shortenrelative(x,y) tries to reduce the size of x, such that one
// wouldn't notice it when adding x to y. y must be /= 0. More precisely,
// this returns a float approximation of x, such that 1 ulp(x) < 1 ulp(y).
extern const cl_LF cl_LF_shortenrelative (const cl_LF& x, const cl_LF& y);

// cl_LF_shortenwith(x,y) tries to reduce the size of x, such that still
// 1 ulp(x) < y. y must be >0.
extern const cl_LF cl_LF_shortenwith (const cl_LF& x, const cl_LF& y);

}  // namespace cln

#endif /* _CL_LF_H */
