// cl_RA internals

#ifndef _CL_RA_H
#define _CL_RA_H

#include "cln/number.h"
#include "cln/rational.h"
#include "base/cl_macros.h"
#include "cln/malloc.h"
#include "integer/cl_I.h"

namespace cln {

struct cl_heap_ratio : cl_heap {
	cl_I numerator;
	cl_I denominator;
};

inline cl_heap_ratio* TheRatio (cl_heap_ratio* p)
	{ return p; }
inline cl_heap_ratio* TheRatio (const cl_number& obj)
	{ return (cl_heap_ratio*)(obj.pointer); }

inline cl_heap_ratio* allocate_ratio (const cl_I& num, const cl_I& den)
{
	cl_heap_ratio* p = (cl_heap_ratio*) malloc_hook(sizeof(cl_heap_ratio));
	p->refcount = 1;
	p->type = &cl_class_ratio;
	p->numerator.pointer = num.pointer;	cl_inc_refcount(num);
	p->denominator.pointer = den.pointer;	cl_inc_refcount(den);
	return p;
}

// Private constructor.
// ptr should be the result of some allocate_ratio() call.
inline cl_RA::cl_RA (cl_heap_ratio* ptr)
	: cl_R ((cl_private_thing) ptr) {}

// Both work, but the first definition results in less compiler-generated
// temporaries.
#if 1
  #define Ratio  cl_heap_ratio*
#else
  #define Ratio  cl_RA
#endif

// Type tests.
inline bool rationalp (const cl_RA& x)
	{ unused x; return true; }
inline bool integerp (const cl_RA& x)
{
	if (!x.pointer_p())
		return true;
	else
		if (x.pointer_type() == &cl_class_bignum)
			return true;
	return false;
}
inline bool ratiop (const cl_RA& x)
{
	if (!x.pointer_p())
		return false;
	else
		if (x.pointer_type() == &cl_class_bignum)
			return false;
	return true;
}


// A ratio (cl_RT) is a rational number which is not an integer (cl_I).

// typedef
class cl_RT : public cl_RA {
public:
};

inline bool integerp (const cl_RT& x)
	{ unused x; return false; }
inline bool ratiop (const cl_RT& x)
	{ unused x; return true; }

// Access numerator and denominator.
inline const cl_I& numerator (const cl_RT& x)
	{ return TheRatio(x)->numerator; }
inline const cl_I& denominator (const cl_RT& x)
	{ return TheRatio(x)->denominator; }


// Sign test:

// (MINUSP x) == (< x 0)
inline bool minusp (const cl_RT& x)
	{ return minusp(numerator(x)); }
inline bool minusp (const cl_RA& x)
{
	if (ratiop(x)) {
		DeclareType(cl_RT,x);
		return minusp(x);
	} else {
		DeclareType(cl_I,x);
		return minusp(x);
	}
}

// (ZEROP x) == (= x 0)
inline bool zerop (const cl_RT& x)
	{ unused x; return false; }
inline bool zerop (const cl_RA& x)
{
	return x.word == cl_combine(cl_FN_tag,0);
}

// (EQ x y) == (= x y), assuming y a fixnum
inline bool eq (const cl_RA& x, sint32 y)
{
	return x.word == cl_combine(cl_FN_tag,y);
}

// Liefert zu den Integers a und b mit b>1 und ggT(a,b)=1 den Bruch a/b.
// I_I_to_RT(a,b)
  extern const cl_RA I_I_to_RT (const cl_I& a, const cl_I& b);

// Liefert zu den Integers a und b mit b>0 und ggT(a,b)=1 den Bruch a/b
// (Ratio oder Integer).
// I_I_to_RA(a,b)
  extern const cl_RA I_I_to_RA (const cl_I& a, const cl_I& b);

// Liefert zu den Integers a und b mit b>0 den Bruch a/b (Ratio oder Integer).
// I_posI_div_RA(a,b)
  extern const cl_RA I_posI_div_RA (const cl_I& a, const cl_I& b);

// Liefert zu den Integers a und b den Bruch a/b (Ratio oder Integer).
// I_I_div_RA(a,b)
  extern const cl_RA I_I_div_RA (const cl_I& a, const cl_I& b);

// Liefert den Zähler einer rationalen Zahl.
// numerator(r)
inline const cl_I numerator (const cl_RA& r)
{
	if (integerp(r)) {
		DeclareType(cl_I,r);
		return r;
	} else
		return TheRatio(r)->numerator;
}

// Liefert den Nenner einer rationalen Zahl.
// denominator(r)
inline const cl_I denominator (const cl_RA& r)
{
	if (integerp(r))
		return 1;
	else
		return TheRatio(r)->denominator;
}

// Liefert Zähler und Nenner einer rationalen Zahl.
// RA_numden_I_I(r, num=,den=);
// > r: rationale Zahl
// < num: (numerator r)
// < den: (denominator r)
  #define RA_numden_I_I(r,num_zuweisung,den_zuweisung)  \
    { if (integerp(r))							\
        { num_zuweisung *(const cl_I *)&r;				\
          den_zuweisung 1; /* Zähler = r, Nenner = 1 */			\
        }								\
        else								\
        { num_zuweisung TheRatio(r)->numerator;				\
          den_zuweisung TheRatio(r)->denominator;			\
        }								\
    }

}  // namespace cln

#endif /* _CL_RA_H */
