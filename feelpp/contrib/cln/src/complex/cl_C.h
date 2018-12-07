// cl_C internals

#ifndef _CL_C_H
#define _CL_C_H

#include "cln/number.h"
#include "cln/complex.h"
#include "cln/sfloat_class.h"
#include "cln/ffloat_class.h"
#include "cln/dfloat_class.h"
#include "cln/lfloat_class.h"
#include "base/cl_macros.h"
#include "cln/malloc.h"

namespace cln {

struct cl_heap_complex : cl_heap {
	cl_R realpart;
	cl_R imagpart;
};

inline cl_heap_complex* TheComplex (cl_heap_complex* p)
	{ return p; }
inline cl_heap_complex* TheComplex (const cl_number& obj)
	{ return (cl_heap_complex*)(obj.pointer); }

inline cl_heap_complex* allocate_complex (const cl_R& real, const cl_R& imag)
{
	cl_heap_complex* p = (cl_heap_complex*) malloc_hook(sizeof(cl_heap_complex));
	p->refcount = 1;
	p->type = &cl_class_complex;
	p->realpart.pointer = real.pointer;	cl_inc_refcount(real);
	p->imagpart.pointer = imag.pointer;	cl_inc_refcount(imag);
	return p;
}

// Private constructor.
// ptr should be the result of some allocate_complex() call.
inline cl_N::cl_N (cl_heap_complex* ptr)
	: cl_number ((cl_private_thing) ptr) {}

// Both work, but the first definition results in less compiler-generated
// temporaries.
#if 1
  #define Complex  cl_heap_complex*
#else
  #define Complex  cl_N
#endif

// Type tests
inline bool realp (const cl_N& x)
{
	if (x.pointer_p())
		if (x.heappointer->type == &cl_class_complex)
			return false;
	return true;
}
inline bool complexp (const cl_N& x)
{
	if (x.pointer_p())
		if (x.heappointer->type == &cl_class_complex)
			return true;
	return false;
}

// Comparison with a fixnum.
inline bool eq (const cl_N& x, sint32 y)
{
	return x.word == cl_combine(cl_FN_tag,y);
}

inline bool exact_zerop (const cl_N& x)
{
	return eq(x,0);
}


// A complex (cl_C) is a number which is not a real number (cl_R).

// typedef
class cl_C : public cl_N {
public:
};

inline bool realp (const cl_C& x)
	{ unused x; return false; }
inline bool complexp (const cl_C& x)
	{ unused x; return true; }


// Liefert zu reellen Zahlen a und b /= Fixnum 0 die komplexe Zahl a+bi.
// complex_C(a,b)
extern const cl_N complex_C (const cl_R& a, const cl_R& b);

// realpart(x) liefert den Realteil der Zahl x.
// imagpart(x) liefert den Imaginärteil der Zahl x.
inline const cl_R& realpart (const cl_C& x)
{
	return TheComplex(x)->realpart;
}
inline const cl_R& imagpart (const cl_C& x)
{
	return TheComplex(x)->imagpart;
}


// Primitive forms of complex numbers with restricted part type.

// typedef
struct cl_C_SF {
	cl_SF realpart;
	cl_SF imagpart;
	cl_C_SF (const cl_SF& re, const cl_SF& im) : realpart(re), imagpart(im) {}
};
inline const cl_N complex_C (const cl_C_SF& c)
	{ return complex_C(c.realpart,c.imagpart); }

// typedef
struct cl_C_FF {
	cl_FF realpart;
	cl_FF imagpart;
	cl_C_FF (const cl_FF& re, const cl_FF& im) : realpart(re), imagpart(im) {}
};
inline const cl_N complex_C (const cl_C_FF& c)
	{ return complex_C(c.realpart,c.imagpart); }

// typedef
struct cl_C_DF {
	cl_DF realpart;
	cl_DF imagpart;
	cl_C_DF (const cl_DF& re, const cl_DF& im) : realpart(re), imagpart(im) {}
};
inline const cl_N complex_C (const cl_C_DF& c)
	{ return complex_C(c.realpart,c.imagpart); }

// typedef
struct cl_C_LF {
	cl_LF realpart;
	cl_LF imagpart;
	cl_C_LF (const cl_LF& re, const cl_LF& im) : realpart(re), imagpart(im) {}
};
inline const cl_N complex_C (const cl_C_LF& c)
	{ return complex_C(c.realpart,c.imagpart); }

// cl_C_recip(a,b) liefert 1/(a+bi), wo a und b Short-Floats sind.
extern const cl_C_SF cl_C_recip (const cl_SF& a, const cl_SF& b);

// cl_C_recip(a,b) liefert 1/(a+bi), wo a und b Single-Floats sind.
extern const cl_C_FF cl_C_recip (const cl_FF& a, const cl_FF& b);

// cl_C_recip(a,b) liefert 1/(a+bi), wo a und b Double-Floats sind.
extern const cl_C_DF cl_C_recip (const cl_DF& a, const cl_DF& b);

// cl_C_recip(a,b) liefert 1/(a+bi), wo a und b gleichlange Long-Floats sind.
extern const cl_C_LF cl_C_recip (const cl_LF& a, const cl_LF& b);


// cl_C_hypot(a,b) liefert abs(a+bi), wo a und b Short-Floats sind.
extern const cl_SF cl_hypot (const cl_SF& a, const cl_SF& b);

// cl_C_hypot(a,b) liefert abs(a+bi), wo a und b Single-Floats sind.
extern const cl_FF cl_hypot (const cl_FF& a, const cl_FF& b);

// cl_C_hypot(a,b) liefert abs(a+bi), wo a und b Double-Floats sind.
extern const cl_DF cl_hypot (const cl_DF& a, const cl_DF& b);

// cl_C_hypot(a,b) liefert abs(a+bi), wo a und b gleichlange Long-Floats sind.
extern const cl_LF cl_hypot (const cl_LF& a, const cl_LF& b);

// cl_C_hypot(a,b) liefert abs(a+bi), wo a und b reelle Zahlen sind.
extern const cl_R cl_hypot (const cl_R& a, const cl_R& b);

// Liefert (abs x), wo x eine nicht-reelle Zahl ist.
extern const cl_R abs (const cl_C& x);


// typedef
struct cl_C_R {
	cl_R realpart;
	cl_R imagpart;
	cl_C_R () : realpart(0), imagpart(0) {}
	cl_C_R (const cl_R& re, const cl_R& im) : realpart(re), imagpart(im) {}
};

// Hilfsfunktion für atanh und atan: u+iv := artanh(x+iy). Liefert cl_C_R(u,v).
extern const cl_C_R atanh (const cl_R& x, const cl_R& y);

// Hilfsfunktion für asinh und asin: u+iv := arsinh(x+iy). Liefert cl_C_R(u,v).
extern const cl_C_R asinh (const cl_R& x, const cl_R& y);

}  // namespace cln

#endif /* _CL_C_H */
