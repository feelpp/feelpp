// Simple vectors.

#ifndef _CL_SV_H
#define _CL_SV_H

#include "cln/object.h"
#include "cln/V.h"
#include "cln/exception.h"
#include <cstdlib>
#include <cstddef>

namespace cln {

// A simple vector has the same operations as a vector, but it can store
// _only_ cl_gcobject's.
// This class is here because the general vectors always need a function
// call for getting/setting the element of a vector. Our main application
// of the general vectors are the bit vectors, needed for implementing
// polynomials over modular integer rings. I don't want that polynomials
// over other rings (in particular cl_I) be penalized by the mere existence
// of polynomials over modular integer rings.

// When the vectors were implemented like this:
//
//    cl_GV<cl_I>  -->  cl_GV<cl_RA>  -->  cl_GV<cl_R>  -->  cl_GV<cl_N>
//
// a bit/byte-vector (of integers with limited range) could actually be
// treated correctly by all the functions which manipulate vectors of cl_N.
// This is not crucial, however. Here, we'll have disjoint sets
//
//    cl_SV<cl_I>  -->  cl_SV<cl_RA>  -->  cl_SV<cl_R>  -->  cl_SV<cl_N>
//
//    cl_GV<cl_I>
//
// i.e. the functions which manipulate a (simple!) vector of cl_N cannot
// deal with a bit/byte-vector.
// (This is the same issue as UPGRADED-ARRAY-ELEMENT-TYPE in Common Lisp.)

template <class T> class cl_SV_inner;

template <class T>
class cl_SV_inner {
protected:
	std::size_t len; // number of elements
private:
//	T data[]; // the elements
	T * data() { return (T *) (this+1); }
	const T * data() const { return (const T *) (this+1); }
public:
	std::size_t size() const { return len; } // number of elements
	const T & operator[] (unsigned long index) const
	{
		#ifndef CL_SV_NO_RANGECHECKS
		if (!(index < size())) throw runtime_exception();
		#endif
		return data()[index];
	}
	T & operator[] (unsigned long index)
	{
		#ifndef CL_SV_NO_RANGECHECKS
		if (!(index < size())) throw runtime_exception();
		#endif
		return data()[index];
	}
	// New ANSI C++ compilers also want the following.
	const T & operator[] (unsigned int index) const
	{ return operator[]((unsigned long)index); }
	T & operator[] (unsigned int index)
	{ return operator[]((unsigned long)index); }
	const T & operator[] (long index) const
	{ return operator[]((unsigned long)index); }
	T & operator[] (long index)
	{ return operator[]((unsigned long)index); }
	const T & operator[] (int index) const
	{ return operator[]((unsigned long)index); }
	T & operator[] (int index)
	{ return operator[]((unsigned long)index); }
	#if long_bitsize < pointer_bitsize
	const T & operator[] (unsigned long long index) const
	{
		#ifndef CL_SV_NO_RANGECHECKS
		if (!(index < size())) throw runtime_exception();
		#endif
		return data()[index];
	}
	T & operator[] (unsigned long long index)
	{
		#ifndef CL_SV_NO_RANGECHECKS
		if (!(index < size())) throw runtime_exception();
		#endif
		return data()[index];
	}
	const T & operator[] (long long index) const
	{ return operator[]((unsigned long long)index); }
	T & operator[] (long long index)
	{ return operator[]((unsigned long long)index); }
	#endif
public: /* ugh */
	// Constructor.
	cl_SV_inner (std::size_t l) : len (l) {}
public:
	// Destructor.
	~cl_SV_inner ();
	// Ability to place an object at a given address.
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
private:
// No default constructor, copy constructor, assignment operator, new.
	cl_SV_inner ();
	cl_SV_inner (const cl_SV_inner&);
	cl_SV_inner& operator= (const cl_SV_inner&);
	void* operator new (size_t size);
};

// All member functions are inline.

template <class T>
inline cl_SV_inner<T>::~cl_SV_inner ()
{
	std::size_t i = len;
	while (i > 0) {
		i--;
		data()[i].~T();
	}
}


// In memory, a simple vector looks like this:

template <class T>
struct cl_heap_SV : cl_heap {
	cl_SV_inner<T> v;
	// here room for the elements
};

template <class T, class BASE>
struct cl_SV : public BASE {
public:
	// Length.
	std::size_t size() const
	{
		return ((const cl_heap_SV<T> *) this->pointer)->v.size();
	}
	// Reference. Forbid modification of `const cl_SV&' arguments.
	const T & operator[] (unsigned long index) const
	{
		return ((const cl_heap_SV<T> *) this->pointer)->v[index];
	}
	T & operator[] (unsigned long index)
	{
		return ((cl_heap_SV<T> *) this->pointer)->v[index];
	}
	// New ANSI C++ compilers also want the following.
	const T & operator[] (unsigned int index) const
	{ return operator[]((unsigned long)index); }
	T & operator[] (unsigned int index)
	{ return operator[]((unsigned long)index); }
	const T & operator[] (long index) const
	{ return operator[]((unsigned long)index); }
	T & operator[] (long index)
	{ return operator[]((unsigned long)index); }
	const T & operator[] (int index) const
	{ return operator[]((unsigned long)index); }
	T & operator[] (int index)
	{ return operator[]((unsigned long)index); }
	#if long_bitsize < pointer_bitsize
	const T & operator[] (unsigned long long index) const
	{
		return ((const cl_heap_SV<T> *) this->pointer)->v[index];
	}
	T & operator[] (unsigned long long index)
	{
		return ((cl_heap_SV<T> *) this->pointer)->v[index];
	}
	const T & operator[] (long long index) const
	{ return operator[]((unsigned long long)index); }
	T & operator[] (long long index)
	{ return operator[]((unsigned long long)index); }
	#endif
	// Constructors.
	cl_SV (const cl_SV&);
	// Assignment operators.
	cl_SV& operator= (const cl_SV&);
	// Private pointer manipulations.
	cl_SV (cl_heap_SV<T>* p) : BASE ((cl_private_thing)p) {}
	cl_SV (cl_private_thing p) : BASE (p) {}
protected:
	// Forbid use of default constructor.
	cl_SV ();
};
#define CL_SV(T,BASE) cl_SV<T,BASE>
// Define copy constructor.
template <class T, class BASE>
	_CL_DEFINE_COPY_CONSTRUCTOR2(CL_SV(T,BASE),cl_SV,BASE)
// Define assignment operator.
template <class T, class BASE>
	CL_DEFINE_ASSIGNMENT_OPERATOR(CL_SV(T,BASE),CL_SV(T,BASE))
#undef CL_SV

// The "generic" simple vector type.

typedef cl_heap_SV<cl_gcobject> cl_heap_SV_any;
typedef cl_SV<cl_gcobject,cl_V_any> cl_SV_any;

// Copy a simple vector.
extern const cl_SV_any copy (const cl_SV_any&);


// Hack section.

// Conversions to subtypes without checking:
  #define The(type)  *(const type *) & cl_identity
// This inline function is for type checking purposes only.
  inline const cl_SV_any& cl_identity (const cl_SV_any& x) { return x; }

}  // namespace cln

#endif /* _CL_SV_H */
