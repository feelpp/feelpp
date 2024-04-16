// General vectors.

#ifndef _CL_GV_H
#define _CL_GV_H

#include "cln/object.h"
#include "cln/V.h"
#include "cln/exception.h"
#include <cstdlib>
#include <cstddef>

namespace cln {

// A vector is a structure having the following interface:
//     v.size()        returns the number of elements
//     v[i]              returns the i-th element (0<=i<length), as a
//                       pseudo-lvalue (you can assign to it, but not take its
//                       address - exactly what you want for bit-vectors)
// This is implemented by letting v[i] be of a special "vector index" type.

template <class T> class cl_GV_inner;
template <class T> class cl_GV_index;
template <class T> class cl_GV_constindex;
template <class T> struct cl_GV_vectorops;

template <class T>
class cl_GV_inner {
protected:
	std::size_t len; // number of elements
public:
	std::size_t size() const; // number of elements
	cl_GV_vectorops<T>* vectorops; // get/set element
	const cl_GV_index<T> operator[] (unsigned long index);
	const cl_GV_constindex<T> operator[] (unsigned long index) const;
	const cl_GV_index<T> operator[] (long index);
	const cl_GV_constindex<T> operator[] (long index) const;
	const cl_GV_index<T> operator[] (unsigned int index);
	const cl_GV_constindex<T> operator[] (unsigned int index) const;
	const cl_GV_index<T> operator[] (int index);
	const cl_GV_constindex<T> operator[] (int index) const;
	#if long_bitsize < pointer_bitsize
	const cl_GV_index<T> operator[] (unsigned long long index);
	const cl_GV_constindex<T> operator[] (unsigned long long index) const;
	const cl_GV_index<T> operator[] (long long index);
	const cl_GV_constindex<T> operator[] (long long index) const;
	#endif
public: /* ugh */
	// Constructor.
	cl_GV_inner (std::size_t l, cl_GV_vectorops<T>* ops) : len (l), vectorops (ops) {}
public:
	// Destructor.
	~cl_GV_inner ();
	// Ability to place an object at a given address.
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
private:
// No default constructor, copy constructor, assignment operator, new.
	cl_GV_inner ();
	cl_GV_inner (const cl_GV_inner&);
	cl_GV_inner& operator= (const cl_GV_inner&);
	void* operator new (size_t size)
		{ (void)size; return (void*)1; } // SGI CC needs this definition
// Friend declarations. They are for the compiler. Just ignore them.
	friend class cl_GV_index<T>;
	friend class cl_GV_constindex<T>;
};

template <class T>
class cl_GV_index {
	// This is the class of objects created by accessing a non-const vector
	// through [].
public:
	cl_GV_inner<T>* vec;
	std::size_t index;
	operator T () const;
	// Constructor:
	cl_GV_index (cl_GV_inner<T>* v, std::size_t i) : vec (v), index (i) {}
	// Assignment operator.
	void operator= (const T& x) const;
#if (defined(__sparc__) || defined(__sparc64__) || defined(__mips__) || defined(__mips64__)) && !defined(__GNUC__) // maybe an SGI CC and Sun CC bug?
	void operator= (const cl_GV_index<T>&) const;
	void operator= (const cl_GV_constindex<T>&) const;
#else
private:
	// No assignment operator.
	cl_GV_index& operator= (const cl_GV_index&);
#endif
private:
// No default constructor.
	cl_GV_index ();
};

template <class T>
class cl_GV_constindex {
	// This is the class of objects created by accessing a const vector
	// through []. It lacks the assignment operator.
public:
	const cl_GV_inner<T>* vec;
	std::size_t index;
	operator T () const;
	// Constructor:
	cl_GV_constindex (const cl_GV_inner<T>* v, std::size_t i) : vec (v), index (i) {}
private:
// No default constructor, assignment operator.
	cl_GV_constindex ();
	cl_GV_constindex& operator= (const cl_GV_constindex&);
};

template <class T>
struct cl_GV_vectorops {
	const T (*element) (const cl_GV_inner<T>* vec, std::size_t index);
	void (*set_element) (cl_GV_inner<T>* vec, std::size_t index, const T& x);
	void (*do_delete) (cl_GV_inner<T>* vec);
	void (*copy_elements) (const cl_GV_inner<T>* srcvec, std::size_t srcindex, cl_GV_inner<T>* destvec, std::size_t destindex, std::size_t count);
};

// All member functions are inline.

template <class T>
inline std::size_t cl_GV_inner<T>::size() const
{
	return len;
}

template <class T>
inline const cl_GV_index<T> cl_GV_inner<T>::operator[] (unsigned long index)
{
	return cl_GV_index<T>(this,index);
}

template <class T>
inline const cl_GV_constindex<T> cl_GV_inner<T>::operator[] (unsigned long index) const
{
	return cl_GV_constindex<T>(this,index);
}

template <class T>
inline const cl_GV_index<T> cl_GV_inner<T>::operator[] (long index)
{
	return operator[]((unsigned long)index);
}

template <class T>
inline const cl_GV_constindex<T> cl_GV_inner<T>::operator[] (long index) const
{
	return operator[]((unsigned long)index);
}

template <class T>
inline const cl_GV_index<T> cl_GV_inner<T>::operator[] (unsigned int index)
{
	return operator[]((unsigned long)index);
}

template <class T>
inline const cl_GV_constindex<T> cl_GV_inner<T>::operator[] (unsigned int index) const
{
	return operator[]((unsigned long)index);
}

template <class T>
inline const cl_GV_index<T> cl_GV_inner<T>::operator[] (int index)
{
	return operator[]((unsigned long)index);
}

template <class T>
inline const cl_GV_constindex<T> cl_GV_inner<T>::operator[] (int index) const
{
	return operator[]((unsigned long)index);
}

#if long_bitsize < pointer_bitsize

template <class T>
inline const cl_GV_index<T> cl_GV_inner<T>::operator[] (unsigned long long index)
{
	return cl_GV_index<T>(this,index);
}

template <class T>
inline const cl_GV_constindex<T> cl_GV_inner<T>::operator[] (unsigned long long index) const
{
	return cl_GV_constindex<T>(this,index);
}

template <class T>
inline const cl_GV_index<T> cl_GV_inner<T>::operator[] (long long index)
{
	return operator[]((unsigned long)index);
}

template <class T>
inline const cl_GV_constindex<T> cl_GV_inner<T>::operator[] (long long index) const
{
	return operator[]((unsigned long)index);
}

#endif

template <class T>
inline cl_GV_inner<T>::~cl_GV_inner ()
{
	vectorops->do_delete(this);
}

template <class T>
inline cl_GV_index<T>::operator T () const
{
	#ifndef CL_GV_NO_RANGECHECKS
	if (!(index < vec->len)) throw runtime_exception();
	#endif
	return vec->vectorops->element(vec,index);
}

template <class T>
inline void cl_GV_index<T>::operator= (const T& x) const
{
	#ifndef CL_GV_NO_RANGECHECKS
	if (!(index < vec->len)) throw runtime_exception();
	#endif
	vec->vectorops->set_element(vec,index,x);
}

template <class T>
inline cl_GV_constindex<T>::operator T () const
{
	#ifndef CL_GV_NO_RANGECHECKS
	if (!(index < vec->len)) throw runtime_exception();
	#endif
	return vec->vectorops->element(vec,index);
}

#if (defined(__sparc__) || defined(__mips__) || defined(__mips64__)) && !defined(__GNUC__) // maybe an SGI CC and Sun CC bug? handle "y[j] = x[i];"
template <class T>
inline void cl_GV_index<T>::operator= (const cl_GV_index<T>& x) const
{ operator= ((T) x); }
template <class T>
inline void cl_GV_index<T>::operator= (const cl_GV_constindex<T>& x) const
{ operator= ((T) x); }
#endif


// In memory, a vector looks like this:

template <class T>
struct cl_heap_GV : cl_heap {
	cl_GV_inner<T> v;
	// here room for the elements
};

// And a reference to a vector always looks like this:

template <class T, class BASE>
struct cl_GV : public BASE {
public:
	// Length.
	std::size_t size() const
	{
		return ((const cl_heap_GV<T> *) this->pointer)->v.size();
	}
	// Reference. Forbid modification of `const cl_GV&' arguments.
	const cl_GV_constindex<T> operator[] (unsigned long index) const
	{
		return ((const cl_heap_GV<T> *) this->pointer)->v[index];
	}
	const cl_GV_index<T> operator[] (unsigned long index)
	{
		return ((cl_heap_GV<T> *) this->pointer)->v[index];
	}
	const cl_GV_constindex<T> operator[] (long index) const
	{ return operator[]((unsigned long)index); }
	const cl_GV_index<T> operator[] (long index)
	{ return operator[]((unsigned long)index); }
	const cl_GV_constindex<T> operator[] (unsigned int index) const
	{ return operator[]((unsigned long)index); }
	const cl_GV_index<T> operator[] (unsigned int index)
	{ return operator[]((unsigned long)index); }
	const cl_GV_constindex<T> operator[] (int index) const
	{ return operator[]((unsigned long)index); }
	const cl_GV_index<T> operator[] (int index)
	{ return operator[]((unsigned long)index); }
	#if long_bitsize < pointer_bitsize
	const cl_GV_constindex<T> operator[] (unsigned long long index) const
	{
		return ((const cl_heap_GV<T> *) this->pointer)->v[index];
	}
	const cl_GV_index<T> operator[] (unsigned long long index)
	{
		return ((cl_heap_GV<T> *) this->pointer)->v[index];
	}
	const cl_GV_constindex<T> operator[] (long long index) const
	{ return operator[]((unsigned long long)index); }
	const cl_GV_index<T> operator[] (long long index)
	{ return operator[]((unsigned long long)index); }
	#endif
	// Copy constructor.
	cl_GV (const cl_GV&);
	// Assignment operator.
	cl_GV& operator= (const cl_GV&);
	// Copy a piece of a vector into another vector.
	// (Both vectors must be of the same type. Overlapping not allowed.)
	static void copy_elements (const cl_GV& src, std::size_t srcindex, cl_GV& dest, std::size_t destindex, std::size_t count)
	{
		const cl_heap_GV<T> * hsrc = (const cl_heap_GV<T> *) src.pointer;
		cl_heap_GV<T> * hdest = (cl_heap_GV<T> *) dest.pointer;
		if (!(hsrc->v.vectorops == hdest->v.vectorops))
			throw runtime_exception();
		hsrc->v.vectorops->copy_elements(&hsrc->v,srcindex,&hdest->v,destindex,count);
	}
	// Private pointer manipulations.
	operator cl_heap_GV<T>* () const;
	cl_GV (cl_heap_GV<T>* p) : BASE ((cl_private_thing) p) {}
	cl_GV (cl_private_thing p) : BASE (p) {}
protected:
	// Forbid use of default constructor.
	cl_GV ();
};
#define CL_GV(T,BASE) cl_GV<T,BASE>
// Define copy constructor.
template <class T, class BASE>
	_CL_DEFINE_COPY_CONSTRUCTOR2(CL_GV(T,BASE),cl_GV,BASE)
// Define assignment operator.
template <class T, class BASE>
	CL_DEFINE_ASSIGNMENT_OPERATOR(CL_GV(T,BASE),CL_GV(T,BASE))
// Private pointer manipulations. Never throw away a `struct cl_heap_GV<T> *'!
template <class T, class BASE>
inline CL_GV(T,BASE)::operator cl_heap_GV<T>* () const
{
	cl_heap_GV<T>* hpointer = (cl_heap_GV<T>*)this->pointer;
	cl_inc_refcount(*this);
	return hpointer;
}
#undef CL_GV

// The "generic" general vector type.

typedef cl_heap_GV<cl_gcobject> cl_heap_GV_any;
typedef cl_GV<cl_gcobject,cl_V_any> cl_GV_any;


// Hack section.

// Conversions to subtypes without checking:
  #define The(type)  *(const type *) & cl_identity
// This inline function is for type checking purposes only.
  inline const cl_GV_any& cl_identity (const cl_GV_any& x) { return x; }

}  // namespace cln

#endif /* _CL_GV_H */
