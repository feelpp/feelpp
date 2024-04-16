// cl_rcpointer hash tables

#ifndef _CL_RCPOINTER2_HASHWEAK_RCPOINTER_H
#define _CL_RCPOINTER2_HASHWEAK_RCPOINTER_H

#include "cln/object.h"

#include "base/hash/cl_hash2weak.h"

namespace cln {

// Equality.
static inline bool equal (const cl_rcpointer& x, const cl_rcpointer& y)
{ return (x.pointer == y.pointer); }

// Hash code. Luckily objects don't move around in memory.
inline uintptr_t hashcode (const cl_rcpointer& x1, const cl_rcpointer& x2)
{
	var uintptr_t hashcode1 = (uintptr_t)x1.pointer;
	var uintptr_t hashcode2 = (uintptr_t)x2.pointer;
	hashcode2 = (hashcode2 << 5) | (hashcode2 >> (long_bitsize-5)); // rotate
	return hashcode1 ^ hashcode2;
}

typedef cl_htentry2<cl_rcpointer,cl_rcpointer,cl_rcpointer> cl_htentry_from_rcpointer2_to_rcpointer;

typedef cl_heap_weak_hashtable_2<cl_rcpointer,cl_rcpointer,cl_rcpointer> cl_heap_weak_hashtable_from_rcpointer2_to_rcpointer;

typedef _cl_hashtable_iterator<cl_htentry_from_rcpointer2_to_rcpointer> cl_hashtable_from_rcpointer2_to_rcpointer_iterator;

struct cl_wht_from_rcpointer2_to_rcpointer : public cl_rcpointer {
	// Constructors.
	cl_wht_from_rcpointer2_to_rcpointer (bool (*maygc_htentry) (const cl_htentry_from_rcpointer2_to_rcpointer&));
	cl_wht_from_rcpointer2_to_rcpointer (const cl_wht_from_rcpointer2_to_rcpointer&);
	// Assignment operators.
	cl_wht_from_rcpointer2_to_rcpointer& operator= (const cl_wht_from_rcpointer2_to_rcpointer&);
	// Iterator.
	cl_hashtable_from_rcpointer2_to_rcpointer_iterator iterator () const
	{ return ((cl_heap_weak_hashtable_from_rcpointer2_to_rcpointer*)pointer)->iterator(); }
	// Lookup.
	cl_rcpointer * get (const cl_rcpointer& x, const cl_rcpointer& y) const;
	// Store.
	void put (const cl_rcpointer& x, const cl_rcpointer& y, const cl_rcpointer& z) const;
};

}  // namespace cln

#endif /* _CL_RCPOINTER2_HASHWEAK_RCPOINTER_H */
