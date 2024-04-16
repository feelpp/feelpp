// cl_rcpointer hash tables

#ifndef _CL_RCPOINTER_HASHWEAK_RCPOINTER_H
#define _CL_RCPOINTER_HASHWEAK_RCPOINTER_H

#include "cln/object.h"

#include "base/hash/cl_hash1weak.h"

namespace cln {

// Equality.
static inline bool equal (const cl_rcpointer& x, const cl_rcpointer& y)
{ return (x.pointer == y.pointer); }

// Hash code. Luckily objects don't move around in memory.
inline uintptr_t hashcode (const cl_rcpointer& x)
{ return (uintptr_t)x.pointer; }

typedef cl_htentry1<cl_rcpointer,cl_rcpointer> cl_htentry_from_rcpointer_to_rcpointer;

typedef cl_heap_weak_hashtable_1<cl_rcpointer,cl_rcpointer> cl_heap_weak_hashtable_from_rcpointer_to_rcpointer;

typedef _cl_hashtable_iterator<cl_htentry_from_rcpointer_to_rcpointer> cl_hashtable_from_rcpointer_to_rcpointer_iterator;

struct cl_wht_from_rcpointer_to_rcpointer : public cl_rcpointer {
	// Constructors.
	cl_wht_from_rcpointer_to_rcpointer (bool (*maygc_htentry) (const cl_htentry_from_rcpointer_to_rcpointer&));
	cl_wht_from_rcpointer_to_rcpointer (const cl_wht_from_rcpointer_to_rcpointer&);
	// Assignment operators.
	cl_wht_from_rcpointer_to_rcpointer& operator= (const cl_wht_from_rcpointer_to_rcpointer&);
	// Iterator.
	cl_hashtable_from_rcpointer_to_rcpointer_iterator iterator () const
	{ return ((cl_heap_weak_hashtable_from_rcpointer_to_rcpointer*)pointer)->iterator(); }
	// Lookup.
	cl_rcpointer * get (const cl_rcpointer& x) const;
	// Store.
	void put (const cl_rcpointer& x, const cl_rcpointer& y) const;
};

}  // namespace cln

#endif /* _CL_RCPOINTER_HASHWEAK_RCPOINTER_H */
