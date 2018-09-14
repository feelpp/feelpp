// cl_I hash tables

#ifndef _CL_I_HASH_GCPOINTER_H
#define _CL_I_HASH_GCPOINTER_H

#include "cln/number.h"
#include "cln/integer.h"
#include "base/hash/cl_hash1.h"

namespace cln {

typedef cl_htentry1<cl_I,cl_gcpointer> cl_htentry_from_integer_to_gcpointer;

typedef cl_heap_hashtable_1<cl_I,cl_gcpointer> cl_heap_hashtable_from_integer_to_gcpointer;

typedef _cl_hashtable_iterator<cl_htentry_from_integer_to_gcpointer> cl_hashtable_from_integer_to_gcpointer_iterator;

struct cl_ht_from_integer_to_gcpointer : public cl_gcpointer {
	// Constructors.
	cl_ht_from_integer_to_gcpointer ();
	cl_ht_from_integer_to_gcpointer (const cl_ht_from_integer_to_gcpointer&);
	// Assignment operators.
	cl_ht_from_integer_to_gcpointer& operator= (const cl_ht_from_integer_to_gcpointer&);
	// Iterator.
	cl_hashtable_from_integer_to_gcpointer_iterator iterator () const
	{ return ((cl_heap_hashtable_from_integer_to_gcpointer*)pointer)->iterator(); }
	// Lookup.
	cl_gcpointer * get (const cl_I& x) const;
	// Store.
	void put (const cl_I& x, const cl_gcpointer& y) const;
};

}  // namespace cln

#endif /* _CL_I_HASH_GCPOINTER_H */
