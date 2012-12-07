// cl_I hash tables

#ifndef _CL_I_HASH_GCOBJECT_H
#define _CL_I_HASH_GCOBJECT_H

#include "cln/number.h"
#include "cln/integer.h"
#include "base/hash/cl_hash1.h"

namespace cln {

typedef cl_htentry1<cl_I,cl_gcobject> cl_htentry_from_integer_to_gcobject;

typedef cl_heap_hashtable_1<cl_I,cl_gcobject> cl_heap_hashtable_from_integer_to_gcobject;

typedef _cl_hashtable_iterator<cl_htentry_from_integer_to_gcobject> cl_hashtable_from_integer_to_gcobject_iterator;

struct cl_ht_from_integer_to_gcobject : public cl_gcpointer {
	// Constructors.
	cl_ht_from_integer_to_gcobject ();
	cl_ht_from_integer_to_gcobject (const cl_ht_from_integer_to_gcobject&);
	// Assignment operators.
	cl_ht_from_integer_to_gcobject& operator= (const cl_ht_from_integer_to_gcobject&);
	// Iterator.
	cl_hashtable_from_integer_to_gcobject_iterator iterator () const
	{ return ((cl_heap_hashtable_from_integer_to_gcobject*)pointer)->iterator(); }
	// Lookup.
	cl_gcobject * get (const cl_I& x) const;
	// Store.
	void put (const cl_I& x, const cl_gcobject& y) const;
};

}  // namespace cln

#endif /* _CL_I_HASH_GCOBJECT_H */
