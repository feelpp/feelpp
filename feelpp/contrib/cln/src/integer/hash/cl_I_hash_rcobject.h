// cl_I hash tables

#ifndef _CL_I_HASH_RCOBJECT_H
#define _CL_I_HASH_RCOBJECT_H

#include "cln/number.h"
#include "cln/integer.h"
#include "base/hash/cl_hash1.h"

namespace cln {

typedef cl_htentry1<cl_I,cl_rcobject> cl_htentry_from_integer_to_rcobject;

typedef cl_heap_hashtable_1<cl_I,cl_rcobject> cl_heap_hashtable_from_integer_to_rcobject;

typedef _cl_hashtable_iterator<cl_htentry_from_integer_to_rcobject> cl_hashtable_from_integer_to_rcobject_iterator;

struct cl_ht_from_integer_to_rcobject : public cl_gcpointer {
	// Constructors.
	cl_ht_from_integer_to_rcobject ();
	cl_ht_from_integer_to_rcobject (const cl_ht_from_integer_to_rcobject&);
	// Assignment operators.
	cl_ht_from_integer_to_rcobject& operator= (const cl_ht_from_integer_to_rcobject&);
	// Iterator.
	cl_hashtable_from_integer_to_rcobject_iterator iterator () const
	{ return ((cl_heap_hashtable_from_integer_to_rcobject*)pointer)->iterator(); }
	// Lookup.
	cl_rcobject * get (const cl_I& x) const;
	// Store.
	void put (const cl_I& x, const cl_rcobject& y) const;
};

}  // namespace cln

#endif /* _CL_I_HASH_RCOBJECT_H */
