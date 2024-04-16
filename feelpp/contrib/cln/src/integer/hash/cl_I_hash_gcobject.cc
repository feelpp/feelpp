// class cl_ht_from_integer_to_gcobject.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/hash/cl_I_hash_gcobject.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/hash/cl_hash1.h"

namespace cln {

static void cl_hashtable_from_integer_to_gcobject_destructor (cl_heap* pointer)
{
	(*(cl_heap_hashtable_from_integer_to_gcobject*)pointer).~cl_heap_hashtable_from_integer_to_gcobject();
}

cl_class cl_class_hashtable_from_integer_to_gcobject = {
	cl_hashtable_from_integer_to_gcobject_destructor,
	0
};

// These are not inline, because they tend to duplicate a lot of template code.

cl_ht_from_integer_to_gcobject::cl_ht_from_integer_to_gcobject ()
{
	var cl_heap_hashtable_from_integer_to_gcobject* ht = new cl_heap_hashtable_from_integer_to_gcobject ();
	ht->refcount = 1;
	ht->type = &cl_class_hashtable_from_integer_to_gcobject;
	pointer = ht;
}

cl_gcobject * cl_ht_from_integer_to_gcobject::get (const cl_I& x) const
{
	return ((cl_heap_hashtable_from_integer_to_gcobject*)pointer)->get(x);
}

void cl_ht_from_integer_to_gcobject::put (const cl_I& x, const cl_gcobject& y) const
{
	((cl_heap_hashtable_from_integer_to_gcobject*)pointer)->put(x,y);
}

}  // namespace cln
