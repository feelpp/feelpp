// class cl_ht_from_integer_to_rcobject.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/hash/cl_I_hash_rcobject.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/hash/cl_hash1.h"

namespace cln {

static void cl_hashtable_from_integer_to_rcobject_destructor (cl_heap* pointer)
{
	(*(cl_heap_hashtable_from_integer_to_rcobject*)pointer).~cl_heap_hashtable_from_integer_to_rcobject();
}

cl_class cl_class_hashtable_from_integer_to_rcobject = {
	cl_hashtable_from_integer_to_rcobject_destructor,
	0
};

// These are not inline, because they tend to duplicate a lot of template code.

cl_ht_from_integer_to_rcobject::cl_ht_from_integer_to_rcobject ()
{
	var cl_heap_hashtable_from_integer_to_rcobject* ht = new cl_heap_hashtable_from_integer_to_rcobject ();
	ht->refcount = 1;
	ht->type = &cl_class_hashtable_from_integer_to_rcobject;
	pointer = ht;
}

cl_rcobject * cl_ht_from_integer_to_rcobject::get (const cl_I& x) const
{
	return ((cl_heap_hashtable_from_integer_to_rcobject*)pointer)->get(x);
}

void cl_ht_from_integer_to_rcobject::put (const cl_I& x, const cl_rcobject& y) const
{
	((cl_heap_hashtable_from_integer_to_rcobject*)pointer)->put(x,y);
}

}  // namespace cln
