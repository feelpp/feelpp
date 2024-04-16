// class cl_wht_from_integer_to_rcpointer.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/hash/cl_I_hashweak_rcpointer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/hash/cl_hash1weak.h"

namespace cln {

static void cl_weak_hashtable_from_integer_to_rcpointer_destructor (cl_heap* pointer)
{
	(*(cl_heap_weak_hashtable_from_integer_to_rcpointer*)pointer).~cl_heap_weak_hashtable_from_integer_to_rcpointer();
}

cl_class cl_class_weak_hashtable_from_integer_to_rcpointer = {
	cl_weak_hashtable_from_integer_to_rcpointer_destructor,
	0
};

// These are not inline, because they tend to duplicate a lot of template code.

cl_wht_from_integer_to_rcpointer::cl_wht_from_integer_to_rcpointer (bool (*maygc_htentry) (const cl_htentry_from_integer_to_rcpointer&))
{
	var cl_heap_weak_hashtable_from_integer_to_rcpointer* ht = new cl_heap_weak_hashtable_from_integer_to_rcpointer (maygc_htentry);
	ht->refcount = 1;
	ht->type = &cl_class_weak_hashtable_from_integer_to_rcpointer;
	pointer = ht;
}

cl_rcpointer * cl_wht_from_integer_to_rcpointer::get (const cl_I& x) const
{
	return ((cl_heap_weak_hashtable_from_integer_to_rcpointer*)pointer)->get(x);
}

void cl_wht_from_integer_to_rcpointer::put (const cl_I& x, const cl_rcpointer& y) const
{
	((cl_heap_weak_hashtable_from_integer_to_rcpointer*)pointer)->put(x,y);
}

}  // namespace cln
