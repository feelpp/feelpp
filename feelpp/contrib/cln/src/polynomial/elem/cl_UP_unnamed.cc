// find_univpoly_ring().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/univpoly.h"


// Implementation.

#include "polynomial/cl_UP.h"

// The table of univariate polynomial rings without named variable.
// A weak hash table cl_ring -> cl_univpoly_ring.
// (It could also be a weak hashuniq table cl_ring -> cl_univpoly_ring.)

#include "base/hash/cl_rcpointer_hashweak_rcpointer.h"

namespace cln {

// An entry can be collected when the value (the ring) isn't referenced any more
// except from the hash table, and when the key (the base ring) isn't referenced
// any more except from the hash table and the ring. Note that the ring contains
// exactly one reference to the base ring.

static bool maygc_htentry (const cl_htentry_from_rcpointer_to_rcpointer& entry)
{
	if (!entry.key.pointer_p() || (entry.key.heappointer->refcount == 2))
		if (!entry.val.pointer_p() || (entry.val.heappointer->refcount == 1))
			return true;
	return false;
}

class univpoly_ring_cache
{
	static cl_wht_from_rcpointer_to_rcpointer* univpoly_ring_table;
	static int count;
public:
	inline cl_univpoly_ring* get_univpoly_ring(const cl_ring& r)
	{
		return (cl_univpoly_ring*) univpoly_ring_table->get(r);
	}
	inline void store_univpoly_ring(const cl_univpoly_ring& R)
	{
		univpoly_ring_table->put(R->basering(), R);
	}
	univpoly_ring_cache();
	~univpoly_ring_cache();
};

cl_wht_from_rcpointer_to_rcpointer* univpoly_ring_cache::univpoly_ring_table = 0;
int univpoly_ring_cache::count = 0;

univpoly_ring_cache::univpoly_ring_cache()
{
	if (count++ == 0)
		univpoly_ring_table = new cl_wht_from_rcpointer_to_rcpointer(maygc_htentry);
}

univpoly_ring_cache::~univpoly_ring_cache()
{
	if (--count == 0)
		delete univpoly_ring_table;
}

const cl_univpoly_ring find_univpoly_ring (const cl_ring& r)
{
	static univpoly_ring_cache cache;
	var cl_univpoly_ring* ring_in_table = cache.get_univpoly_ring(r);
	if (!ring_in_table) {
		var cl_univpoly_ring R = cl_make_univpoly_ring(r);
		cache.store_univpoly_ring(R);
		ring_in_table = cache.get_univpoly_ring(r);
		if (!ring_in_table)
			throw runtime_exception();
	}
	return *ring_in_table;
}

}  // namespace cln

