// find_univpoly_ring().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/univpoly.h"


// Implementation.

#include "polynomial/cl_UP.h"

namespace cln {

// Create a new univariate polynomial ring with a named variable.

static inline cl_heap_univpoly_ring* cl_make_univpoly_ring (const cl_ring& r, const cl_symbol& varname)
{
	cl_heap_univpoly_ring* UPR = cl_make_univpoly_ring(r);
	UPR->add_property(new cl_varname_property(cl_univpoly_varname_key,varname));
	return UPR;
}

}  // namespace cln

// The table of univariate polynomial rings with named variable.
// A weak hash table (cl_ring,cl_symbol) -> cl_univpoly_ring.

#include "base/hash/cl_rcpointer2_hashweak_rcpointer.h"

namespace cln {

// An entry can be collected when the value (the ring) isn't referenced any more
// except from the hash table, and when the keys (the base ring and the name)
// are't referenced any more except from the hash table and the ring. Note that
// the ring contains exactly one reference to the base ring and exactly one
// reference to the name (on the property list).

static bool maygc_htentry (const cl_htentry_from_rcpointer2_to_rcpointer& entry)
{
	if (!entry.key1.pointer_p() || (entry.key1.heappointer->refcount == 2))
		if (!entry.key2.pointer_p() || (entry.key2.heappointer->refcount == 2))
			if (!entry.val.pointer_p() || (entry.val.heappointer->refcount == 1))
				return true;
	return false;
}

class named_univpoly_ring_cache
{
	static cl_wht_from_rcpointer2_to_rcpointer* univpoly_ring_table;
	static int count;
public:
	named_univpoly_ring_cache();
	~named_univpoly_ring_cache();

	inline cl_univpoly_ring* get_univpoly_ring(const cl_ring& r, const cl_symbol& v)
	{
		return (cl_univpoly_ring*) univpoly_ring_table->get(r,v);
	}
	inline void store_univpoly_ring(const cl_univpoly_ring& R)
	{
		univpoly_ring_table->put(R->basering(),
			                 ((cl_varname_property*)(R->get_property(cl_univpoly_varname_key)))->varname,
					 R);
	}
};

cl_wht_from_rcpointer2_to_rcpointer* named_univpoly_ring_cache::univpoly_ring_table = 0;
int named_univpoly_ring_cache::count = 0;

named_univpoly_ring_cache::named_univpoly_ring_cache()
{
	if (count++ == 0)
		univpoly_ring_table = new cl_wht_from_rcpointer2_to_rcpointer(maygc_htentry);
}


named_univpoly_ring_cache::~named_univpoly_ring_cache()
{
	if (--count == 0)
		delete univpoly_ring_table;
}

const cl_univpoly_ring find_univpoly_ring (const cl_ring& r, const cl_symbol& varname)
{
	static named_univpoly_ring_cache cache;
	var cl_univpoly_ring* ring_in_table = cache.get_univpoly_ring(r,varname);
	if (!ring_in_table) {
		var cl_univpoly_ring R = cl_make_univpoly_ring(r,varname);
		cache.store_univpoly_ring(R);
		ring_in_table = cache.get_univpoly_ring(r,varname);
		if (!ring_in_table)
			throw runtime_exception();
	}
	return *ring_in_table;
}

}  // namespace cln

