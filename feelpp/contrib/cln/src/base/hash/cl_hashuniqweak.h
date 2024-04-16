// Weak hash tables for making objects unique

#ifndef _CL_HASHUNIQWEAK_H
#define _CL_HASHUNIQWEAK_H

#include "base/hash/cl_hashuniq.h"

namespace cln {

// This is a hashuniq table in which an entry can be removed when the
// value is not referenced any more.
// Best example: string -> symbol uniquification. When a symbol is not
// referenced any more (except from the hash table itself), the (string,symbol)
// pair may be removed from the hash table.
// We don't remove unused entries immediately, only when the hash table
// wants to grow. This way the hash table also serves as a cache.

// Requirements:
// - same as for hashuniq,
// - value_type must be a subclass of cl_[gc|rc][object|pointer].
// Note that since the reference counts are compared against 1, it doesn't
// make sense to have more than one weak hash table for the same value_type.

template <class key1_type, class value_type>
struct cl_heap_weak_hashtable_uniq : public cl_heap_hashtable_uniq <key1_type,value_type> {
	// Allocation.
	void* operator new (size_t size) { return malloc_hook(size); }
	// Deallocation.
	void operator delete (void* ptr) { free_hook(ptr); }
public:
	// Constructor.
	cl_heap_weak_hashtable_uniq ()
		: cl_heap_hashtable_uniq <key1_type,value_type> ()
	{
		this->_garcol_fun = garcol;
	}
private:
	// Garbage collection.
	// Before growing the table, we check whether we can remove unused
	// entries.
	static bool garcol (cl_heap* _ht)
	{
		var cl_heap_weak_hashtable_uniq* ht = (cl_heap_weak_hashtable_uniq*)_ht;
		// Now ht->_garcol_fun = garcol.
		// It is not worth doing a garbage collection if the table
		// is small, say, has fewer than 100 entries.
		if (ht->_count < 100)
			return false;
		// Do a garbage collection.
		var intptr_t removed = 0;
		for (intptr_t i = 0; i < ht->_size; i++)
		    if (ht->_entries[i].next >= 0) {
			var value_type& v = ht->_entries[i].entry.val;
			if (!v.pointer_p() || (v.heappointer->refcount == 1)) {
				// This is hairy. We remove the entry and
				// free the value after its refcount has
				// dropped to zero. But in order to protect
				// against too early destruction (depending on
				// how the C++ compiler optimizes hashkey())
				// we have to temporarily increase the refcount.
				if (v.pointer_p())
					v.inc_pointer_refcount();
				ht->remove(hashkey(v));
				if (v.pointer_p()) {
					var cl_heap* p = v.heappointer;
					if (!(--p->refcount == 0)) throw runtime_exception();
					cl_free_heap_object(p);
				}
				removed++;
			}
		    }
		if (removed == 0)
			// Unsuccessful. Let the table grow immediately.
			return false;
		else if (2*removed < ht->_count) {
			// Table shrank by less than a factor of 1/1.5.
			// Don't expand the table now, but expand it next time.
			ht->_garcol_fun = garcol_nexttime;
			return true;
		} else {
			// Table shrank much. Don't expand the table now,
			// and try a GC next time.
			return true;
		}
	}
	static bool garcol_nexttime (cl_heap* _ht)
	{
		var cl_heap_weak_hashtable_uniq* ht = (cl_heap_weak_hashtable_uniq*)_ht;
		// Now ht->_garcol_fun = garcol_nexttime.
		ht->_garcol_fun = garcol;
		return false;
	}
};

}  // namespace cln

#endif /* _CL_HASHUNIQWEAK_H */
