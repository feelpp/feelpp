// Weak hash tables with 2 keys and a value

#ifndef _CL_HASH2WEAK_H
#define _CL_HASH2WEAK_H

#include "base/hash/cl_hash2.h"

namespace cln {

// This is a hash table in which an entry can be removed when a user-defined
// condition is fulfilled (e.g. the value is not referenced any more).
// We don't remove unused entries immediately, only when the hash table
// wants to grow. This way the hash table also serves as a cache.

// Requirements:
// - same as for hash2,
// - key1_type and key2_type must be subclasses of cl_gc[object|pointer],
// - value_type must be a subclass of cl_[gc|rc][object|pointer],
// - function maygc_htentry(const cl_htentry2<key1_type,key2_type,value_type>&);
//   must be filled in at runtime.

template <class key1_type, class key2_type, class value_type>
struct cl_heap_weak_hashtable_2 : public cl_heap_hashtable_2 <key1_type,key2_type,value_type> {
	// Allocation.
	void* operator new (size_t size) { return malloc_hook(size); }
	// Deallocation.
	void operator delete (void* ptr) { free_hook(ptr); }
public:
	// Function which tells when an unused entry may be garbage collected.
	bool (* const _maygc_htentry) (const cl_htentry2<key1_type,key2_type,value_type>&);
	// Constructor.
	cl_heap_weak_hashtable_2 (bool (*maygc_htentry) (const cl_htentry2<key1_type,key2_type,value_type>&))
		: cl_heap_hashtable_2 <key1_type,key2_type,value_type> (),
		  _maygc_htentry (maygc_htentry)
	{
		this->_garcol_fun = garcol;
	}
private:
	// Garbage collection.
	// Before growing the table, we check whether we can remove unused
	// entries.
	static bool garcol (cl_heap* _ht)
	{
		var cl_heap_weak_hashtable_2* ht = (cl_heap_weak_hashtable_2*)_ht;
		// Now ht->_garcol_fun = garcol.
		// It is not worth doing a garbage collection if the table
		// is small, say, has fewer than 100 entries.
		if (ht->_count < 100)
			return false;
		// Do a garbage collection.
		var intptr_t removed = 0;
		for (intptr_t i = 0; i < ht->_size; i++)
		    if (ht->_entries[i].next >= 0) {
			var cl_htentry2<key1_type,key2_type,value_type>& entry = ht->_entries[i].entry;
			if (ht->_maygc_htentry(entry)) {
				// This is hairy. We remove the entry and
				// free the value after its refcount has
				// dropped to zero. But in order to protect
				// against too early destruction
				// we have to temporarily increase the refcount.
				if (entry.val.pointer_p())
					entry.val.inc_pointer_refcount();
				ht->remove(entry.key1,entry.key2);
				if (entry.val.pointer_p()) {
					var cl_heap* p = entry.val.heappointer;
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
		var cl_heap_weak_hashtable_2* ht = (cl_heap_weak_hashtable_2*)_ht;
		// Now ht->_garcol_fun = garcol_nexttime.
		ht->_garcol_fun = cl_heap_weak_hashtable_2<key1_type,key2_type,value_type>::garcol;
		return false;
	}
};

}  // namespace cln

#endif /* _CL_HASH2WEAK_H */
