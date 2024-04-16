// Hash sets (hash tables with 1 key and no value)

#ifndef _CL_HASHSET_H
#define _CL_HASHSET_H

#include "base/hash/cl_hash.h"
#include "base/cl_iterator.h"

namespace cln {

// Requirements:
// - function  bool equal (key1_type,key1_type);
// - function  uintptr_t hashcode (key1_type);

template <class key1_type>
struct cl_htsetentry {
    ALLOCATE_ANYWHERE(cl_htsetentry)
    key1_type key;
    cl_htsetentry (const key1_type& k)
        : key (k) {}
};

template <class key1_type>
struct cl_heap_hashtable_set : public cl_heap_hashtable <cl_htsetentry <key1_type> > {
protected:
    // Abbreviations.
    typedef cl_heap_hashtable <cl_htsetentry <key1_type> > inherited;
    typedef typename inherited::htxentry htxentry;
public:
    // Allocation.
    void* operator new (size_t size) { return malloc_hook(size); }
    // Deallocation.
    void operator delete (void* ptr) { free_hook(ptr); }
public:
    // Lookup (htref alias gethash).
    bool get (const key1_type& key)
    {
        var intptr_t index = this->_slots[hashcode(key) % this->_modulus] - 1;
        while (index >= 0) {
            if (!(index < this->_size))
                throw runtime_exception();
            if (equal(key,this->_entries[index].entry.key))
                return true;
            index = this->_entries[index].next - 1;
        }
        return false;
    }
    // Store (htset alias puthash).
    void put (const key1_type& key)
    {
        var uintptr_t hcode = hashcode(key);
        // Search whether it is already there.
        {
            var intptr_t index = this->_slots[hcode % this->_modulus] - 1;
            while (index >= 0) {
                if (!(index < this->_size))
                    throw runtime_exception();
                if (equal(key,this->_entries[index].entry.key))
                    return;
                index = this->_entries[index].next - 1;
            }
        }
        // Put it into the table.
        prepare_store();
        var intptr_t hindex = hcode % this->_modulus; // _modulus may have changed!
        var intptr_t index = this->get_free_index();
        new (&this->_entries[index].entry) cl_htsetentry<key1_type> (key);
        this->_entries[index].next = this->_slots[hindex];
        this->_slots[hindex] = 1+index;
        this->_count++;
    }
    // Remove (htrem alias remhash).
    void remove (const key1_type& key)
    {
        var intptr_t* _index = &this->_slots[hashcode(key) % this->_modulus];
        while (*_index > 0) {
            var intptr_t index = *_index - 1;
            if (!(index < this->_size))
                throw runtime_exception();
            if (equal(key,this->_entries[index].entry.key)) {
                // Remove _entries[index].entry
                *_index = this->_entries[index].next;
                this->_entries[index].~htxentry();
                // The entry is now free.
                this->put_free_index(index);
                // That's it.
                this->_count--;
                return;
            }
            _index = &this->_entries[index].next;
        }
    }
    // Iterate through the table.
    // No stuff should be inserted into the table during the iteration,
    // or you may find yourself iterating over an entry vector which has
    // already been freed!
    // ??
private:
    // Prepare a store operation: make sure that the free list is non-empty.
    // This may change the table's size!
    void prepare_store ()
    {
        if (this->_freelist < -1)
            return;
        // Can we make room?
        if (this->_garcol_fun(this))
            if (this->_freelist < -1)
                return;
        // No! Have to grow the hash table.
        grow();
    }
    void grow ()
    {
        var intptr_t new_size = this->_size + (this->_size >> 1) + 1; // _size*1.5
        var intptr_t new_modulus = inherited::compute_modulus(new_size);
        var void* new_total_vector = malloc_hook(new_modulus*sizeof(intptr_t) + new_size*sizeof(htxentry));
        var intptr_t* new_slots = (intptr_t*) ((char*)new_total_vector + 0);
        var htxentry* new_entries = (htxentry *) ((char*)new_total_vector + new_modulus*sizeof(intptr_t));
        for (var intptr_t hi = new_modulus-1; hi >= 0; hi--)
            new_slots[hi] = 0;
        var intptr_t free_list_head = -1;
        for (var intptr_t i = new_size-1; i >= 0; i--) {
            new_entries[i].next = free_list_head;
            free_list_head = -2-i;
        }
        var htxentry* old_entries = this->_entries;
        for (var intptr_t old_index = 0; old_index < this->_size; old_index++)
            if (old_entries[old_index].next >= 0) {
                var key1_type& key = old_entries[old_index].entry.key;
                var intptr_t hindex = hashcode(key) % new_modulus;
                var intptr_t index = -2-free_list_head;
                free_list_head = new_entries[index].next;
                new (&new_entries[index].entry) cl_htsetentry<key1_type> (key);
                new_entries[index].next = new_slots[hindex];
                new_slots[hindex] = 1+index;
                old_entries[old_index].~htxentry();
            }
        free_hook(this->_total_vector);
        this->_modulus = new_modulus;
        this->_size = new_size;
        this->_freelist = free_list_head;
        this->_slots = new_slots;
        this->_entries = new_entries;
        this->_total_vector = new_total_vector;
    }
};

}  // namespace cln

#endif /* _CL_HASHSET_H */
