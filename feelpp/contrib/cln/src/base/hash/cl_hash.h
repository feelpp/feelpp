// General hashtables

#ifndef _CL_HASH_H
#define _CL_HASH_H

#include "cln/object.h"
#include "cln/malloc.h"
#include "cln/exception.h"
#include "base/cl_iterator.h"

namespace cln {

const long htentry_last = 0; // means that there is no next entry

// These forward declarations are needed for Sun CC 3.0.1 and 4.0.1.
template <class htentry> struct _cl_hashtable_iterator;

template <class htentry>
struct cl_heap_hashtable : public cl_heap {
    friend struct _cl_hashtable_iterator<htentry>;
protected:
    typedef struct htxentry {
        long next;     // > 0: pseudo-list continues at next-1
                       // == 0: end of pseudo-list
                       // == -1: end of pseudo-free-list
                       // < -1: part of pseudo-free-list, continues at -next-2
        htentry entry; // if next >= 0
    } htxentry;
    long _modulus; // size of the primary entry table, > 0
    long _size;  // maximum number of entries
    long _count; // current number of entries
    long _freelist; // start of pseudo-free-list
    long * _slots;  // vector of length _modulus
    htxentry * _entries; // vector of length _size
    void* _total_vector;
    bool (*_garcol_fun) (cl_heap*); // Function to make room in the table.
                               // Putting some intelligent function here turns
                               // a normal hash table into a "weak" hash table.
public:
    // Allocation.
    void* operator new (size_t size) { return malloc_hook(size); }
    // Deallocation.
    void operator delete (void* ptr) { free_hook(ptr); }
    // Constructor: build a new, empty table.
    cl_heap_hashtable (long initial_size = 5) : cl_heap (),
        _size (initial_size), _count (0), _garcol_fun (no_garcol)
    {
        _modulus = compute_modulus(_size);
        _total_vector = malloc_hook(_modulus*sizeof(long) + _size*sizeof(htxentry));
        _slots = (long*) ((char*)_total_vector + 0);
        _entries = (htxentry *) ((char*)_total_vector + _modulus*sizeof(long));
        for (var long hi = _modulus-1; hi >= 0; hi--)
            _slots[hi] = 0;
        var long free_list_head = -1;
        for (var long i = _size-1; i >= 0; i--) {
            _entries[i].next = free_list_head;
            free_list_head = -2-i;
        }
        _freelist = free_list_head;
    }
    // Destructor.
    ~cl_heap_hashtable ()
    {
        for (long i = 0; i < _size; i++)
            if (_entries[i].next >= 0)
                _entries[i].~htxentry();
        free_hook(_total_vector);
    }
    // Count number of entries.
    long num_entries ()
    {
        #if 0
        var long n = 0;
        for (long i = 0; i < _size; i++)
            if (_entries[i].next >= 0)
                n++;
        return n;
        #else
        /* We already have an up-to-date count. */
        return _count;
        #endif
    }
    // Iterator.
    _cl_hashtable_iterator<htentry> iterator ();
protected:
    // Compute the modulus, given the maximum number of entries.
    static long compute_modulus (long size)
    {
        // It should be somewhat greater than size, since we want to
        // avoid collisions.
        // With N = size and M = modulus := k*size, the probability for a
        // * primary slot to be empty is
        //     (1-1/M)^N == exp(-N/M) == exp(-1/k).
        // * primary slot to carry a pseudo-list of length 1 is
        //     N 1/M (1-1/M)^(N-1) == exp(-N/M)*N/M == exp(-1/k)*1/k.
        // * primary slot to carry a pseudo-list of length >1 (collision) is
        //     1 - (1-1/M)^N - N 1/M (1-1/M)^(N-1)
        //     == 1 - exp(-N/M)*(1 + N/M) == 1 - (exp(-1/k)*(1+1/k)).
        // Sample values:
        //              = 0   = 1   > 1
        //   k = 1.0   0.37  0.37  0.26
        //   k = 1.5   0.51  0.34  0.14
        //   k = 2.0   0.61  0.30  0.09
        // I think k = 1.0 is reasonable.
        // Furthermore, we make sure that M is not divisible by 2, 3, 5.
        // Because in some applications, the hash codes are divisible
        // by 2 or 3, and if the modulus were divisible by this number,
        // only every second or every third primary slot would be filled,
        // resulting in many collisions.
        var long m = 1*size;
        // Make sure m is not divisible by 2.
        if ((m % 2) == 0)
            m++;
        // Make sure m is not divisible by 3.
        if ((m % 3) == 0)
            m += 2;
        // Make sure m is not divisible by 5.
        if ((m % 5) == 0) {
            m += 2;
            if ((m % 3) == 0)
                m += 2;
        }
        return m;
    }
    // Return the index of a free entry. Assumes the free list is non-empty.
    long get_free_index ()
    {
        // Check whether there is some in the free list.
        if (_freelist < -1) {
            var long index = -2-_freelist;
            _freelist = _entries[index].next;
            return index;
        }
        #if !(defined(__hppa__) && !defined(__GNUC__)) // workaround HP CC problem
        throw runtime_exception();
        #endif
        return -1; // dummy
    }
    // Put a free index into the free list.
    void put_free_index (long index)
    {
        _entries[index].next = _freelist;
        _freelist = -2-index;
    }
private:
    // Default function to make room in a hash table.
    static bool no_garcol (cl_heap* ht) { unused ht; return false; }
};

template <class htentry>
struct _cl_hashtable_iterator
  #if !(defined(__mips__) && !defined(__GNUC__)) // workaround SGI CC bug
    : cl_abstract_iterator<htentry>
  #endif
{
private:
    typename cl_heap_hashtable<htentry>::htxentry * _entries;
    long _index;
public:
    _cl_hashtable_iterator () : _entries (0), _index (-1) {}
public: /* ugh */
    _cl_hashtable_iterator (typename cl_heap_hashtable<htentry>::htxentry * e, long i)
        : _entries (e), _index (i)
    {
        do { _index--; }
           while (_index >= 0 && _entries[_index].next < 0);
    }
public:
    bool endp () { return (_index < 0); }
    htentry& next ()
    {
        if (_index < 0)
            throw runtime_exception();
        var long old_index = _index;
        do { _index--; }
           while (_index >= 0 && _entries[_index].next < 0);
        return _entries[old_index].entry;
    }
};

template <class htentry>
inline _cl_hashtable_iterator<htentry> cl_heap_hashtable<htentry>::iterator ()
{
    return _cl_hashtable_iterator<htentry>(_entries,_size);
}

}  // namespace cln

#endif /* _CL_HASH_H */
