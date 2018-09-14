// Property lists.

#ifndef _CL_PROPLIST_H
#define _CL_PROPLIST_H

#include "cln/symbol.h"
#include "cln/malloc.h"

namespace cln {

// The only extensible way to extend objects at runtime in an extensible
// and decentralized way (without having to modify the object's class)
// is to add a property table to every object.
// For the moment, only very few properties are planned, so lists should be
// enough. Since properties represent additional information about the object,
// there is no need for removing properties, so singly linked lists will be
// enough.

// This is the base class for all properties.
struct cl_property {
private:
	cl_property* next;
public:
	cl_symbol key;
	// Constructor.
	cl_property (const cl_symbol& k) : next (NULL), key (k) {}
	// Destructor.
	virtual ~cl_property () {}
	// Allocation and deallocation.
	void* operator new (size_t size) { return malloc_hook(size); }
	void operator delete (void* ptr) { free_hook(ptr); }
private:
	virtual void dummy ();
// Friend declarations. They are for the compiler. Just ignore them.
	friend class cl_property_list;
};
#define SUBCLASS_cl_property() \
	void* operator new (size_t size) { return malloc_hook(size); } \
	void operator delete (void* ptr) { free_hook(ptr); }

struct cl_property_list {
private:
	cl_property* list;
public:
	cl_property* get_property (const cl_symbol& key);
	void add_property (cl_property* new_property);
	// Constructor.
	cl_property_list () : list (NULL) {}
	// Destructor.
	~cl_property_list ();
};

}  // namespace cln

#endif /* _CL_PROPLIST_H */
