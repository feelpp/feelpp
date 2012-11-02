// Conditions (a.k.a. exceptions)

#ifndef _CL_CONDITION_H
#define _CL_CONDITION_H

#include "cln/malloc.h"
#include "cln/io.h"

namespace cln {

struct cl_condition {
	// Allocation.
	void* operator new (size_t size) { return malloc_hook(size); }
	// Deallocation.
	void operator delete (void* ptr) { free_hook(ptr); }
	// Name.
	virtual const char * name () const = 0;
	// Print.
	virtual void print (std::ostream&) const = 0;
	// Virtual destructor.
	virtual ~cl_condition () = 0;
private:
	virtual void dummy ();
};
#define SUBCLASS_cl_condition() \
public:									  \
	/* Allocation. */						  \
	void* operator new (size_t size) { return malloc_hook(size); } \
	/* Deallocation. */						  \
	void operator delete (void* ptr) { free_hook(ptr); }

// Functions which want to raise a condition return a `cl_condition*'.
// The caller checks this value. NULL means no condition. The one who
// disposes the condition (handles it without resignalling it) should
// call `delete' on the condition pointer.

}  // namespace cln

#endif /* _CL_CONDITION_H */
