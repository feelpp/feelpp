// cl_alloca().

#ifndef _CL_ALLOCA_H
#define _CL_ALLOCA_H

#include "base/cl_macros.h"
#include <cstdlib>

namespace cln {

// Allocating temporary data of arbitrary size.
// We prefer to allocate it on the stack instead of via malloc(), because
// that's fully inlinable and causes less cache misses. But the global stack
// size of applications is limited (typically 8 MB on Unix, 1 MB on Windows),
// and we don't want users of CLN to need to change these limits. Therefore
// we use stack allocation only for amounts < 64KB, and malloc() for larger
// blocks.
// Usage:
//   {CL_ALLOCA_STACK;
//    ...
//    ... = cl_alloca(...);
//    ...
//    ... = cl_small_alloca(...);
//    ...
//    ... = cl_alloca(...);
//    ...
//   }
// CL_ALLOCA_STACK declares that use of cl_alloca() and cl_small_alloca() is
// possible. Then cl_alloca() and cl_small_alloca() can be used an arbitrary
// number of times to get room.
// The allocated room's extent ends at the end of the { ... } block.
// In every C function CL_ALLOCA_STACK should only called once.
// Because of a gcc bug, functions using these macros shouldn't be declared
// inline.
// cl_alloca(size) fetches a block of size bytes.
// cl_small_alloca(size) fetches a block of size bytes, with size < 65536.
// CL_SMALL_ALLOCA_STACK is similar to CL_ALLOCA_STACK, but allows only
// the use of cl_small_alloca(), not cl_alloca().

// CL_ALLOCA_STACK creates a variable containing a linked list of pointers
// to be freed when the block is exited.

struct cl_alloca_header {
	cl_alloca_header* next;
	intptr_t usable_memory[1]; // "intptr_t" guarantees alignment
};

extern cl_alloca_header* cl_alloc_alloca_header (size_t size);
extern void cl_free_alloca_header (cl_alloca_header* pointer);

class cl_alloca_stack {
	cl_alloca_header* pointer;
public:
	cl_alloca_stack () { pointer = NULL; }
	~cl_alloca_stack () { if (pointer) cl_free_alloca_header(pointer); }
	void* push (cl_alloca_header* p) { p->next = pointer; pointer = p; return &p->usable_memory; }
};

#define CL_ALLOCA_STACK  \
  cl_alloca_stack _alloca_stack

#define CL_ALLOCA_MAX  65536

#if defined(__GNUC__) && !defined(__riscos) && !defined(__convex__)
  #define cl_alloca(size)  ((size) >= CL_ALLOCA_MAX ? _alloca_stack.push(cl_alloc_alloca_header(size)) : __builtin_alloca(size))
  #define cl_small_alloca(size)  __builtin_alloca(size)
  #define CL_SMALL_ALLOCA_STACK
#elif !defined(NO_ALLOCA) && !defined(__sparc__) && !defined(__sparc64__)
  #define cl_alloca(size)  ((size) >= CL_ALLOCA_MAX ? _alloca_stack.push(cl_alloc_alloca_header(size)) : alloca(size))
  #define cl_small_alloca(size)  alloca(size)
  #define CL_SMALL_ALLOCA_STACK
#else
  #define cl_alloca(size)  _alloca_stack.push(cl_alloc_alloca_header(size))
  #define cl_small_alloca(size)  _alloca_stack.push(cl_alloc_alloca_header(size))
  #define CL_SMALL_ALLOCA_STACK  CL_ALLOCA_STACK
#endif

// cl_alloc_array(type,size)
// cl_small_alloc_array(type,size)
// allocate an array with dynamic extent.
  #define cl_alloc_array(arrayeltype,arraysize)  \
    (arrayeltype*)cl_alloca((arraysize)*sizeof(arrayeltype))
  #define cl_small_alloc_array(arrayeltype,arraysize)  \
    (arrayeltype*)cl_small_alloca((arraysize)*sizeof(arrayeltype))

}  // namespace cln

#endif /* _CL_ALLOCA_H */
