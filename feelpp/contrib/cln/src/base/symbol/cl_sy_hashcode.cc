// cln/symbol.hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/symbol.h"


// Implementation.

#include "base/cl_offsetof.h"

namespace cln {

#define declare_alignof(where,type)  \
  struct CONCAT(aligndummy,__LINE__) { char slot1; type slot2; }; \
  const uintptr_t where = offsetof(CONCAT(aligndummy,__LINE__), slot2);

uintptr_t hashcode (const cl_symbol& s)
{
	// Strings don't move in memory, so we can just take the address.
	declare_alignof(string_alignment,cl_heap_string);
	return (uintptr_t)(s.pointer)
	       / (string_alignment & -string_alignment); // divide by power of 2
}

}  // namespace cln
