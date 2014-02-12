// cached_power().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

namespace cln {

// Table: For each base b (2 <= b <= 36), store k and b^k where k is the largest
// integer such that b^k < 2^intDsize, i.e. k == floor(log(2^intDsize-1,b)).
struct power_table_entry {
	uintC k;
	uintD b_to_the_k;
};
extern const power_table_entry power_table [36-2+1];

// Table: contains for each base b (2 <= b <= 36) either NULL or an array of
// lazily computed b^(k*2^i) and maybe 1/b^(k*2^i).
//#define MUL_REPLACES_DIV
struct cached_power_table_entry {
	ALLOCATE_ANYWHERE(cached_power_table_entry)
	cl_I base_pow; // 0 or b^(k*2^i)
#ifdef MUL_REPLACES_DIV
	cl_I inv_base_pow; // if base_pow: floor(2^(2*integer_length(base_pow))/base_pow)
#endif
};

struct cached_power_table {
	cached_power_table_entry element[40];
	// Constructor and destructor - nothing special.
	cached_power_table () {}
	~cached_power_table () {}
	// Allocation and deallocation.
	void* operator new (size_t size) { return malloc_hook(size); }
	void operator delete (void* ptr) { free_hook(ptr); }
};

extern cached_power_table* ctable [36-2+1];

const cached_power_table_entry * cached_power (uintD base, uintL i);

}  // namespace cln
