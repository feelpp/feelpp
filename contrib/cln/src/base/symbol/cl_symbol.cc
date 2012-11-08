// class cl_symbol.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/symbol.h"


// Implementation.

#include "base/hash/cl_hashuniqweak.h"

namespace cln {

inline const cl_string hashkey (const cl_symbol& sym)
{
	return (cl_string)sym;
}

// A symbol points to a string, so to convert cl_string -> cl_symbol, we just
// take the pointer and put it into a cl_symbol.
inline cl_symbol::cl_symbol (struct hashuniq * null, const cl_string& s)
	: cl_rcpointer (as_cl_private_thing(s)) { unused null; }

typedef cl_htuniqentry<cl_string,cl_symbol> cl_htentry_from_string_to_symbol;

typedef cl_heap_weak_hashtable_uniq<cl_string,cl_symbol> cl_heap_hashtable_from_string_to_symbol;

typedef _cl_hashtable_iterator<cl_htentry_from_string_to_symbol> cl_hashtable_from_string_to_symbol_iterator;

static void cl_hashtable_from_string_to_symbol_destructor (cl_heap* pointer)
{
#if (defined(__mips__) || defined(__mips64__)) && !defined(__GNUC__) // workaround SGI CC bug
	(*(cl_heap_hashtable_from_string_to_symbol*)pointer).~cl_heap_weak_hashtable_uniq();
#else
	(*(cl_heap_hashtable_from_string_to_symbol*)pointer).~cl_heap_hashtable_from_string_to_symbol();
#endif
}


struct cl_ht_from_string_to_symbol : public cl_gcpointer {
	// Constructors.
	cl_ht_from_string_to_symbol ();
	cl_ht_from_string_to_symbol (const cl_ht_from_string_to_symbol&);
	// Assignment operators.
	cl_ht_from_string_to_symbol& operator= (const cl_ht_from_string_to_symbol&);
	// Iterator.
	cl_hashtable_from_string_to_symbol_iterator iterator () const
	{ return ((cl_heap_hashtable_from_string_to_symbol*)pointer)->iterator(); }
	// Lookup.
	cl_symbol * get (const cl_string& s) const;
	// Store.
	void put (const cl_string& s) const;
};

// These are not inline, because they tend to duplicate a lot of template code.

cl_ht_from_string_to_symbol::cl_ht_from_string_to_symbol ()
{
	static const cl_class cl_class_hashtable_from_string_to_symbol = {
		cl_hashtable_from_string_to_symbol_destructor,
		0
	};
	var cl_heap_hashtable_from_string_to_symbol* ht = new cl_heap_hashtable_from_string_to_symbol ();
	ht->refcount = 1;
	ht->type = &cl_class_hashtable_from_string_to_symbol;
	pointer = ht;
}

cl_symbol * cl_ht_from_string_to_symbol::get (const cl_string& s) const
{
	return ((cl_heap_hashtable_from_string_to_symbol*)pointer)->get(s);
}

void cl_ht_from_string_to_symbol::put (const cl_string& s) const
{
	((cl_heap_hashtable_from_string_to_symbol*)pointer)->put(s);
}

// The global symbol table.
class global_symbol_table
{
	static int count;
	static cl_ht_from_string_to_symbol* symbol_table;
public:
	inline cl_symbol* get(const cl_string& s)
	{
		return symbol_table->get(s);
	}

	inline void put(const cl_string& s)
	{
		symbol_table->put(s);
	}

	global_symbol_table();
	~global_symbol_table();
};

int global_symbol_table::count = 0;
cl_ht_from_string_to_symbol* global_symbol_table::symbol_table;

global_symbol_table::global_symbol_table()
{
	if (count++ == 0)
		symbol_table = new cl_ht_from_string_to_symbol();
}

global_symbol_table::~global_symbol_table()
{
	if (--count == 0)
		delete symbol_table;
}

// Create or lookup a symbol from its name.
cl_symbol::cl_symbol (const cl_string& s)
{
	static global_symbol_table symbol_table;
	var cl_symbol * sym_in_table;
	sym_in_table = symbol_table.get(s);
	if (!sym_in_table) {
		symbol_table.put(s);
		sym_in_table = symbol_table.get(s);
		if (!sym_in_table)
			throw runtime_exception();
	}
	var cl_heap* p = sym_in_table->heappointer;
	cl_inc_pointer_refcount(p);
	pointer = p;
}

}  // namespace cln

