// Symbols.

#ifndef _CL_SYMBOL_H
#define _CL_SYMBOL_H

#include "cln/string.h"

namespace cln {

// Symbols are just strings, uniquified through a global hash table.

struct cl_symbol : public cl_rcpointer {
public:
	// Conversion to string.
	operator cl_string () const;
	// Constructors.
	cl_symbol (const cl_string&); // create or lookup a symbol from its name
	cl_symbol (const cl_symbol&);
	// Assignment operators.
	cl_symbol& operator= (const cl_symbol&);
	// Private pointer manipulations.
	cl_symbol (cl_private_thing p) : cl_rcpointer (p) {}
public: /* ugh */
	// Create a new symbol given its name.
	cl_symbol (struct hashuniq * null, const cl_string& s);
};
CL_DEFINE_COPY_CONSTRUCTOR2(cl_symbol,cl_rcpointer)
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_symbol,cl_symbol)

// A symbol points to a string, so to convert cl_symbol -> cl_string, we just
// take the pointer and put it into a cl_string.
inline cl_symbol::operator cl_string () const
{
	return cl_string(_as_cl_private_thing());
}

// Comparison.
inline bool equal (const cl_symbol& s1, const cl_symbol& s2)
{
	return (s1.pointer == s2.pointer);
}

// Hash code.
extern uintptr_t hashcode (const cl_symbol& s);

}  // namespace cln

#endif /* _CL_SYMBOL_H */
