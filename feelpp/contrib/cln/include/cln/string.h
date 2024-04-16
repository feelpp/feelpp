// Strings.

#ifndef _CL_STRING_H
#define _CL_STRING_H

#include "cln/object.h"
#include "cln/io.h"
#include "cln/exception.h"
#include <cstring>

namespace cln {

struct cl_string;

// General, reference counted and garbage collected strings.
struct cl_heap_string : public cl_heap {
private:
	unsigned long length;	// length (in characters)
	char data[1];		// the characters, plus a '\0' at the end
	// Standard allocation disabled.
	void* operator new (size_t size) = delete;
	// Standard deallocation disabled.
	void operator delete (void* ptr) = delete;
	// No default constructor.
	cl_heap_string ();
private:
// Friend declarations. They are for the compiler. Just ignore them.
	friend class cl_string;
	friend cl_heap_string* cl_make_heap_string (unsigned long len);
	friend cl_heap_string* cl_make_heap_string (const char * s);
	friend cl_heap_string* cl_make_heap_string (const char * ptr, unsigned long len);
	friend const cl_string operator+ (const cl_string& str1, const cl_string& str2);
	friend const cl_string operator+ (const char* str1, const cl_string& str2);
	friend const cl_string operator+ (const cl_string& str1, const char* str2);
};

struct cl_string : public cl_gcpointer {
public:
	// Conversion to simple string.
	// NOTE! The resulting pointer is valid only as long as the string
	// is live, i.e. you must keep the string in a variable until you
	// are done with the pointer to the characters.
	const char * asciz () const
	{
		return &((cl_heap_string*)pointer)->data[0];
	}
	// Return the length (number of characters).
	unsigned long size() const
	{
		return ((cl_heap_string*)pointer)->length;
	}
	// Return a specific character.
	char operator[] (unsigned long i) const
	{
		if (!(i < size())) throw runtime_exception(); // Range check.
		return ((cl_heap_string*)pointer)->data[i];
	}
	// New ANSI C++ compilers also want the following.
	char operator[] (unsigned int i) const
		{ return operator[]((unsigned long)i); }
	char operator[] (long i) const
		{ return operator[]((unsigned long)i); }
	char operator[] (int i) const
		{ return operator[]((unsigned long)i); }
	// Constructors.
	cl_string ();
	cl_string (const cl_string&);
	cl_string (const char * s);
	cl_string (const char * ptr, unsigned long len);
	// Assignment operators.
	cl_string& operator= (const cl_string&);
	cl_string& operator= (const char *);
	// Private pointer manipulations.
	operator cl_heap_string* () const;
	cl_string (cl_heap_string* str) { pointer = str; }
	cl_string (cl_private_thing p) : cl_gcpointer (p) {}
};
CL_DEFINE_COPY_CONSTRUCTOR2(cl_string,cl_gcpointer)
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_string,cl_string)
inline cl_string::cl_string (const char * s)
{
	extern cl_heap_string* cl_make_heap_string (const char *);
	pointer = cl_make_heap_string(s);
}
inline cl_string& cl_string::operator= (const char * s)
{
	extern cl_heap_string* cl_make_heap_string (const char *);
	cl_heap_string* tmp = cl_make_heap_string(s);
	cl_dec_refcount(*this);
	pointer = tmp;
	return *this;
}

// Length.
inline unsigned long strlen (const cl_string& str)
{
	return str.size();
}
// Conversion to `const char *'.
inline const char * asciz (const char * s) { return s; }
inline const char * asciz (const cl_string& s) { return s.asciz(); }

// Comparison.
inline bool equal (const cl_string& str1, const cl_string& str2)
{
    return str1.size() == str2.size()
           && !strcmp(str1.asciz(), str2.asciz());
}
inline bool equal (const char * str1, const cl_string& str2)
{
    return !strcmp(str1, str2.asciz());
}
inline bool equal (const cl_string& str1, const char * str2)
{
    return !strcmp(str1.asciz(), str2);
}

// Private pointer manipulations. Never throw away a `struct cl_heap_string *'!
inline cl_string::operator cl_heap_string* () const
{
	cl_heap_string* hpointer = (cl_heap_string*)pointer;
	cl_inc_refcount(*this);
	return hpointer;
}
inline cl_string::cl_string ()
{
	static const cl_string cl_null_st(NULL, 0);
	pointer = (cl_heap_string*) cl_null_st;
}

// Hash code.
extern uintptr_t hashcode (const cl_string& str);

// Output.
extern void fprint (std::ostream& stream, const cl_string& str);
CL_DEFINE_PRINT_OPERATOR(cl_string)

// Input.

// Reads a line. Up to delim. The delimiter character is not placed in the
// resulting string. The delimiter character is kept in the input stream.
// If EOF is encountered, the stream's eofbit is set.
extern const cl_string cl_fget (std::istream& stream, char delim = '\n');

// Reads a line. Up to delim. The delimiter character is not placed in the
// resulting string. The delimiter character is extracted from the input stream.
// If EOF is encountered, the stream's eofbit is set.
extern const cl_string cl_fgetline (std::istream& stream, char delim = '\n');

// Like above, but only up to n-1 characters. If n-1 characters were read
// before the delimiter character was seen, the stream's failbit is set.
extern const cl_string cl_fget (std::istream& stream, int n, char delim = '\n');
extern const cl_string cl_fgetline (std::istream& stream, int n, char delim = '\n');

// Skips whitespace and then reads a non-whitespace string.
// If stream.width() is greater than 0, at most stream.width()-1 non-whitespace
// characters are read. When done, stream.width(0) is called.
// If EOF is encountered, the stream's eofbit is set.
extern std::istream& operator>> (std::istream& stream, cl_string& str);

// Runtime typing support.
extern cl_class cl_class_string;

// Debugging support.
#ifdef CL_DEBUG
extern int cl_string_debug_module;
CL_FORCE_LINK(cl_string_debug_dummy, cl_string_debug_module)
#endif

}  // namespace cln

#endif /* _CL_STRING_H */
