// Exception types.

#ifndef _CL_EXCEPTION_H
#define _CL_EXCEPTION_H

#include <string>
#include <stdexcept>

namespace cln {

// Base class of all exception classes thrown by CLN.
class runtime_exception : public std::runtime_error {
public:
	runtime_exception ()
		: std::runtime_error(std::string()) {}
	explicit runtime_exception (const std::string & what)
		: std::runtime_error(what) {}
};

// Thrown when an assertion is violated.
class notreached_exception : public runtime_exception {
public:
	notreached_exception (const char* filename, int lineno);
};

// Thrown when a pole is encountered.
class division_by_0_exception : public runtime_exception {
public:
	division_by_0_exception ();
};

// Thrown when a conversion with As(TYPE)(VALUE) fails.
class as_exception : public runtime_exception {
public:
	as_exception (const class cl_number& obj, const char * typestring, const char * filename, int line);
};

}  // namespace cln

#endif /* _CL_EXCEPTION_H */
