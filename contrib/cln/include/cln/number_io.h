// I/O of numbers.

#ifndef _CL_NUMBER_IO_H
#define _CL_NUMBER_IO_H

#include "cln/io.h"
#include "cln/number.h"
#include "cln/exception.h"

namespace cln {

// Input.

class read_number_exception : public runtime_exception {
public:
	explicit read_number_exception(const std::string & what)
		: runtime_exception(what) {}
};

// Finish with bad syntax.
class read_number_bad_syntax_exception : public read_number_exception {
public:
	read_number_bad_syntax_exception(const char * string, const char * string_limit);
};

// Finish with junk after the number.
class read_number_junk_exception : public read_number_exception {
public:
	read_number_junk_exception(const char * string_rest, const char * string, const char * string_limit);
};

// Finish with premature EOF.
class read_number_eof_exception : public read_number_exception {
public:
	read_number_eof_exception();
};

struct cl_read_flags;

}  // namespace cln

#endif /* _CL_NUMBER_IO_H */
