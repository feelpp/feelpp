// read_real().
// This file contains a slimmed down version of read_complex().
// It does not pull in all the complex function code.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real_io.h"


// Implementation.

#include "cln/io.h"
#include "base/string/cl_spushstring.h"
#include "cln/input.h"

namespace cln {

// We read an entire token (or even more, if it begins with #C) into a
// buffer and then call read_real() on the buffer.

class pushstring_hack : public cl_spushstring {
public:
	char* start_pointer (void) { return buffer; }
	char* end_pointer (void) { return buffer+index; }
};

static bool number_char_p (char c)
{
	if ((c >= '0') && (c <= '9'))
		return true;
	if (((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z')))
		return true;
	switch (c) {
		case '+': case '-': case '.': case '_': case '/':
			return true;
		default:
			return false;
	}
}

const cl_R read_real (std::istream& stream, const cl_read_flags& flags)
{
	// One pre-allocated buffer. This reduces the allocation/free cost.
	static pushstring_hack buffer;

	var int c;
	// Skip whitespace at the beginning.
	loop {
		c = stream.get();
		if (stream.eof() || stream.fail()) goto eof;
		if ((c == ' ') || (c == '\t') || (c == '\n'))
			continue;
		else
			break;
	}
	// Found first non-whitespace character.
	// Numbers cannot cross lines. We can treat EOF and '\n' the same way.
	buffer.reset();
	if (c == '#') {
		if (!(flags.lsyntax & lsyntax_commonlisp))
			goto syntax1;
		buffer.push(c);
		// Read some digits, then a letter, then a token.
		loop {
			c = stream.get();
			if (stream.eof() || stream.fail()) goto eof;
			buffer.push(c);
			if ((c >= '0') && (c <= '9'))
				continue;
			else
				break;
		}
		if (!(((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z'))))
			goto syntax1;
		c = stream.get();
		if (stream.eof() || stream.fail()) goto eof;
	}
	// Read a number token.
	if (!number_char_p(c))
		goto syntax1;
	loop {
		buffer.push(c);
		c = stream.peek();  // Avoid fail state on EOF.
		if (stream.eof() || stream.fail() || !number_char_p(c))
			break;
		c = stream.get();
	}
	// Parse the number.
	return read_real(flags,
	                 buffer.start_pointer(), buffer.end_pointer(),
	                 NULL
	                );

	// Handle syntax error.
syntax1:	buffer.push(c);
	throw read_number_bad_syntax_exception(buffer.start_pointer(),buffer.end_pointer());

	// Handle premature EOF.
eof:	throw read_number_eof_exception();
}

}  // namespace cln
