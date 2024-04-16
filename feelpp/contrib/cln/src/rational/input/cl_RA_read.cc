// read_rational().
// This file contains a slimmed down version of read_real().
// It does not pull in all the floating-point, complex and transcendental
// function code.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational_io.h"


// Implementation.

#include <cstring>
#include <sstream>
#include "cln/input.h"
#include "cln/integer.h"
#include "cln/integer_io.h"
#include "integer/cl_I.h"
#include "cln/exception.h"

namespace cln {

// Step forward over all digits, to the end of string or to the next non-digit.
static const char * skip_digits (const char * ptr, const char * string_limit, unsigned int base)
{
	for ( ; ptr != string_limit; ptr++) {
		var char ch = *ptr;
		if ((ch >= '0') && (ch <= '9'))
			if (ch < '0' + (int)base)
				continue;
			else
				break;
		else {
			if (base <= 10)
				break;
			if (((ch >= 'A') && (ch < 'A'-10+(int)base))
			    || ((ch >= 'a') && (ch < 'a'-10+(int)base))
			   )
				continue;
			else
				break;
		}
	}
	return ptr;
}

#define at_end_of_parse(ptr)  \
  if (end_of_parse)							\
    { *end_of_parse = (ptr); }						\
  else									\
    { if ((ptr) != string_limit) { throw read_number_junk_exception((ptr),string,string_limit); } }

const cl_RA read_rational (const cl_read_flags& flags, const char * string, const char * string_limit, const char * * end_of_parse)
{
	ASSERT((flags.syntax & ~(syntax_rational|syntax_maybe_bad)) == 0);
	// If no string_limit is given, it defaults to the end of the string.
	if (!string_limit)
		string_limit = string + ::strlen(string);
	if (flags.syntax & syntax_rational) {
		// Check for rational number syntax.
		var unsigned int rational_base = flags.rational_base;
		var const char * ptr = string;
		if (flags.lsyntax & lsyntax_commonlisp) {
			if (ptr == string_limit) goto not_rational_syntax;
			if (*ptr == '#') {
				// Check for #b, #o, #x, #nR syntax.
				ptr++;
				if (ptr == string_limit) goto not_rational_syntax;
				switch (*ptr) {
				case 'b': case 'B':
					rational_base = 2; break;
				case 'o': case 'O':
					rational_base = 8; break;
				case 'x': case 'X':
					rational_base = 16; break;
				default:
					var const char * base_end_ptr =
						skip_digits(ptr,string_limit,10);
					if (base_end_ptr == ptr) goto not_rational_syntax;
					if (base_end_ptr == string_limit) goto not_rational_syntax;
					if (!((*base_end_ptr == 'r') || (*base_end_ptr == 'R')))
						goto not_rational_syntax;
					var cl_I base = read_integer(10,0,ptr,0,base_end_ptr-ptr);
					if (!((base >= 2) && (base <= 36))) {
						std::ostringstream buf;
						fprint(buf, "Base must be an integer in the range from 2 to 36, not ");
						fprint(buf, base);
						throw runtime_exception(buf.str());
					}
					rational_base = FN_to_UV(base); ptr = base_end_ptr;
					break;
				}
				ptr++;
			}
		}
		var const char * ptr_after_prefix = ptr;
		var cl_signean sign = 0;
		if (ptr == string_limit) goto not_rational_syntax;
		switch (*ptr) {
			case '-': sign = ~sign; // fallthrough
			case '+': ptr++;
			default: break;
		}
		var const char * ptr_after_sign = ptr;
		if (flags.syntax & syntax_integer) {
			// Check for integer syntax:  {'+'|'-'|} {digit}+ {'.'|}
			// Allow final dot only in Common Lisp syntax if there was no #<base> prefix.
			if ((flags.lsyntax & lsyntax_commonlisp) && (ptr_after_prefix == string)) {
				ptr = skip_digits(ptr_after_sign,string_limit,10);
				if (ptr != ptr_after_sign)
				  if (ptr != string_limit)
				    if (*ptr == '.') {
					ptr++;
					if ((ptr == string_limit) || !(((*ptr >= '0') && (*ptr <= '9')) || ((*ptr >= 'A') && (*ptr <= 'Z') && (*ptr != 'I')) || ((*ptr >= 'a') && (*ptr <= 'z') && (*ptr != 'i')) || (*ptr == '.') || (*ptr == '_') || (*ptr == '/'))) {
						at_end_of_parse(ptr);
						return read_integer(10,sign,ptr_after_sign,0,ptr-ptr_after_sign);
					}
				}
			}
			ptr = skip_digits(ptr_after_sign,string_limit,rational_base);
			if ((ptr == string_limit) || !(((*ptr >= '0') && (*ptr <= '9')) || ((*ptr >= 'A') && (*ptr <= 'Z') && (*ptr != 'I')) || ((*ptr >= 'a') && (*ptr <= 'z') && (*ptr != 'i')) || (*ptr == '.') || (*ptr == '_') || (*ptr == '/'))) {
				at_end_of_parse(ptr);
				return read_integer(rational_base,sign,ptr_after_sign,0,ptr-ptr_after_sign);
			}
		}
		if (flags.syntax & syntax_ratio) {
			// Check for ratio syntax: {'+'|'-'|} {digit}+ '/' {digit}+
			ptr = skip_digits(ptr_after_sign,string_limit,rational_base);
			if (ptr != ptr_after_sign)
			  if (ptr != string_limit)
			    if (*ptr == '/') {
				var const char * ptr_at_slash = ptr;
				ptr = skip_digits(ptr_at_slash+1,string_limit,rational_base);
				if (ptr != ptr_at_slash+1)
				  if ((ptr == string_limit) || !(((*ptr >= '0') && (*ptr <= '9')) || ((*ptr >= 'A') && (*ptr <= 'Z') && (*ptr != 'I')) || ((*ptr >= 'a') && (*ptr <= 'z') && (*ptr != 'i')) || (*ptr == '.') || (*ptr == '_') || (*ptr == '/'))) {
					at_end_of_parse(ptr);
					return read_rational(rational_base,sign,ptr_after_sign,0,ptr_at_slash-ptr_after_sign,ptr-ptr_after_sign);
				}
			}
		}
	}
not_rational_syntax:
	if (flags.syntax & syntax_maybe_bad) {
		ASSERT(end_of_parse);
		*end_of_parse = string;
		return 0; // dummy return
	}
	throw read_number_bad_syntax_exception(string,string_limit);
}

}  // namespace cln
