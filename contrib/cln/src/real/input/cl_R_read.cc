// read_real().
// This file contains a slimmed down version of read_complex().
// It does not pull in all the complex function code.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real_io.h"


// Implementation.

#include <cstring>
#include <sstream>
#include "cln/input.h"
#include "cln/rational_io.h"
#include "cln/integer_io.h"
#include "cln/float_io.h"
#include "cln/integer.h"
#include "integer/cl_I.h"
#include "float/cl_F.h"
#include "cln/exception.h"

#undef floor
#include <cmath>
#define floor cln_floor


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

const cl_R read_real (const cl_read_flags& flags, const char * string, const char * string_limit, const char * * end_of_parse)
{
	ASSERT((flags.syntax & ~(syntax_real|syntax_maybe_bad)) == 0);
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
			case '-': sign = ~sign;
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
	if (flags.syntax & syntax_float) {
		// Check for floating-point number syntax:
		// {'+'|'-'|} {digit}+ {'.' {digit}* | } expo {'+'|'-'|} {digit}+
		// {'+'|'-'|} {digit}* '.' {digit}+ expo {'+'|'-'|} {digit}+
		// {'+'|'-'|} {digit}* '.' {digit}+
		var const char * ptr = string;
		var const unsigned int float_base = 10;
		var cl_signean sign = 0;
		if (ptr == string_limit) goto not_float_syntax;
		switch (*ptr) {
			case '-': sign = ~sign;
			case '+': ptr++;
			default: break;
		}
		var const char * ptr_after_sign = ptr;
		var const char * ptr_after_intpart = skip_digits(ptr_after_sign,string_limit,float_base);
		var bool have_dot = false;
		var const char * ptr_before_fracpart = ptr_after_intpart;
		var const char * ptr_after_fracpart = ptr_after_intpart;
		ptr = ptr_after_intpart;
		if (ptr != string_limit)
		  if (*ptr == '.') {
			have_dot = true;
			ptr_before_fracpart = ptr+1;
			ptr_after_fracpart = skip_digits(ptr_before_fracpart,string_limit,float_base);
		}
		ptr = ptr_after_fracpart;
		var char exponent_marker;
		var bool have_exponent;
		var const char * ptr_in_exponent = ptr;
		var const char * ptr_after_exponent = ptr;
		if ((ptr == string_limit) || !(((*ptr >= '0') && (*ptr <= '9')) || ((*ptr >= 'A') && (*ptr <= 'Z') && (*ptr != 'I')) || ((*ptr >= 'a') && (*ptr <= 'z') && (*ptr != 'i')) || (*ptr == '.') || (*ptr == '/'))) {
			// No exponent.
			have_exponent = false;
			// Must have at least one fractional part digit.
			if (ptr_after_fracpart == ptr_before_fracpart) goto not_float_syntax;
			exponent_marker = 'E';
		} else {
			have_exponent = true;
			// Must have at least one digit.
			if (ptr_after_sign == ptr_after_intpart)
				if (ptr_after_fracpart == ptr_before_fracpart)
					goto not_float_syntax;
			exponent_marker = ((*ptr >= 'a') && (*ptr <= 'z') ? *ptr - 'a' + 'A' : *ptr);
			switch (exponent_marker) {
				case 'E':
				case 'S': case 'F': case 'D': case 'L':
					break;
				default:
					goto not_float_syntax;
			}
		}
		if (have_exponent) {
			ptr++;
			if (ptr == string_limit) goto not_float_syntax;
			switch (*ptr) {
				case '-':
				case '+': ptr++;
				default: break;
			}
			ptr_in_exponent = ptr;
			ptr_after_exponent = skip_digits(ptr_in_exponent,string_limit,10);
			if (ptr_after_exponent == ptr_in_exponent) goto not_float_syntax;
		}
		ptr = ptr_after_exponent;
		var const char * ptr_after_prec = ptr;
		var float_format_t prec;
		if ((ptr != string_limit) && (*ptr == '_')) {
			ptr++;
			ptr_after_prec = skip_digits(ptr,string_limit,10);
			if (ptr_after_prec == ptr) goto not_float_syntax;
			var cl_I prec1 = digits_to_I(ptr,ptr_after_prec-ptr,10);
			var uintC prec2 = cl_I_to_ulong(prec1);
			prec = (float_base==10 ? float_format(prec2)
			                       : (float_format_t)((uintC)((1+prec2)*::log((double)float_base)*1.442695041)+1)
			       );
		} else {
			switch (exponent_marker) {
				case 'S': prec = float_format_sfloat; break;
				case 'F': prec = float_format_ffloat; break;
				case 'D': prec = float_format_dfloat; break;
				case 'L': prec = flags.float_flags.default_lfloat_format; break;
				case 'E': prec = flags.float_flags.default_float_format; break;
				default: NOTREACHED
			}
			if (flags.float_flags.mantissa_dependent_float_format) {
				// Count the number of significant digits.
				ptr = ptr_after_sign;
				while (ptr < ptr_after_fracpart && (*ptr == '0' || *ptr == '.')) ptr++;
				var uintC num_significant_digits =
				  (ptr_after_fracpart - ptr) - (ptr_before_fracpart > ptr ? 1 : 0);
				var uintC prec2 = (num_significant_digits>=2 ? num_significant_digits-2 : 0);
				var float_format_t precx =
				  (float_base==10 ? float_format(prec2)
				                  : (float_format_t)((uintC)((1+prec2)*::log((double)float_base)*1.442695041)+1)
				  );
				if ((uintC)precx > (uintC)prec)
					prec = precx;
			}
		}
		floatformatcase(prec
		,	if (!(flags.syntax & syntax_sfloat)) goto not_float_syntax;
		,	if (!(flags.syntax & syntax_ffloat)) goto not_float_syntax;
		,	if (!(flags.syntax & syntax_dfloat)) goto not_float_syntax;
		,	unused len;
			if (!(flags.syntax & syntax_lfloat)) goto not_float_syntax;
		);
		at_end_of_parse(ptr_after_prec);
		return read_float(float_base,prec,sign,ptr_after_sign,0,ptr_after_fracpart-ptr_after_sign,ptr_after_exponent-ptr_after_sign,ptr_before_fracpart-ptr_after_sign);
	}
not_float_syntax:
	if (flags.syntax & syntax_maybe_bad) {
		ASSERT(end_of_parse);
		*end_of_parse = string;
		return 0; // dummy return
	}
	throw read_number_bad_syntax_exception(string,string_limit);
}

}  // namespace cln
