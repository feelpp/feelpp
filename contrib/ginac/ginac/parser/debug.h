/** @file debug.h
 *
 *  Debugging facilities for parser. */

/*
 *  GiNaC Copyright (C) 1999-2011 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef GINAC_PARSER_DEBUG_H
#define GINAC_PARSER_DEBUG_H

#include "compiler.h"

#include <iosfwd>
#include <sstream>
#include <stdexcept>

#ifndef __GNUC__
#if __STDC_VERSION__ < 199901L
#define __PRETTY_FUNCTION__ "<unknown>"
#else
#define __PRETTY_FUNCTION__ __func__
#endif
#endif

#define bail_out(exception, message) \
do { \
	std::ostringstream err; \
	err << __PRETTY_FUNCTION__ << "(" << __FILE__ << ':' << __LINE__ << ": "; \
	err << message; \
	throw exception(err.str()); \
} while (0)

#define Parse_error_(message) \
do { \
	std::ostringstream err; \
	err << "GiNaC: parse error at line " << scanner->line_num << \
		", column " << scanner->column << ": "; \
	err << message << std::endl; \
	err << '[' << __PRETTY_FUNCTION__ << "(" << __FILE__ << ':' << __LINE__ << ")]" << std::endl; \
	throw parse_error(err.str(), scanner->line_num, scanner->column); \
} while (0)

#define Parse_error(message) \
	Parse_error_(message << ", got: " << scanner->tok2str(token))

#define bug(message) bail_out(std::logic_error, message)

#define dout(condition, message) \
do { \
	if (unlikely(condition)) { \
		std::cerr << __PRETTY_FUNCTION__ \
			<< " (" << __FILE__ << ':' << __LINE__ << "): " \
			<< message << std::endl; \
	} \
} while (0)

#endif // GINAC_PARSER_DEBUG_H
