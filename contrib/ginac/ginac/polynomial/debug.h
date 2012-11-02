/** @file debug.h
 *
 *  Utility macros and functions for debugging. */

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

#ifndef GINAC_MOD_GCD_DEBUG_H
#define GINAC_MOD_GCD_DEBUG_H

#include "compiler.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#define DEBUG_PREFIX __func__ << ':' << __LINE__ << ": "
#define EXCEPTION_PREFIX std::string(__func__) + std::string(": ") +

#define Dout2(stream, msg)                                    \
do {                                                          \
	stream << DEBUG_PREFIX << msg << std::endl << std::flush; \
} while (0)
#define Dout(msg) Dout2(std::cout, msg)

#define bug3_on(condition, the_exception, msg)  \
do {                                            \
	if (unlikely(condition)) {                  \
		std::ostringstream err_stream;          \
		Dout2(err_stream, "BUG: " << msg);      \
		throw the_exception(err_stream.str());  \
	}                                           \
} while (0)

#define bug_on(condition, msg) bug3_on(condition, std::logic_error, msg)

#endif // ndef GINAC_MOD_GCD_DEBUG_H
