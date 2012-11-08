/** @file error_report.h
 *
 *  Macro for additional debugging output.
 */

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

#ifndef GINAC_CHECK_ERROR_REPORT_H
#define GINAC_CHECK_ERROR_REPORT_H

#include <sstream>
#include <stdexcept>

#define cbug_on(cond, what)				\
do {							\
if (cond) {						\
	std::ostringstream err_stream;			\
	err_stream << __FILE__ << ':' << __LINE__	\
		   << what;				\
	throw std::logic_error(err_stream.str());	\
}							\
} while (0)

#endif // ndef GINAC_CHECK_ERROR_REPORT_H
