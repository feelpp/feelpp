/** @file parse_context.cpp
 *
 *  Implementation of the parser context. */

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

#include "parse_context.h"

#include <sstream>
#include <stdexcept>

namespace GiNaC {

ex
find_or_insert_symbol(const std::string& name, symtab& syms, const bool strict)
{
	symtab::const_iterator p = syms.find(name);
	if (p != syms.end())
		return p->second;

	if (strict)
		throw std::invalid_argument(
				std::string("find_or_insert_symbol: symbol \"") 
				+ name + "\" not found");

	const symbol sy(name);
	syms[name] = sy;
	return sy;
}

} // namespace GiNaC
