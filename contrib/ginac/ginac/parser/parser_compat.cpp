/** @file parser_compat.cpp
 *
 * Parser interface compatible with the old (bison/flex based) parser. */

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

#include "ex.h"
#include "idx.h"
#include "lst.h"
#include "parser.h"

#include <iostream>
#include <string>

namespace GiNaC {

static symtab make_symtab(const ex& l);

ptr<basic> ex::construct_from_string_and_lst(const std::string &s, const ex &l)
{
	static const bool strict = true;
	symtab syms = make_symtab(l);
	parser reader(syms, strict); 
	ex parsed_ex = reader(s);
	return parsed_ex.bp;
}

static std::string get_symbol_name(const ex & s);

static symtab make_symtab(const ex& l)
{
	symtab syms;
	if (is_exactly_a<lst>(l)) {
		for (std::size_t i = 0; i < l.nops(); i++) {
			const ex &o = l.op(i);
			if (is_a<symbol>(o) || (is_a<idx>(o) && is_a<symbol>(o.op(0))))
				syms[get_symbol_name(o)] = o;
		}
	}
	return syms;
}

static std::string get_symbol_name(const ex & s)
{
	if (is_a<symbol>(s))
		return ex_to<symbol>(s).get_name();
	else if (is_a<idx>(s) && is_a<symbol>(s.op(0)))
		return ex_to<symbol>(s.op(0)).get_name();
	else
		throw (std::invalid_argument("get_symbol_name(): unexpected expression type"));
}

} // namespace GiNaC
