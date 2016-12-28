/** @file wildcard.cpp
 *
 *  Implementation of GiNaC's wildcard objects. */

/*
 *  GiNaC Copyright (C) 1999-2016 Johannes Gutenberg University Mainz, Germany
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

#include "wildcard.h"
#include "archive.h"
#include "utils.h"
#include "hash_seed.h"

#include <iostream>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(wildcard, basic,
  print_func<print_context>(&wildcard::do_print).
  print_func<print_tree>(&wildcard::do_print_tree).
  print_func<print_python_repr>(&wildcard::do_print_python_repr))

//////////
// default constructor
//////////

wildcard::wildcard() : label(0)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// other constructors
//////////

wildcard::wildcard(unsigned l) : label(l)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

//////////
// archiving
//////////

void wildcard::read_archive(const archive_node& n, lst& sym_lst)
{
	inherited::read_archive(n, sym_lst);
	n.find_unsigned("label", label);
	setflag(status_flags::evaluated | status_flags::expanded);
}
GINAC_BIND_UNARCHIVER(wildcard);

void wildcard::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_unsigned("label", label);
}

//////////
// functions overriding virtual functions from base classes
//////////

int wildcard::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<wildcard>(other));
	const wildcard &o = static_cast<const wildcard &>(other);

	if (label == o.label)
		return 0;
	else
		return label < o.label ? -1 : 1;
}

void wildcard::do_print(const print_context & c, unsigned level) const
{
	c.s << "$" << label;
}

void wildcard::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << "(" << label << ")" << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << std::endl;
}

void wildcard::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << '(' << label << ')';
}

unsigned wildcard::calchash() const
{
	// this is where the schoolbook method
	// (golden_ratio_hash(typeid(*this).name()) ^ label)
	// is not good enough yet...
	unsigned seed = make_hash_seed(typeid(*this));
	hashvalue = golden_ratio_hash(seed ^ label);
	setflag(status_flags::hash_calculated);
	return hashvalue;
}

bool wildcard::match(const ex & pattern, exmap& repl_lst) const
{
	// Wildcards must match each other exactly (this is required for
	// subs() to work properly because in the final step it substitutes
	// all wildcards by their matching expressions)
	return is_equal(ex_to<basic>(pattern));
}

bool haswild(const ex & x)
{
	if (is_a<wildcard>(x))
		return true;
	for (size_t i=0; i<x.nops(); ++i)
		if (haswild(x.op(i)))
			return true;
	return false;
}

} // namespace GiNaC
