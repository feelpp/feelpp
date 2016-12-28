/** @file fderivative.cpp
 *
 *  Implementation of abstract derivatives of functions. */

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

#include "fderivative.h"
#include "operators.h"
#include "archive.h"
#include "utils.h"

#include <iostream>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(fderivative, function,
  print_func<print_context>(&fderivative::do_print).
  print_func<print_csrc>(&fderivative::do_print_csrc).
  print_func<print_tree>(&fderivative::do_print_tree))

//////////
// default constructor
//////////

fderivative::fderivative()
{
}

//////////
// other constructors
//////////

fderivative::fderivative(unsigned ser, unsigned param, const exvector & args) : function(ser, args)
{
	parameter_set.insert(param);
}

fderivative::fderivative(unsigned ser, const paramset & params, const exvector & args) : function(ser, args), parameter_set(params)
{
}

fderivative::fderivative(unsigned ser, const paramset & params, exvector && v) : function(ser, std::move(v)), parameter_set(params)
{
}

//////////
// archiving
//////////

void fderivative::read_archive(const archive_node& n, lst& sym_lst)
{
	inherited::read_archive(n, sym_lst);
	unsigned i = 0;
	while (true) {
		unsigned u;
		if (n.find_unsigned("param", u, i))
			parameter_set.insert(u);
		else
			break;
		++i;
	}
}
GINAC_BIND_UNARCHIVER(fderivative);

void fderivative::archive(archive_node &n) const
{
	inherited::archive(n);
	auto i = parameter_set.begin(), end = parameter_set.end();
	while (i != end) {
		n.add_unsigned("param", *i);
		++i;
	}
}


//////////
// functions overriding virtual functions from base classes
//////////

void fderivative::print(const print_context & c, unsigned level) const
{
	// class function overrides print(), but we don't want that
	basic::print(c, level);
}

void fderivative::do_print(const print_context & c, unsigned level) const
{
	c.s << "D[";
	auto i = parameter_set.begin(), end = parameter_set.end();
	--end;
	while (i != end) {
		c.s << *i++ << ",";
	}
	c.s << *i << "](" << registered_functions()[serial].name << ")";
	printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());
}

void fderivative::do_print_csrc(const print_csrc & c, unsigned level) const
{
	c.s << "D_";
	auto i = parameter_set.begin(), end = parameter_set.end();
	--end;
	while (i != end)
		c.s << *i++ << "_";
	c.s << *i << "_" << registered_functions()[serial].name;
	printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());
}

void fderivative::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << class_name() << " "
	    << registered_functions()[serial].name << " @" << this
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << ", nops=" << nops()
	    << ", params=";
	auto i = parameter_set.begin(), end = parameter_set.end();
	--end;
	while (i != end)
		c.s << *i++ << ",";
	c.s << *i << std::endl;
	for (auto & i : seq)
		i.print(c, level + c.delta_indent);
	c.s << std::string(level + c.delta_indent, ' ') << "=====" << std::endl;
}

ex fderivative::eval() const
{
	// No parameters specified? Then return the function itself
	if (parameter_set.empty())
		return function(serial, seq);

	// If the function in question actually has a derivative, return it
	if (registered_functions()[serial].has_derivative() && parameter_set.size() == 1)
		return pderivative(*(parameter_set.begin()));

	return this->hold();
}

/** The series expansion of derivatives falls back to Taylor expansion.
 *  @see basic::series */
ex fderivative::series(const relational & r, int order, unsigned options) const
{
	return basic::series(r, order, options);
}

ex fderivative::thiscontainer(const exvector & v) const
{
	return fderivative(serial, parameter_set, v);
}

ex fderivative::thiscontainer(exvector && v) const
{
	return fderivative(serial, parameter_set, std::move(v));
}

/** Implementation of ex::diff() for derivatives. It applies the chain rule.
 *  @see ex::diff */
ex fderivative::derivative(const symbol & s) const
{
	ex result;
	for (size_t i=0; i<seq.size(); i++) {
		ex arg_diff = seq[i].diff(s);
		if (!arg_diff.is_zero()) {
			paramset ps = parameter_set;
			ps.insert(i);
			result += arg_diff * fderivative(serial, ps, seq);
		}
	}
	return result;
}

int fderivative::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<fderivative>(other));
	const fderivative & o = static_cast<const fderivative &>(other);

	if (parameter_set != o.parameter_set)
		return parameter_set < o.parameter_set ? -1 : 1;
	else
		return inherited::compare_same_type(o);
}

bool fderivative::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<fderivative>(other));
	const fderivative & o = static_cast<const fderivative &>(other);

	if (parameter_set != o.parameter_set)
		return false;
	else
		return inherited::is_equal_same_type(o);
}

bool fderivative::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<fderivative>(other));
	const fderivative & o = static_cast<const fderivative &>(other);

	return parameter_set == o.parameter_set && inherited::match_same_type(other);
}

/** Expose this object's derivative structure.
 *
 *  Parameter numbers occurring more than once stand for repeated
 *  differentiation with respect to that parameter. If a symbolic function
 *  f(x,y) is differentiated with respect to x, this method will return {0}.
 *  If f(x,y) is differentiated twice with respect to y, it will return {1,1}.
 *  (This corresponds to the way this object is printed.)
 *
 *  @return  multiset of function's parameter numbers that are abstractly
 *  differentiated. */
const paramset& fderivative::derivatives() const
{
	return parameter_set;
}


} // namespace GiNaC
