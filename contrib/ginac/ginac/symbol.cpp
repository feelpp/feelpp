/** @file symbol.cpp
 *
 *  Implementation of GiNaC's symbolic objects. */

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

#include "symbol.h"
#include "lst.h"
#include "archive.h"
#include "tostring.h"
#include "utils.h"
#include "hash_seed.h"
#include "inifcns.h"

#include <map>
#include <stdexcept>
#include <string>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(symbol, basic,
  print_func<print_context>(&symbol::do_print).
  print_func<print_latex>(&symbol::do_print_latex).
  print_func<print_tree>(&symbol::do_print_tree).
  print_func<print_python_repr>(&symbol::do_print_python_repr))

//////////
// default constructor
//////////

// symbol

symbol::symbol() : serial(next_serial++), name(""), TeX_name("")
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

// realsymbol

realsymbol::realsymbol() : symbol() { }

// possymbol

possymbol::possymbol() : realsymbol() { }

//////////
// other constructors
//////////

// public

// symbol

symbol::symbol(const std::string & initname) : serial(next_serial++),
	name(initname), TeX_name("")
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

symbol::symbol(const std::string & initname, const std::string & texname) :
	serial(next_serial++), name(initname), TeX_name(texname)
{
	setflag(status_flags::evaluated | status_flags::expanded);
}

// realsymbol
	
realsymbol::realsymbol(const std::string & initname) : symbol(initname) { }

realsymbol::realsymbol(const std::string & initname, const std::string & texname)
	: symbol(initname, texname) { }

// possymbol
	
possymbol::possymbol(const std::string & initname) : realsymbol(initname) { }

possymbol::possymbol(const std::string & initname, const std::string & texname) 
	: realsymbol(initname, texname) { }

//////////
// archiving
//////////

/** Read object from archive_node. */
void symbol::read_archive(const archive_node &n, lst &sym_lst)
{
	inherited::read_archive(n, sym_lst);
	serial = next_serial++;
	std::string tmp_name;
	n.find_string("name", tmp_name);

	// If symbol is in sym_lst, return the existing symbol
	for (lst::const_iterator it = sym_lst.begin(); it != sym_lst.end(); ++it) {
		if (is_a<symbol>(*it) && (ex_to<symbol>(*it).name == tmp_name)) {
			*this = ex_to<symbol>(*it);
			// XXX: This method is responsible for reading realsymbol
			// and possymbol objects too. But
			// basic::operator=(const basic& other)
			// resets status_flags::evaluated if other and *this are
			// of different types. Usually this is a good idea, but
			// doing this for symbols is wrong (for one, nothing is
			// going to set status_flags::evaluated, evaluation will
			// loop forever). Therefore we need to restore flags.
			setflag(status_flags::evaluated | status_flags::expanded);
			return;
		}
	}
	name = tmp_name;
	if (!n.find_string("TeXname", TeX_name))
		TeX_name = std::string("");
	setflag(status_flags::evaluated | status_flags::expanded);

	setflag(status_flags::dynallocated);
	sym_lst.append(*this);
}

/** Archive the object. */
void symbol::archive(archive_node &n) const
{
	inherited::archive(n);
	// XXX: we should not archive anonymous symbols.
	if (!name.empty())
		n.add_string("name", name);
	if (!TeX_name.empty())
		n.add_string("TeX_name", TeX_name);
}

//////////
// functions overriding virtual functions from base classes
//////////

/** Return default TeX name for symbol. This recognizes some greek letters. */
static const std::string& get_default_TeX_name(const std::string& name);

// public

std::string symbol::get_name() const
{
	if (name.empty()) {
		std::ostringstream s;
		s << "symbol" << serial;
		name = s.str();
	}
	return name;
}

// protected

void symbol::do_print(const print_context & c, unsigned level) const
{
	c.s << get_name();
}

void symbol::do_print_latex(const print_latex & c, unsigned level) const
{
	if (!TeX_name.empty())
		c.s << TeX_name;
	else if (!name.empty())
		c.s << get_default_TeX_name(name);
	else
		c.s << "symbol" << serial;
}

void symbol::do_print_tree(const print_tree & c, unsigned level) const
{
	c.s << std::string(level, ' ') << name << " (" << class_name() << ")" << " @" << this
	    << ", serial=" << serial
	    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
	    << ", domain=" << get_domain()
	    << std::endl;
}

void symbol::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << "('";
	if (!name.empty())
		c.s << name;
	else
		c.s << "symbol" << serial;
	if (!TeX_name.empty())
		c.s << "','" << TeX_name;
	c.s << "')";
}

bool symbol::info(unsigned inf) const
{
	switch (inf) {
		case info_flags::symbol:
		case info_flags::polynomial:
		case info_flags::integer_polynomial: 
		case info_flags::cinteger_polynomial: 
		case info_flags::rational_polynomial: 
		case info_flags::crational_polynomial: 
		case info_flags::rational_function: 
		case info_flags::expanded:
			return true;
		case info_flags::real:
			return get_domain() == domain::real || get_domain() == domain::positive;
		case info_flags::positive:
		case info_flags::nonnegative:
			return get_domain() == domain::positive;
		case info_flags::has_indices:
			return false;
	}
	return inherited::info(inf);
}

ex symbol::conjugate() const
{
	return conjugate_function(*this).hold();
}

ex symbol::real_part() const
{
	return real_part_function(*this).hold();
}

ex symbol::imag_part() const
{
	return imag_part_function(*this).hold();
}

bool symbol::is_polynomial(const ex & var) const
{
	return true;
}

// protected

/** Implementation of ex::diff() for single differentiation of a symbol.
 *  It returns 1 or 0.
 *
 *  @see ex::diff */
ex symbol::derivative(const symbol & s) const
{
	if (compare_same_type(s))
		return _ex0;
	else
		return _ex1;
}

int symbol::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<symbol>(other));
	const symbol *o = static_cast<const symbol *>(&other);
	if (serial==o->serial) return 0;
	return serial < o->serial ? -1 : 1;
}

bool symbol::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<symbol>(other));
	const symbol *o = static_cast<const symbol *>(&other);
	return serial==o->serial;
}

unsigned symbol::calchash() const
{
	unsigned seed = make_hash_seed(typeid(*this));
	hashvalue = golden_ratio_hash(seed ^ serial);
	setflag(status_flags::hash_calculated);
	return hashvalue;
}

//////////
// virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

/** Return default TeX name for symbol. This recognizes some greek letters. */
static const std::string& get_default_TeX_name(const std::string& name)
{
	static std::map<std::string, std::string> standard_names;
	static bool names_initialized = false;
	if (!names_initialized) {
		standard_names["alpha"] = std::string("\\alpha");
		standard_names["beta"] = std::string("\\beta");;
		standard_names["gamma"] = std::string("\\gamma");;
		standard_names["delta"] = std::string("\\delta");;
		standard_names["epsilon"] = std::string("\\epsilon");
		standard_names["varepsilon"] = std::string("\\varepsilon");
		standard_names["zeta"] = std::string("\\zeta");
		standard_names["eta" ] = std::string("\\eta" );
		standard_names["theta"] = std::string("\\theta");
		standard_names["vartheta"] = std::string("\\vartheta");
		standard_names["iota"] = std::string("\\iota");
		standard_names["kappa"] = std::string("\\kappa");
		standard_names["lambda"] = std::string("\\lambda");
		standard_names["mu"] = std::string("\\mu");
		standard_names["nu"] = std::string("\\nu");
		standard_names["xi"] = std::string("\\xi");
		standard_names["omicron"] = std::string("\\omicron");
		standard_names["pi"] = std::string("\\pi");
		standard_names["varpi"] = std::string("\\varpi");
		standard_names["rho"] = std::string("\\rho");
		standard_names["varrho"] = std::string("\\varrho");
		standard_names["sigma"] = std::string("\\sigma");
		standard_names["varsigma"] = std::string("\\varsigma");
		standard_names["tau"] = std::string("\\tau");
		standard_names["upsilon"] = std::string("\\upsilon");
		standard_names["phi"] = std::string("\\phi");
		standard_names["varphi"] = std::string("\\varphi");
		standard_names["chi"] = std::string("\\chi");
		standard_names["psi"] = std::string("\\psi");
		standard_names["omega"] = std::string("\\omega");
		standard_names["Gamma"] = std::string("\\Gamma");
		standard_names["Delta"] = std::string("\\Delta");
		standard_names["Theta"] = std::string("\\Theta");
		standard_names["Lambda"] = std::string("\\Lambda");
		standard_names["Xi"] = std::string("\\Xi");
		standard_names["Pi"] = std::string("\\Pi");
		standard_names["Sigma"] = std::string("\\Sigma");
		standard_names["Upsilon"] = std::string("\\Upsilon");
		standard_names["Phi"] = std::string("\\Phi");
		standard_names["Psi"] = std::string("\\Psi");
		standard_names["Omega"] = std::string("\\Omega");
		names_initialized = true;
	}
	std::map<std::string, std::string>::const_iterator it = standard_names.find(name);
	if (it != standard_names.end())
		return it->second;
	else
		return name;
}

GINAC_BIND_UNARCHIVER(symbol);
GINAC_BIND_UNARCHIVER(realsymbol);
GINAC_BIND_UNARCHIVER(possymbol);

//////////
// static member variables
//////////

// private

unsigned symbol::next_serial = 0;

} // namespace GiNaC
