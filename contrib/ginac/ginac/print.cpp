/** @file print.cpp
 *
 *  Implementation of helper classes for expression output. */

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

#include "print.h"

#include <iostream>

namespace GiNaC {

/** Next free ID for print_context types. */
unsigned next_print_context_id = 0;


GINAC_IMPLEMENT_PRINT_CONTEXT(print_context, void)
GINAC_IMPLEMENT_PRINT_CONTEXT(print_dflt, print_context)
GINAC_IMPLEMENT_PRINT_CONTEXT(print_latex, print_context)
GINAC_IMPLEMENT_PRINT_CONTEXT(print_python, print_context)
GINAC_IMPLEMENT_PRINT_CONTEXT(print_python_repr, print_context)
GINAC_IMPLEMENT_PRINT_CONTEXT(print_tree, print_context)
GINAC_IMPLEMENT_PRINT_CONTEXT(print_csrc, print_context)
GINAC_IMPLEMENT_PRINT_CONTEXT(print_csrc_float, print_csrc)
GINAC_IMPLEMENT_PRINT_CONTEXT(print_csrc_double, print_csrc)
GINAC_IMPLEMENT_PRINT_CONTEXT(print_csrc_cl_N, print_csrc)

print_context::print_context()
	: s(std::cout), options(0) {}
print_context::print_context(std::ostream & os, unsigned opt)
	: s(os), options(opt) {}

print_dflt::print_dflt()
	: print_context(std::cout) {}
print_dflt::print_dflt(std::ostream & os, unsigned opt)
	: print_context(os, opt) {}

print_latex::print_latex()
	: print_context(std::cout) {}
print_latex::print_latex(std::ostream & os, unsigned opt)
	: print_context(os, opt) {}

print_python::print_python()
	: print_context(std::cout) {}
print_python::print_python(std::ostream & os, unsigned opt)
	: print_context(os, opt) {}

print_python_repr::print_python_repr()
	: print_context(std::cout) {}
print_python_repr::print_python_repr(std::ostream & os, unsigned opt)
	: print_context(os, opt) {}

print_tree::print_tree()
	: print_context(std::cout), delta_indent(4) {}
print_tree::print_tree(unsigned d)
	: print_context(std::cout), delta_indent(d) {}
print_tree::print_tree(std::ostream & os, unsigned opt, unsigned d)
	: print_context(os, opt), delta_indent(d) {}

print_csrc::print_csrc()
	: print_context(std::cout) {}
print_csrc::print_csrc(std::ostream & os, unsigned opt)
	: print_context(os, opt) {}

print_csrc_float::print_csrc_float()
	: print_csrc(std::cout) {}
print_csrc_float::print_csrc_float(std::ostream & os, unsigned opt)
	: print_csrc(os, opt) {}

print_csrc_double::print_csrc_double()
	: print_csrc(std::cout) {}
print_csrc_double::print_csrc_double(std::ostream & os, unsigned opt)
	: print_csrc(os, opt) {}

print_csrc_cl_N::print_csrc_cl_N()
	: print_csrc(std::cout) {}
print_csrc_cl_N::print_csrc_cl_N(std::ostream & os, unsigned opt)
	: print_csrc(os, opt) {}

} // namespace GiNaC
