//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 17 Sep 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_PYEXPR_HPP
#define FEELPP_PYEXPR_HPP 1

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <pybind11/embed.h>
#pragma clang diagnostic pop
#include <iostream>
#include <map>
#include <vector>
#include <string>

namespace Feel {


namespace py = pybind11;
using namespace py::literals;

//!
//! evaluate python code \p pycode and retrieve dictionary of Feel++ expressions
//! expressions are of the form : \c f(x,y,z):x:y:z
//!
std::map<std::string,std::string> 
pyexpr( std::string const& pycode, std::vector<std::string> const& vars, std::map<std::string,std::string> const& locals );

//!
//! evaluate python code  from file \p pyfilename using local variables \p locals
//! and retrieve dictionary of Feel++ expressions expressions are of the form : \c f(x,y,z):x:y:z
//!
//! @c locals allows to set some input variables for the python code and stores the output variables
//! that are required by the c++ code.
//! @note that if the python expression to be retrieved is a sympy expression then the expression is
//! transformed into an expression that is then compiled and loaded as a plugin.
//!
//! we start with ta test without any input variable (eg @p locals is empty)
//! @code
//! // no locals
//! std::map<std::string,string> locals={};
//! pyexprFromFile( filename, locals );
//! @endcode
//! We define now some input data
//! @code
//! // no locals
//! std::map<std::string,string> locals={{"dim","3"}};
//! pyexprFromFile( filename, locals );
//! for( auto d : locals ) cout << d.first << ":" << d.second << std::endl;
//! @endcode
//!
void
pyexprFromFile( std::string const& pyfilename, std::map<std::string,std::string> & locals );

void
pyexprFromFile( std::string const& pyfilename, std::map<std::string,std::map<std::string,std::string>> & locals );

}

#endif
