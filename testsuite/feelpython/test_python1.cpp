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

#include <pybind11/embed.h>
#include <iostream>
#include <sstream>
#define BOOST_TEST_MODULE submesh testsuite
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpython/pyexpr.hpp>

namespace py = pybind11;
using namespace py::literals;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( python )

BOOST_AUTO_TEST_CASE( t1 )
{
    py::scoped_interpreter guard{};

    auto locals = py::dict("f"_a="[x**2,y**2,z**2]");
    std::ostringstream ostr;
    ostr <<  "from sympy2ginac import *\n"
         << "f=[x**2,y**2,z**2];\n"
         << "laplacian_f=toginac(laplacian(f,[x,y,z]),[x,y,z]);\n"
         << "print('lap(f)=', laplacian_f)\n";
    
    py::exec(ostr.str().c_str(), py::globals(), locals);

    
    auto laplacian_f = locals["laplacian_f"].cast<std::string>();
    BOOST_CHECK_EQUAL( laplacian_f, "{2,2,2}:x:y:z" );
    BOOST_TEST_MESSAGE( "laplacian_f=" << laplacian_f );
}

BOOST_AUTO_TEST_CASE( t2 )
{
    py::scoped_interpreter guard{};
    std::ostringstream ostr;
    ostr <<  "from sympy2ginac import *\n"
         << "f=[x**2,y**2,z**2];\n"
         << "k=4;\n"
         << "kgrad_f=k*grad(f,[x,y,z]);\n"
         << "print('k*grad(f)=', kgrad_f)\n"
         << "div_kgrad_f=div(k*grad(f,[x,y,z]),[x,y,z]);\n"
         << "print('div(k*grad(f))=', div_kgrad_f)\n";

    std::map<std::string,std::string> locals;
    auto m = Feel::pyexpr( ostr.str(), {"kgrad_f","div_kgrad_f"}, locals );

    BOOST_CHECK_EQUAL( m.at("kgrad_f"), "{8*x,0,0,0,8*y,0,0,0,8*z}:x:y:z" );
    BOOST_CHECK_EQUAL( m.at("div_kgrad_f"), "{8,8,8}:x" );
    BOOST_TEST_MESSAGE( "div_kgrad_f=" << m.at("div_kgrad_f") );
}

BOOST_AUTO_TEST_SUITE_END()
