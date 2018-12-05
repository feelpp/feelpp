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
//! @date 06 Aug 2017
//! @copyright 2017 Feel++ Consortium
//!

#define BOOST_TEST_MODULE gmsh
#include <testsuite.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/domain.hpp>

BOOST_AUTO_TEST_SUITE( gmsh )

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace Feel;
    auto m = unitSquare();
}

BOOST_AUTO_TEST_SUITE_END()
