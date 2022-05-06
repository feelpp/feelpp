//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
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
//! @date 16 Apr 2017
//! @copyright 2017 Feel++ Consortium
//!
//!
#define BOOST_TEST_MODULE mesh
#include <feel/feelcore/testsuite.hpp>

#if !defined(USE_BOOST_TEST)
#include <feel/feelfilters/exporter.hpp>
#endif

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/makemesh.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( mesh )

BOOST_AUTO_TEST_CASE(test_shared_mesh)
{
    auto m1 = makeSharedMesh<Simplex<1>>();
    BOOST_CHECK( m1 );
    auto m2 = makeSharedMesh<Simplex<2>>();
    BOOST_CHECK( m2 );
    auto m3 = makeSharedMesh<Simplex<3>>();
    BOOST_CHECK( m3 );
    auto m4 = makeSharedMesh<Hypercube<1>>();
    BOOST_CHECK( m4 );
    auto m5 = makeSharedMesh<Hypercube<2>>();
    BOOST_CHECK( m5 );
    auto m6 = makeSharedMesh<Hypercube<3>>();
    BOOST_CHECK( m6 );
}

BOOST_AUTO_TEST_CASE(test_unique_mesh)
{
    auto m1 = makeUniqueMesh<Simplex<1>>();
    BOOST_CHECK( m1 );
    auto m2 = makeUniqueMesh<Simplex<2>>();
    BOOST_CHECK( m2 );
    auto m3 = makeUniqueMesh<Simplex<3>>();
    BOOST_CHECK( m3 );
    auto m4 = makeUniqueMesh<Hypercube<1>>();
    BOOST_CHECK( m4 );
    auto m5 = makeUniqueMesh<Hypercube<2>>();
    BOOST_CHECK( m5 );
    auto m6 = makeUniqueMesh<Hypercube<3>>();
    BOOST_CHECK( m6 );
}

BOOST_AUTO_TEST_SUITE_END()
