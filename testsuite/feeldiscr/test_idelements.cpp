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
//! @date 17 Jul 2017
//! @copyright 2017 Feel++ Consortium
//!
#define BOOST_TEST_MODULE test_idelements
#include <testsuite.hpp>

#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelmath/randint.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( testsuite_idelements )

BOOST_AUTO_TEST_CASE( listelements )
{
    auto mesh=unitCircle();
    auto rangeElt = elements(mesh);
    size_type nEltForTest = std::min( size_type(3), nelements( rangeElt ) );
    std::vector<size_type> l;
    if ( nEltForTest > 0 )
    {
        // take random element in the range of active element
        l = randint( nEltForTest, size_type(0), nEltForTest-1 );
        auto itbegin = begin( rangeElt );
        for ( size_type & k : l )
            k = unwrap_ref( *std::next( itbegin, k ) ).id();
    }
    LOG(INFO) << "l=" << l;
    // maybe some duplicated ids have been removed
    nEltForTest =  nelements( idelements( mesh, l ) );
    auto newmesh = createSubmesh( mesh, idelements( mesh, l ) );
    size_type nGlobalEltForTest = nEltForTest;
    mpi::all_reduce(Environment::worldComm(), mpi::inplace(nGlobalEltForTest), std::plus<size_type>());
    BOOST_CHECK_EQUAL( newmesh->numGlobalElements(), nGlobalEltForTest );
}

BOOST_AUTO_TEST_CASE( space_fromlistelements )
{
    auto mesh=unitCircle();
    auto rangeElt = elements(mesh);
    size_type nEltForTest = std::min( size_type(3), nelements( rangeElt ) );
    std::vector<size_type> l;
    if ( nEltForTest > 0 )
    {
        // take random element in the range of active element
        l = randint( nEltForTest, size_type(0), nEltForTest-1 );
        auto itbegin = begin( rangeElt );
        for ( size_type & k : l )
            k = unwrap_ref( *std::next( itbegin, k ) ).id();
    }
    // maybe some duplicated ids have been removed
    nEltForTest = nelements( idelements( mesh, l ) );
    auto newmesh = createSubmesh( mesh, idelements( mesh, l ) );
    auto Xh=Pdh<0>(newmesh);
    auto v = Xh->element(cst(1.));
    auto s = v.sum();
    size_type nGlobalEltForTest = nEltForTest;
    mpi::all_reduce(Environment::worldComm(), mpi::inplace(nGlobalEltForTest), std::plus<size_type>());
    BOOST_CHECK_EQUAL( s, nGlobalEltForTest );
}

BOOST_AUTO_TEST_SUITE_END()

