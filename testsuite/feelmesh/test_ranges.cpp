/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 17 août 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#define BOOST_TEST_MODULE range testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelmesh/ranges.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelcore/enumerate.hpp>
#include <feel/feelfilters/detail/meshcontiguousnumberingmapping.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

using namespace Feel;

BOOST_AUTO_TEST_SUITE( ranges_testsuite )

BOOST_AUTO_TEST_CASE( test_range1 )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_range1" );
    
    auto mesh = unitSquare();
    auto elts = elements(mesh);
    BOOST_CHECK_EQUAL( is_filter_v<decay_type<decltype(elts)>>, true );
    BOOST_CHECK_EQUAL( is_range_v<decay_type<decltype(elts)>>, true );
    for( auto const& ew : elts )
    {
        auto e = unwrap_ref(ew);
        BOOST_TEST_MESSAGE( "elt " << e.id() );
    }
    BOOST_CHECK(nelements(elts)>0);

    auto thefaces = faces(mesh);
    BOOST_CHECK_EQUAL( is_filter_v<decay_type<decltype(thefaces)>>, true );
    BOOST_CHECK_EQUAL( is_range_v<decay_type<decltype(thefaces)>>, true );
    for( auto const& fw : thefaces )
    {
        auto f = unwrap_ref(fw);
        BOOST_TEST_MESSAGE( "face " << f.id() );
    }
    BOOST_CHECK(nelements(thefaces)>0);
    
    auto bfaces = boundaryfaces(mesh);
    BOOST_CHECK_EQUAL( is_filter_v<decay_type<decltype(bfaces)>>, true );
    BOOST_CHECK_EQUAL( is_range_v<decay_type<decltype(bfaces)>>, true );
    for( auto const& fw : bfaces )
    {
        auto f = unwrap_ref(fw);
        BOOST_TEST_MESSAGE( "bface " << f.id() );
    }
    BOOST_CHECK(nelements(bfaces)>0);

    auto r = fragmentationMarkedElements(mesh);
    for ( auto const& [fragment, data] : r )
    {
        BOOST_TEST_MESSAGE( "fragment " << fragment );
        auto const& [elts, mark, name ] = data;
        BOOST_CHECK(nelements(elts)>0);

        for ( auto const& eltw : elts )
        {
            auto elt = unwrap_ref(eltw);
            BOOST_TEST_MESSAGE( "elt " << elt.id() << " mark " << mark );
        }
    }

    Feel::detail::MeshContiguousNumberingMapping<decay_type<decltype(mesh)>,float> mcnm( mesh.get() );
    BOOST_TEST_MESSAGE( "test_range1 done" );

}
#if 1
BOOST_AUTO_TEST_CASE( test_range2 )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_range2" );
    
    auto mesh = unitSquare();
    BOOST_TEST_MESSAGE("number of elements " << nelements(elements(mesh)) );
    double area = integrate( _range=elements(mesh), _expr=expr("1")).evaluate()(0,0);
    BOOST_CHECK_CLOSE( area, 1, 1e-10 );
    BOOST_TEST_MESSAGE( "test_range1 done" );
}
BOOST_AUTO_TEST_CASE( test_range3 )
{
    using namespace Feel;

    auto mesh = unitSquare();
    BOOST_TEST_MESSAGE("number of elements " << nelements(elements(mesh)) );
    for( int nparts : {1,2,4,8,16} )
    {
        auto vparts = partitionRange(elements(mesh),nparts);
        double area = integrate( _range=elements(mesh), _expr=expr("1")).evaluate()(0,0);
        std::vector<double> areas(nparts);
        double sum = 0;
        for ( auto [j, rj] : enumerate(vparts) )
        {
            areas[j] = integrate( _range=rj, _expr=expr("1")).evaluate()(0,0);
            sum+=areas[j];
        }
        BOOST_MESSAGE( fmt::format("nparts={}, area={}, areas={}", nparts, area, areas ) );
        BOOST_CHECK_CLOSE( area, sum, 1e-10 );
    }
}
#endif
BOOST_AUTO_TEST_SUITE_END()

