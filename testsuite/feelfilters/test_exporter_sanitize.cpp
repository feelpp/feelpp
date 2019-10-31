/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 15 janv. 2016
 
 Copyright (C) 2016 Feel++ Consortium
 
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
#define BOOST_TEST_MODULE test_exporter_sanitize
#include <feel/feelcore/testsuite.hpp>


#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feeldiscr/pch.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Exporter  Options" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_exporter" ,
                     "test_exporter" ,
                     "0.2",
                     "nD(n=2,3) test exporter",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@cemosis.fr", "" );
    return about;
}



 FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
 BOOST_AUTO_TEST_SUITE( inner_suite )


BOOST_AUTO_TEST_CASE( test_1 )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>(mesh);
    auto u = Xh->element(Px());
    auto v = Xh->element(Py());
    auto e = exporter(_mesh=mesh);
    e->add( "u ;)[ ", u );
    e->add( "v**v ", v );

    auto it = e->step(0)->beginNodal();
    auto en = e->step(0)->endNodal();
    BOOST_CHECK( std::distance( it,en ) == 2 );
    auto n = it->second.second[0][0]->name() ;
    BOOST_MESSAGE( "1st name : " << n );
    BOOST_CHECK_EQUAL(  n, "u" );
    n = (++it)->second.second[0][0]->name() ;
    BOOST_MESSAGE( "2nd name : " << n );
    BOOST_CHECK_EQUAL(  n, "v_v" );

    e->save();
}

BOOST_AUTO_TEST_SUITE_END()
