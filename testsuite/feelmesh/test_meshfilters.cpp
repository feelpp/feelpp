/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: April-05-2022

 Copyright (C) 2022 Feel++ Consortium

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
#define BOOST_TEST_MODULE test_meshfilters
#include <feel/feelcore/testsuite.hpp>


#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feelvf/measure.hpp>
#include <feel/feelmesh/ranges.hpp>

/** use Feel namespace */
using namespace Feel;

inline po::options_description makeOptions()
{
    po::options_description options( "Test Update Marker  Options" );
    return options;
}

inline AboutData
makeAbout()
{
    AboutData about( "test_meshfilters",
                     "test_meshfilters",
                     "0.2",
                     "nD(n=2,3) test meshfilters",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( meshfilters_suite )

//using dim_t = boost::mpl::list<boost::mpl::int_<2>, boost::mpl::int_<3>>;
using dim_t = boost::mpl::list<boost::mpl::int_<2>>;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_elements, T, dim_t )
{
    BOOST_MESSAGE( "test_elements starts for dim=" << T::value );
    auto mesh = unitHypercube<T::value>();

    int flag = 1;
    if ( T::value == 3 )
        flag = 15;
    auto len1 = measure( _range = elements( mesh ) );
    BOOST_CHECK_CLOSE( len1, 1, 1e-12 );

    BOOST_MESSAGE( "test_elements ends for dim=" << T::value );
}

BOOST_AUTO_TEST_SUITE_END()
