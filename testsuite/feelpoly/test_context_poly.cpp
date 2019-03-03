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
//! @date 02 Mar 2019
//! @copyright 2019 Feel++ Consortium
//!


#define BOOST_TEST_MODULE test_context_poly
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/context.hpp>

/** use Feel namespace */
using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "test_context_poly" ,
                     "test_context_poly" ,
                     "0.2",
                     "nD(n=2,3) test inner",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( context_suite )

BOOST_AUTO_TEST_CASE( t0 )
{
    using namespace vm;
    const size_type v = POINT|JACOBIAN;
    BOOST_CHECK_EQUAL((has_value_v<v,POINT>),true);
    BOOST_CHECK_EQUAL((has_point_v<v,POINT>),true);
    BOOST_CHECK_EQUAL((has_value_v<v,JACOBIAN>),true);
    BOOST_CHECK_EQUAL((has_value_v<v,JACOBIAN|POINT>),true);
    BOOST_CHECK_EQUAL((has_value_v<v,HESSIAN>),false);

}

BOOST_AUTO_TEST_CASE( test_1 )
{

}

BOOST_AUTO_TEST_SUITE_END()





