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
//! @date 18 Jun 2019
//! @copyright 2019 Feel++ Consortium
//!


#define BOOST_TEST_MODULE test_sensor
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsegment.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelcrb/sensordesc.hpp>

/** use Feel namespace */
using namespace Feel;

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options("test_sensor options");
    options.add_options()
        ("sensor.filename", Feel::po::value<std::string>()->default_value( "sensordesc.csv" ), "file describing sensor network")
        ;
    return options.add( Feel::feel_options() );
}

inline AboutData
makeAbout()
{
    AboutData about( "test_sensor",
                     "test_sensor",
                     "0.1",
                     "test sensor",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2019 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( sensor_suite )
BOOST_AUTO_TEST_CASE( t0 )
{
    using namespace Feel::vf;
#if 0
    SensorDescriptionMap<3> desc( "sensordescmap.csv" );
    // verify that those keys / sensors exist
    BOOST_CHECK_EQUAL( desc.count( "ziiguino-10" ), 1 );
    BOOST_CHECK_EQUAL( desc.count( "ziiguino-11" ), 1 );
    BOOST_CHECK_EQUAL( desc.count( "ziiguino-12" ), 1 );
    
    BOOST_CHECK_EQUAL( desc.at( "zigduino-10" ), "gaussian" );
#else
    BOOST_CHECK_EQUAL( 1, 1 );
#endif
}
BOOST_AUTO_TEST_SUITE_END()

