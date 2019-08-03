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
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelcrb/sensordesc.hpp>
#include <feel/feelcrb/sensors.hpp>

/** use Feel namespace */
using namespace Feel;

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options("test_sensor options");
    options.add_options()
        ("sensor.filename", Feel::po::value<std::string>()->default_value( "$cfgdir/sensordesc.csv" ), "file describing sensor network")
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

    //SensorDescriptionMap<3> desc( Environment::expand( soption( "sensor.filename" ) ) );
    SensorDescriptionMap<3> desc( Environment::expand( "$cfgdir/sensordesc.csv" ) );
    // verify that those keys / sensors exist
    BOOST_CHECK_EQUAL( desc.count( "zigduino-1" ), 1 );
    BOOST_CHECK_EQUAL( desc.count( "zigduino-2" ), 1 );
    BOOST_CHECK_EQUAL( desc.count( "zigduino-12" ), 0 );

    BOOST_CHECK_EQUAL( desc.at( "zigduino-1" ).type(), "gaussian" );
}
BOOST_AUTO_TEST_CASE( t1 )
{
    using namespace Feel;
    using namespace Feel::vf;
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3>> );
    auto Vh = Pch<1>( mesh );

    //SensorDescriptionMap<3> desc( Environment::expand( soption( "sensor.filename" ) ) );
    SensorDescriptionMap<3> desc( Environment::expand( "$cfgdir/sensordesc_bis.csv" ) );
    SensorMap<Pch_type<Mesh<Simplex<3>>,1>> sensors( Vh, desc );

    auto v = Vh->element();
    v.on(_range=elements(mesh), _expr=cst(1.));
    if ( sensors.count("zigduino-35") )
    {
        double zig10_v = sensors.at("zigduino-35")->operator()(v);
        BOOST_TEST_MESSAGE( "v : " << zig10_v );
    }
        
    // BOOST_CHECK_CLOSE( zig10_v, zig10_v_exact, 1e-10 );
    // if close to 0, use _SMALL
    // BOOST_CHECK_SMALL( zig10_v, 1e-10 );

}

BOOST_AUTO_TEST_SUITE_END()

