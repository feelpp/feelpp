/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2011-07-09

   Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define USE_BOOST_TEST 1

#define BOOST_TEST_MODULE integration_opt testsuite

#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>


using namespace Feel;
using namespace Feel::vf;

double hsize = 2;
int straighten = 1;
int nlevels = 3;

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_integration_opt" ,
                           "test_integration_opt" ,
                           "0.1",
                           "integration tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2006-2010 Universite Joseph Fourier (Grenoble I)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

po::options_description
makeOptions()
{
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc( "Allowed options" );
    desc.add_options()
    ( "help", "produce help message" )
    ( "hsize", po::value<double>( &hsize )->default_value( 2 ), "h size" )
    ( "straight", po::value<int>( &straighten )->default_value( 1 ), "straighten" )
    ( "nlevels", po::value<int>( &nlevels )->default_value( 3 ), "number of refinement levels" )
    ( "shape", po::value<std::string>()->default_value( "pie" ), "pie,circle" )
    ;
    return desc.add( Feel::feel_options() );
}
#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( integration_opt )

typedef boost::mpl::list<boost::mpl::pair<mpl::int_<2>,mpl::int_<2> >,
        boost::mpl::pair<mpl::int_<2>,mpl::int_<3> >,
        boost::mpl::pair<mpl::int_<2>,mpl::int_<4> >,
        boost::mpl::pair<mpl::int_<3>,mpl::int_<2> >,
        boost::mpl::pair<mpl::int_<3>,mpl::int_<3> >
        > dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( integration_opt, T, dim_types )
{
    BOOST_TEST_MESSAGE( "============================================================\n"
                        << "Dim: " << T::first::value << "  Order: " << T::second::value << "\n" );

    using namespace Feel;
    using namespace Feel::vf;

    typedef Mesh<Simplex<T::first::value,T::second::value> > mesh_type;
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name=( boost::format( "ellipsoid-%1%-%2%" ) % T::first::value % T::second::value ).str() ,
                                        _usenames=true,
                                        _shape="ellipsoid",
                                        _dim=T::first::value,
                                        _order=T::second::value,
                                        _h=hsize ),
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                                _straighten=straighten );
    saveGMSHMesh( _mesh=mesh, _filename=( boost::format( "ellipsoid-%1%-%2%-%3%-saved.msh" ) % T::first::value % T::second::value % straighten ).str() );

    boost::timer ti;
    auto i1 = integrate( _range=elements( mesh ), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_HO ).evaluate().norm();
    std::cout << "Ho: " << ti.elapsed() << "s\n";
    ti.restart();
    auto i2 = integrate( _range=elements( mesh ), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_OPT ).evaluate().norm();
    std::cout << "Opt: " << ti.elapsed() << "s\n";
    ti.restart();
    auto i3 = integrate( _range=elements( mesh ),  _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_O1 ).evaluate().norm();
    std::cout << "P1: " << ti.elapsed() << "s\n";
    ti.restart();

    BOOST_CHECK_CLOSE( i1, i2, 1e-12 );
    BOOST_TEST_MESSAGE( "ho = " << std::scientific << std::setprecision( 16 ) << i1 << "\n" <<
                        "opt = " << std::scientific << std::setprecision( 16 ) << i2 << "\n"
                        "p1 = " << std::scientific << std::setprecision( 16 ) << i3 << "\n"
                        << "\n" );


}

//typedef boost::mpl::list<mpl::int_<1>,mpl::int_<2>, mpl::int_<3>, mpl::int_<4>, mpl::int_<5> > order_types;
typedef boost::mpl::list<mpl::int_<1>,mpl::int_<2>, mpl::int_<3>, mpl::int_<4>> order_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_pie, T, order_types )
{
    BOOST_TEST_MESSAGE( "============================================================\n"
                        << "Order: " << T::value << "\n" );


    po::variables_map vm = Environment::vm();

    using namespace Feel;
    using namespace Feel::vf;
    typedef Mesh<Simplex<2,T::value> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    mesh_ptrtype mesh;
    std::string shape =  vm["shape"].template as<std::string>();
    std::map<std::string,std::vector<boost::tuple<double,double,double,double> > > ho;
    double exact=1;

    for ( int l = 1; l <= nlevels; ++l )
    {
        double meshSize= hsize/std::pow( 2.,l );

        if ( shape == "pie" )

        {
            exact = M_PI/8;
            GeoTool::Node C( 0,0 );
            GeoTool::Node A( 1,0. );
            GeoTool::Node B( std::sqrt( 2. )/2,std::sqrt( 2. )/2 ); // 45D

            GeoTool::Pie Pie( meshSize,"pie", C, A, B );
            Pie.setMarker(_type="line",_name="Boundary",_markerAll=true);
            Pie.setMarker(_type="surface",_name="Omega",_markerAll=true);

            mesh = Pie.createMesh(_mesh=new mesh_type,
                                  _name=( boost::format( "pie-%1%-%2%" ) % T::value % l ).str() );
        }

        if ( shape == "circle" )
        {
            exact = M_PI;
            GeoTool::Node C( 0,0 );
            GeoTool::Node A( 1,0. );
            GeoTool::Circle Circle( meshSize,"circle", C, A );
            Circle.setMarker(_type="line",_name="Boundary",_markerAll=true);
            Circle.setMarker(_type="surface",_name="Omega",_markerAll=true);

            mesh = Circle.createMesh(_mesh=new mesh_type,
                                     _name=( boost::format( "circle-%1%-%2%" ) % T::value % l ).str() );
        }


        boost::timer ti;
        double i1 =  integrate( _range=elements( mesh ), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_HO ).evaluate().norm();
        ho["ho"].push_back( boost::make_tuple( meshSize, i1, math::abs( i1-exact ), ti.elapsed() ) );
        ti.restart();
        double i2 = integrate( _range=elements( mesh ), _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_OPT ).evaluate().norm();
        ho["opt"].push_back( boost::make_tuple( meshSize, i2, math::abs( i2-exact ), ti.elapsed() ) );
        ti.restart();
        BOOST_CHECK_CLOSE( i1, i2, 1e-12 );
        //std::cout << "Opt: " << ti.elapsed() << "s\n";
        double i3 = integrate( _range=elements( mesh ),  _expr=cst( 1. ), _geomap=GeomapStrategyType::GEOMAP_O1 ).evaluate().norm();
        ho["p1"].push_back( boost::make_tuple( meshSize, i3, math::abs( i3-exact ), ti.elapsed() ) );
        //std::cout << "P1: " << ti.elapsed() << "s\n";
    }


    for ( auto it = ho.begin(); it != ho.end(); ++ it )
    {
        BOOST_TEST_MESSAGE( std::setw( 10 ) << std::right << "levels" <<
                            std::setw( 10 ) << std::right << "h"  <<
                            std::setw( 15 ) << std::right << it->first <<
                            std::setw( 15 ) << std::right << "error" <<
                            std::setw( 15 ) << std::right << "ROC" <<
                            std::setw( 10 ) << std::right << "Time(s)"
                            << "\n" );

        for ( int l = 1; l <= nlevels; ++l )
        {
            auto data = it->second;
            double roc = 1;

            if ( l > 1 )
                roc = std::log10( data[l-2].template get<2>()/data[l-1].template get<2>() )/std::log10( data[l-2].template get<0>()/data[l-1].template get<0>() );

            BOOST_TEST_MESSAGE( std::right << std::setw( 10 ) << l <<
                                std::right << std::setw( 10 ) << data[l-1].template get<0>() <<
                                std::right << std::setw( 15 ) << std::scientific << std::setprecision( 5 ) << data[l-1].template get<1>() <<
                                std::right << std::setw( 15 ) << std::scientific << std::setprecision( 5 ) << data[l-1].template get<2>() <<
                                std::right << std::setw( 15 ) << std::scientific << std::setprecision( 5 ) << roc <<
                                std::right << std::setw( 10 ) << std::scientific << std::setprecision( 2 ) << data[l-1].template get<3>() << "\n" );
        }
    }

}
BOOST_AUTO_TEST_SUITE_END()

#else

int
main( int argc, char** argv )
{

}

#endif
