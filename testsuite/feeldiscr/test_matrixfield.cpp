/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 18 Apr 2015

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
#if 0
#define BOOST_TEST_MODULE test_matrixfield
#include <testsuite.hpp>
#endif

#include <feel/options.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>


#include <iostream>

using namespace Feel;

namespace test_matrixfield
{

typedef Application Application_type;
typedef std::shared_ptr<Application_type> Application_ptrtype;

/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_matrixfield options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ;
    return desc_options.add( Feel::feel_options() );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "Test_Matrixfield" ,
                     "Test_Matrixfield" ,
                     "0.1",
                     "test matrix fields",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}


}

#if 0
FEELPP_ENVIRONMENT_WITH_OPTIONS( test_matrixfield::makeAbout(), test_matrixfield::makeOptions() )

BOOST_AUTO_TEST_SUITE( interp_matrixfield )

BOOST_AUTO_TEST_CASE( interp_matrixfield )
{

    using namespace test_matrixfield;


    auto test_app = Application_ptrtype( new Application_type );

    test_app->changeRepository( boost::format( "testsuite/feeldiscr/%1%/" )
                                % test_app->about().appName()
                              );

    typedef Mesh<Simplex<2,1,2> > mesh_type;
    typedef std::shared_ptr<  mesh_type > mesh_ptrtype;

    double meshSize = test_app->vm()["hsize"].as<double>();

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );

    auto mesh = R.createMesh(_mesh=new mesh_type,_name= "domain" );
    auto Mh = Pchm<1>( mesh) ;

}

BOOST_AUTO_TEST_SUITE_END()

#else

int main(int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="test_matrixfield",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    typedef Mesh<Simplex<2,1,2> > mesh_type;
    typedef std::shared_ptr<  mesh_type > mesh_ptrtype;

    double meshSize = doption("gmsh.hsize");

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );

    auto mesh = R.createMesh(_mesh=new mesh_type,_name= "domain" );
    auto Mh = Pchm<1>( mesh) ;
    auto v = Mh->element();

    //
    v.setOnes();
    auto M0 = mat<2,2>( cst(1.), cst(1.), cst(1.), cst(1.)  );
    auto i0 = integrate( _range=elements(mesh), _expr=idv(v) ).evaluate();
    auto i0_res = integrate( _range=elements(mesh), _expr=M0 ).evaluate();
    std::cout << "i0 = " << i0 << std::endl;
    std::cout << "i0_res = " << i0_res << std::endl;
    CHECK( (i0-i0_res).norm() < 1e-13 ) << "Invalid integral i0: expected " << i0_res << " got " << i0;

    auto M = mat<2,2>( cst(1.), cst(2.), cst(3.), cst(4.)  );
    v.on(_range=elements(mesh), _expr=M ) ;
    v.printMatlab("v.m");
    auto i1 = integrate( _range=elements(mesh), _expr=idv(v) ).evaluate();
    auto i1_res = integrate( _range=elements(mesh), _expr=M ).evaluate();
    std::cout << "i1 = " << i1 << std::endl;
    std::cout << "i1_res = " << i1_res << std::endl;
    CHECK( (i1-i1_res).norm() < 1e-13 ) << "Invalid integral: expected " << i1_res << " got " << i1;
}

#endif
