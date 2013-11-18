/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake1@e.ujf-grenoble.fr>
       Date: 2011-08-20

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

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
/**
   \file test_trace.cpp
   \author Abdoulaye Samake <abdoulaye.samake.@imag.fr>
   \date 2011-08-20
 */

#define BOOST_TEST_MODULE test_hypercubeinterpolation
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatortrace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>

namespace test_hypercubeinterpolation
{
using namespace Feel;
using namespace Feel::vf;

typedef Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;


inline
po::options_description
makeOptions()
{
    po::options_description testoptions( "Test options" );
    testoptions.add_options()
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ;
    return testoptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_hypercubeinterpolation" ,
                     "test_hypercubeinterpolation" ,
                     "0.2",
                     "nD(n=1,2,3) test operatorinterpolation for hypercube domain",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Abdoulaye Samake", "developer", "abdoulaye.samake@imag.fr", "" );
    return about;

}

template<int Dim, int Order>
void
run( Application_ptrtype & theApp )
{

    typedef Hypercube<Dim,1,Dim> convex_type;
    typedef Mesh< convex_type > mesh_type;
    typedef FunctionSpace<mesh_type, bases<Lagrange<Order,Scalar> > > space_type;

    double meshSize = theApp->vm()["hsize"].as<double>();

    if ( !theApp->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/P%4%G%4%/h_%5%/" )
                                       % theApp->about().appName()
                                       % convex_type::name()
                                       % Dim
                                       % Order
                                       % meshSize );

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name=( boost::format( "hypercube-%1%" ) % Dim ).str(),
                                              _addmidpoint=false,
                                              _usenames=false,
                                              _shape="hypercube",
                                              _dim=Dim,
                                              _h=meshSize,
                                              _convex="Hypercube",
                                              //_convex=convex_type::name(),
                                              _xmin=0.,
                                              _xmax=1.,
                                              _ymin=0.,
                                              _ymax=1.,
                                              _substructuring=true
                                              ),
                                _structured=2);

    auto mesh2 = createGMSHMesh( _mesh=new mesh_type,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                 _desc=domain( _name=( boost::format( "hypercube2-%1%" ) % Dim ).str(),
                                               _addmidpoint=false,
                                               _usenames=false,
                                               _shape="hypercube",
                                               _dim=Dim,
                                               _h=0.075,
                                               _convex="Hypercube",
                                               //_convex=convex_type::name(),
                                               _xmin=0.,
                                               _xmax=1.,
                                               _ymin=0.,
                                               _ymax=1.,
                                               _substructuring=true
                                               ),
                                 _structured=2);


    auto backend = backend_type::build( theApp->vm() );
    auto Xh = space_type::New(_mesh=mesh);
    auto TXh = Xh->trace( markedfaces( mesh,"NORTH" ) ) ;

    auto opI = opInterpolation( _domainSpace=Xh,
                                _imageSpace=TXh,
                                _range=elements(TXh->mesh()),
                                _backend=backend
                                );
    opI->matPtr()->printMatlab("opIM.m");


    auto u = vf::project( _space=Xh,
                          _range=markedfaces(mesh,"NORTH"),
                          //_expr=cos( pi*Px() )*sin( pi*Py() )
                          _expr=cst(1.) );

    auto tu = TXh->element();
    opI->apply( u,tu );
    u.printMatlab( "u1.m" );
    tu.printMatlab( "tu1.m" );

    auto l2error = normL2( _range=elements( TXh->mesh() ), _expr=idv( u )-idv( tu) );

    // BOOST_TEST_MESSAGE( "l2 error 1 = " << l2error );
    // BOOST_CHECK_SMALL( l2error, 1e-13 );
    std::cout<< "l2 error 1 = " << l2error <<"\n";

    u = vf::project( _space=Xh,
                     _range=markedfaces(mesh,"NORTH"),
                     _expr=sin( pi*Px() )*cos( pi*Py() ) );
    opI->apply( u,tu );
    u.printMatlab( "u2.m" );
    tu.printMatlab( "tu2.m" );

    l2error = normL2( _range=elements( TXh->mesh() ), _expr=idv( u )-idv( tu) );

    //BOOST_TEST_MESSAGE( "l2 error 2 = " << l2error );
    //BOOST_CHECK_SMALL( l2error, 1e-13 );
    std::cout<< "l2 error 2 = " << l2error <<"\n";

    auto Xh2 = space_type::New(_mesh=mesh2);
    auto opI2 = opInterpolation( _domainSpace=Xh,
                                 _imageSpace=Xh2,
                                 _range=elements(Xh2->mesh()),
                                 _backend=backend
                                 //_type=InterpolationConforme()
                                 );
    opI2->matPtr()->printMatlab("opIM2.m");

    u = vf::project( _space=Xh,
                     _range=elements(mesh),
                     _expr=cos( pi*Px() )*sin( pi*Py() ) );
                     //_expr=cst(1.) );

    auto u2 = Xh2->element();
    opI2->apply( u,u2 );
    u.printMatlab( "u2.m" );
    u2.printMatlab( "u3.m" );

    l2error = normL2( _range=elements( Xh2->mesh() ), _expr=idv( u )-idv( u2 ) );

    // BOOST_TEST_MESSAGE( "l2 error 3 = " << l2error );
    // BOOST_CHECK_SMALL( l2error, 1e-13 );
    std::cout<< "l2 error 3 = " << l2error <<"\n";

    auto opI3 = opInterpolation( _domainSpace=Xh2,
                                 _imageSpace=Xh,
                                 _range=elements(Xh->mesh()),
                                 _backend=backend
                                 //_type=InterpolationConforme()
                                 );
    opI3->matPtr()->printMatlab("opIM3.m");

    u.zero();
    u2.zero();

    u2 = vf::project( _space=Xh2,
                      _range=elements(mesh2),
                      _expr=cos( pi*Px() )*sin( pi*Py() ) );
    //_expr=cst(1.) );

    opI3->apply( u2,u );
    u.printMatlab( "u4.m" );
    u2.printMatlab( "u5.m" );

    //l2error = normL2( _range=elements( Xh2->mesh() ), _expr=idv( u )-idv( u2 ) );
    double s_1 = normL2( _range=elements( Xh2->mesh() ), _expr=idv( u )*idv( u ) );
    double s_2 = normL2( _range=elements( Xh2->mesh() ), _expr=idv( u2 )*idv( u2 ) );

    //BOOST_TEST_MESSAGE( "l2 error 4 = " << l2error );
    //BOOST_CHECK_SMALL( l2error, 1e-13 );
    //std::cout<< "l2 error 4 = " << l2error <<"\n";
    std::cout<< "s1 = " << s_1 <<"\n";
    std::cout<< "s2 = " << s_2 <<"\n";

} // Test::run

} // test_hypercubeinterpolation

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_hypercubeinterpolation::makeAbout(),
                                 test_hypercubeinterpolation::makeOptions() )

BOOST_AUTO_TEST_SUITE( hypercubeinterpolation )

typedef Feel::Application Application_type;
typedef boost::shared_ptr<Application_type> Application_ptrtype;


BOOST_AUTO_TEST_CASE( hypercubeinterpolation1 )
{
    auto theApp = Application_ptrtype( new Application_type );
    test_hypercubeinterpolation::run<2,1>( theApp );

}
BOOST_AUTO_TEST_SUITE_END()


//#if 0
// int
// main( int argc, char** argv )
// {
//     typedef Feel::Application Application_type;
//     typedef boost::shared_ptr<Application_type> Application_ptrtype;

//     auto theApp = Application_ptrtype( new Application_type( argc,argv,
//                                                              test_hypercubeinterpolation::makeAbout(),
//                                                              test_hypercubeinterpolation::makeOptions() ) );

//     // if ( theApp->vm().count( "help" ) )
//     // {
//     //     std::cout << theApp->optionsDescription() << "\n";
//     //     exit( 0 );
//     // }

//     test_hypercubeinterpolation::run<2,1>( theApp );

// }
//#endif

