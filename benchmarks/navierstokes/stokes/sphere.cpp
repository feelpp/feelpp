/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-07-19

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
/**
   \file sphere.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-07-19
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

template<typename MeshType>
void myexport( std::string const& name, boost::shared_ptr<MeshType> mesh )
{
    using namespace Feel;
    using namespace Feel::vf;

    typedef Exporter<MeshType> exporter_type;
    auto exporter = exporter_type::New( "gmsh", name );
    exporter->step( 0 )->setMesh( mesh );
    exporter->save();
}

int main( int argc, char** argv )
{
    double hsphere = 1;
    double hcube = 1;
    double R = 1;
    double U = 1;
    double mu = 1;
    bool straighten = true;
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc( "Allowed options" );
    desc.add_options()
    ( "help", "produce help message" )
    ( "geomap", po::value<int>()->default_value( 0 ), "straighten mesh (0=ho,1=p1,2=opt)" )
    ( "straighten", po::value<bool>( &straighten )->default_value( true ), "straighten mesh" )
    ( "hcube", po::value<double>( &hcube )->default_value( 1 ), "h size for cube" )
    ( "hsphere", po::value<double>( &hsphere )->default_value( 1 ), "h size for sphere" )
    ( "R", po::value<double>( &R )->default_value( 1 ), "sphere radius" )
    ( "U", po::value<double>( &U )->default_value( 1 ), "characteristic fluid velocity" )
    ( "mu", po::value<double>( &U )->default_value( 1 ), "fluid viscosity" )
    ;

    po::variables_map vm;
    po::store( po::parse_command_line( argc, argv, desc ), vm );
    po::notify( vm );

    const int Order=3;
    const int OrderP=3;
    using namespace Feel;
    using namespace Feel::vf;
    Feel::Environment env( argc, argv );
    typedef Mesh<Simplex<3,Order> > mesh_type;

    GeoTool::Cube R1( hcube,"R1",GeoTool::Node( -3,-3,-3 ),GeoTool::Node( 3,3,3 ) );
    R1.setMarker( _type="surface",_name="Cube",_markerAll=true );
    R1.setMarker( _type="volume",_name="Omega",_markerAll=true );

    GeoTool::Sphere C1( hsphere,"C1",GeoTool::Node( 0,0,0 ),GeoTool::Node( R,0,0 ) );
    C1.setMarker( _type="surface",_name="Sphere",_markerAll=true );

    auto R1mC1mesh = ( R1-C1 ).createMesh( _mesh=new mesh_type, _name = ( boost::format( "sphere-%1%-%2%-%3%" )% Order % hcube % hsphere ).str() );

    //myexport<mesh_type>("R1-C1",R1mC1mesh);

    auto mesh = R1mC1mesh;

    if ( straighten )
        mesh = straightenMesh( R1mC1mesh, Environment::worldComm()  );


    typedef FunctionSpace<mesh_type,bases<Lagrange<Order,Vectorial> > > Vh_t;
    typedef FunctionSpace<mesh_type,bases<Lagrange<Order,Scalar> > > Ph_t;
    auto Vh = Vh_t::New( mesh );
    auto Ph = Ph_t::New( mesh );
    typedef FunctionSpace<mesh_type,bases<Lagrange<OrderP,Vectorial> > > VhP_t;
    typedef FunctionSpace<mesh_type,bases<Lagrange<OrderP,Scalar> > > PhP_t;
    auto VhP = VhP_t::New( mesh );
    auto PhP = PhP_t::New( mesh );

    auto r = sqrt( Px()*Px()+Py()*Py()+Pz()*Pz() );
    auto rxy = sqrt( Px()*Px()+Py()*Py() );
    auto sintheta = sin( acos( Pz()/r ) );
    auto costheta = Pz()/r;

    auto er=vec( Px()/r,Py()/r,Pz()/r );
    auto etheta=vec( Px()*Pz()/( r*rxy ),Py()*Pz()/( r*rxy ),-rxy/r );
    auto ephi=vec( -Py()/rxy,Px()/rxy,cst( 0. ) );

    auto phi = vf::project( _space=PhP, _expr=-R*2*U*0.5*pow( ( r/R )*sintheta,2 )*( 1-1.5*R/r+0.5*pow( R/r,3 ) ) );
    auto w = vf::project( _space=VhP, _expr=1.5*( U/R )*pow( R/r,2 )*sintheta*ephi );
    auto ur = -2*U*( 1-1.5*( R/r )+.5*pow( R/r,3 ) )*costheta/R;
#if 0
    auto utheta = ( 2*U*r*sintheta*sintheta*( 1-1.5*( R/r )+.5*pow( R/r,3 ) ) +
                    U*( r*sintheta )*( r*sintheta )*( 1.5*( R/( r*r ) )-1.5*pow( R,3 )/( r*r*r*r ) ) )/( R*r*sintheta );
#else
    auto utheta = ( 2*U*sintheta*( 1-1.5*( R/r )+.5*pow( R/r,3 ) ) +
                    U*( r*sintheta )*( 1.5*( R/( r*r ) )-1.5*pow( R,3 )/( r*r*r*r ) ) )/( R );
#endif
    auto u = vf::project( _space=VhP, _expr=ur*er+utheta*etheta );
    auto p = vf::project( _space=PhP, _expr=1.5*( mu*U/R )*( R/r )*( R/r )*costheta );
    auto st = vf::project( _space=PhP, _expr=sintheta );
    auto uur = vf::project( _space=PhP, _expr=ur );
    auto uut = vf::project( _space=PhP, _expr=utheta );

    auto u_1 = U*( ( 3*R*Px()*Px()/( 4*r*r*r ) )*( R*R/( r*r ) -1 ) -( R/4*r )*( 3 + R*R/( r*r ) ) +1 );
    auto u_2 = U*( ( 3*R*Px()*Py()/( 4*r*r*r ) )*( R*R/( r*r ) -1 ) );
    auto u_3 = U*( ( 3*R*Px()*Pz()/( 4*r*r*r ) )*( R*R/( r*r ) -1 ) );
    auto p_goncalo = vf::project( _space=PhP, _expr=3*mu*R*U*Px()/( 2*r*r*r ) );
    auto u_goncalo = vf::project( _space=VhP, _expr=vec( u_1,u_2,u_3 ) );

    Eigen::Vector3d force_exact_x;
    force_exact_x << 6*M_PI*mu*U*R, 0, 0;
    std::cout << "force_exacte = " << force_exact_x << "\n";

    auto sigmaN_goncalo=-idv( p_goncalo )*N()+mu*( gradv( u_goncalo ) )*N();
    boost::timer ti;
    auto force_goncalo = integrate( _range=markedfaces( mesh,"Sphere" ), _expr=sigmaN_goncalo,_quad=_Q<15>()  ).evaluate();
    std::cout << "force_goncalo = " << force_goncalo << "\n";
    std::cout << "error_goncalo=" << ( force_exact_x-force_goncalo ).norm() << " time: " << ti.elapsed() << "\n";


    auto sigmaN=-idv( p )*N()+mu*( gradv( u ) )*N();
    auto volume = integrate( _range=elements( mesh ), _expr=cst( 1.0 ),_quad=_Q<6>(), _geomap=( GeomapStrategyType )vm["geomap"].as<int>()  ).evaluate();
    auto surface = integrate( _range=markedfaces( mesh,"Sphere" ), _expr=cst( 1.0 ),_quad=_Q<6>() ).evaluate();
    auto force = integrate( _range=markedfaces( mesh,"Sphere" ), _expr=sigmaN,_quad=_Q<15>() ).evaluate();

    std::cout << "volume = " << volume << "\n";
    double ve = 6*6*6-4*M_PI*R*R*R/3;
    std::cout << "volume exact = " << ve << "\n";
    std::cout << "error volume = " << std::scientific << math::abs( volume( 0,0 )-ve ) << "\n";

    std::cout << "surface = " << surface << "\n";
    std::cout << "surface exact = " << 4*M_PI*R*R << "\n";
    std::cout << "error surface = " << math::abs( surface( 0,0 )-4*M_PI*R*R ) << "\n";
    Eigen::Vector3d force_exact_z;
    force_exact_z << 0,0,6*M_PI*mu*U*R;
    std::cout << "force = " << force << "\n";
    std::cout << "error=" << ( force-force_exact_z ).norm() << "\n";

    auto v_phi = vf::project( _space=Ph, _expr=idv( phi ) );
    auto v_w = vf::project( _space=Vh, _expr=idv( w ) );
    auto v_u = vf::project( _space=Vh, _expr=idv( u ) );
    auto v_st = vf::project( _space=Ph, _expr=idv( st ) );
    auto v_ur = vf::project( _space=Ph, _expr=idv( uur ) );
    auto v_ut = vf::project( _space=Ph, _expr=idv( uut ) );
    auto v_u_goncalo = vf::project( _space=Vh, _expr=idv( u_goncalo ) );
    auto v_p = vf::project( _space=Ph, _expr=idv( p ) );
    auto v_p_goncalo = vf::project( _space=Ph, _expr=idv( p_goncalo ) );

    typedef Exporter<mesh_type,Order> exporter_type;
    auto exporter = exporter_type::New( "gmsh", "sphere" );
    exporter->step( 0 )->setMesh( mesh );
    exporter->step( 0 )->add( "phi", v_phi );
    exporter->step( 0 )->add( "w", v_w );
    exporter->step( 0 )->add( "u", v_u );
    exporter->step( 0 )->add( "sin(theta)", v_st );
    exporter->step( 0 )->add( "ur", v_ur );
    exporter->step( 0 )->add( "utheta", v_ut );
    exporter->step( 0 )->add( "ux", v_u.comp( X ) );
    exporter->step( 0 )->add( "uy", v_u.comp( Y ) );
    exporter->step( 0 )->add( "uy", v_u.comp( Z ) );
    exporter->step( 0 )->add( "p", v_p );
    exporter->step( 0 )->add( "u_1", v_u_goncalo );
    exporter->step( 0 )->add( "u_1x", v_u_goncalo.comp( X ) );
    exporter->step( 0 )->add( "u_2y", v_u_goncalo.comp( Y ) );
    exporter->step( 0 )->add( "u_3y", v_u_goncalo.comp( Z ) );
    exporter->step( 0 )->add( "p_goncalo", v_p_goncalo );

    exporter->save();

}


