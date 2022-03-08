/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
   Date: 2011-06-19

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
/**
   \file test_form_interpolation.cpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-06-19
*/

// give a name to the testsuite
#define BOOST_TEST_MODULE form_interpolation testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>


using namespace Feel;
using namespace Feel::vf;

typedef Application Application_type;
typedef std::shared_ptr<Application_type> Application_ptrtype;


namespace test_form_interpolation
{

/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test form interpolation options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.08 ), "mesh size" )
    ( "hsize-bis", po::value<double>()->default_value( 0.04 ), "mesh bis size" )

    ;
    return desc_options.add( Feel::feel_options() );//.add( backend_options( "pressure" ) );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "test_form_interpolation" ,
                     "test_form_interpolation" ,
                     "0.1",
                     "test linear and bilinear form with interpolation",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Vincent Chabannes", "developer", "vincent.chabannes@imag.fr", "" );
    return about;
}


template <uint16_type OrderPoly>
void run()
{

    Environment::changeRepository( _directory=boost::format( "testsuite/feeldiscr/%1%/P%2%/" )
                                   % Environment::about().appName()
                                   % OrderPoly );


    /* change parameters below */
    const int nDim = 2;
    const int nOrderPoly = OrderPoly;
    const int nOrderGeo = 1;

    //------------------------------------------------------------------------------//

    typedef Mesh< Simplex<nDim, 1,nDim> > mesh_type;
    typedef Mesh< Simplex<nDim, 1,nDim> > mesh_bis_type;

    double meshSize = doption(_name="hsize");
    double meshSizeBis = doption(_name="hsize-bis");

    // mesh
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Rectangle Omega( meshSize,"Omega",x1,x2 );
    Omega.setMarker( _type="line",_name="Entree",_marker4=true );
    Omega.setMarker( _type="line",_name="Sortie",_marker2=true );
    Omega.setMarker( _type="line",_name="Paroi",_marker1=true,_marker3=true );
    Omega.setMarker( _type="surface",_name="Fluid",_markerAll=true );
    auto mesh = Omega.createMesh(_mesh=new mesh_type,_name="omega_"+ mesh_type::shape_type::name() );
    LOG(INFO) << "created mesh\n";

    GeoTool::Rectangle OmegaBis( meshSizeBis,"Omega",x1,x2 );
    OmegaBis.setMarker( _type="line",_name="Entree",_marker4=true );
    OmegaBis.setMarker( _type="line",_name="Sortie",_marker2=true );
    OmegaBis.setMarker( _type="line",_name="Paroi",_marker1=true,_marker3=true );
    OmegaBis.setMarker( _type="surface",_name="Fluid",_markerAll=true );
    auto meshBis = OmegaBis.createMesh(_mesh=new mesh_bis_type,_name="omegaBis_"+ mesh_type::shape_type::name() );

    //auto meshBis= mesh->createP1mesh();
    LOG(INFO) << "created meshBis\n";

    typedef Lagrange<nOrderPoly,Scalar,Continuous,PointSetFekete> basis_type;
    typedef FunctionSpace<mesh_type, bases<basis_type> > space_type;

    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();
    auto u2 = Xh->element();
    LOG(INFO) << "created space and elements\n";

    auto mybackend = backend();

    //--------------------------------------------------------------//

    auto A = mybackend->newMatrix( _test=Xh, _trial=Xh );
    auto F = mybackend->newVector( Xh );
    auto pi = cst( M_PI );
    auto u_exact = cos( pi*Px() )*sin( pi*Py() );
    auto dudx = -pi*sin( pi*Px() )*sin( pi*Py() );
    auto dudy = pi*cos( pi*Px() )*cos( pi*Py() );
    auto grad_u_uexact = vec( dudx,dudy );
    auto lap = -2*pi*pi*cos( pi*Px() )*sin( pi*Py() );
    //auto lap = -pi*pi*cos(pi*Px())*sin(pi*Py())
    //    -pi*pi*cos(pi*Px())*sin(pi*Py());
    LOG(INFO) << "created exact data and matrix/vector\n";

    auto f = -lap;//cst(1.);
    double gammabc=10;

    // assemblage
    form2( _test=Xh, _trial=Xh, _matrix=A, _init=true ) =
        integrate( _range=elements( mesh ), //_Q<15>(),
                   _expr=+ gradt( u )*trans( grad( v ) ) );

    form2( _test=Xh, _trial=Xh, _matrix=A ) +=
        integrate( _range=boundaryfaces( mesh ),
                   _expr=- gradt( u )*N()*id( v )
                   + gammabc*idt( u )*id( v )/hFace() );

    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( _range=elements( mesh ), // _Q<10>(),
                   _expr=trans( f )*id( v ) );

    form1( _test=Xh, _vector=F ) +=
        integrate( _range=boundaryfaces( mesh ),
                   _expr=+ gammabc*u_exact*id( v )/hFace() );

    LOG(INFO) << "A,F assembled\n";

    //form2( _test=Xh, _trial=Xh, _matrix=A ) +=
    //    on( boundaryfaces(mesh) , u, F, u_exact );

    // solve system
    mybackend->solve( _matrix=A,_solution=u,_rhs=F );
    LOG(INFO) << "A u = F solved\n";
    //--------------------------------------------------------------//

    auto a2 = form2( _test=Xh, _trial=Xh );
    auto f2 = form1( _test=Xh );
    LOG(INFO) << "created form2 a2 and form1 F2\n";

    // assemblage

    a2 = integrate( _range=elements( meshBis ),
                    _expr=+ gradt( u2 )*trans( grad( v ) ),
                    _quad=_Q<10>() );
    LOG(INFO) << "a2 grad.grad term\n";
    a2 += integrate( _range=boundaryfaces( meshBis ),
                     _expr=- gradt( u2 )*N()*id( v )
                     + gammabc*idt( u2 )*id( v )/hFace(),
                     _quad=_Q<10>() );
    LOG(INFO) << "a2 weak dirichlet terms\n";

    f2 = integrate( _range=elements( meshBis ),
                    _expr=trans( f )*id( v ),
                    _quad=_Q<10>() );
    LOG(INFO) << "f2 source term\n";
    f2 += integrate( _range=boundaryfaces( meshBis ),
                     _expr=+ gammabc*u_exact*id( v )/hFace(),
                     _quad=_Q<10>() );
    LOG(INFO) << "f2 dirichlet terms\n";

    LOG(INFO) << "a2,f2 assembled\n";

    //form2( _test=Xh, _trial=Xh, _matrix=A2 ) +=
    //     on( boundaryfaces(mesh) , u2, F2, u_exact );



#if 0

    for ( size_type i = 0 ; i< F->size() ; ++i )
    {
        auto errOnF = std::abs( ( *F )( i )-( *F2 )( i ) );

        if ( errOnF > 1e-8 )
            std::cout << "\nerrOnF : " << errOnF;
    }

    std::cout << "\nFin errOnF !!!!\n";
#endif

    // solve system
    a2.solve( _rhs=f2, _solution=u2 );


    auto diff = std::sqrt( integrate( _range=elements( mesh ), _expr=( idv( u )-idv( u2 ) )*( idv( u )-idv( u2 ) ) ).evaluate()( 0,0 ) );
#if 0
    auto int1 = integrate( _range=elements( mesh ), _expr=abs( idv( u ) ) ).evaluate()( 0,0 );
    auto int2 = integrate( _range=elements( mesh ), _expr=abs( idv( u2 ) ) ).evaluate()( 0,0 );

    std::cout << "\nThe diff : " << diff
              << " int1 :" << int1
              << " int2 :" << int2
              << "\n";
#endif
#if USE_BOOST_TEST
    BOOST_CHECK_SMALL( diff,1e-2 );
#endif


    //--------------------------------------------------------------//

    if ( boption(_name="exporter.export") )
    {
        // export
        auto ex = exporter( _mesh=mesh );
        ex->add( "u", u );
        ex->add( "ubis", u2 );
        ex->save();
    }
}

} // namespace test_form_interpolation

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_form_interpolation::makeAbout(),
                                 test_form_interpolation::makeOptions() )

BOOST_AUTO_TEST_SUITE( form_interpolation )

BOOST_AUTO_TEST_CASE( form_interpolation )
{

    test_form_interpolation::run<2>();
}
BOOST_AUTO_TEST_SUITE_END()
