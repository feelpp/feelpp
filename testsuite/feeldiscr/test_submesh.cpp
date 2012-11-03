/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-10-21

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file test_submesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-10-21
 */
#define USE_BOOST_TEST 1
// Boost.Test

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE submesh testsuite
#include <testsuite/testsuite.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelvf/vf.hpp>

const double DEFAULT_MESH_SIZE=0.1;

namespace Feel
{
template<typename T, int Dim, int Order = 1>
struct imesh
{
    typedef Simplex<Dim, Order> convex_type;
    typedef Mesh<convex_type, T > type;
    typedef boost::shared_ptr<type> ptrtype;
};

template<typename value_type = double, int Dim=2>
struct test_submesh: public Application
{
    typedef typename imesh<value_type,Dim>::convex_type convex_type;
    typedef typename imesh<value_type,Dim>::type mesh_type;
    typedef typename imesh<value_type,Dim>::ptrtype mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<3, Scalar> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename mesh_type::location_element_const_iterator location_element_const_iterator;
    typedef Backend<value_type> backend_type;

    test_submesh()
        :
        Application(),
        backend( Backend<double>::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() ),
        mesh()
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                             _usenames=true,
                                             _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                             _shape=shape,
                                             _dim=Dim,
                                             _h=meshSize ),
                               _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    }
    void operator()()
    {
        location_element_const_iterator it,en;
        boost::tie( it,en ) = mesh->boundaryElements( 0 );
        mesh_ptrtype meshbdy( new mesh_type );
        mesh->createSubmesh( *meshbdy, it, en );
        BOOST_CHECK_EQUAL( meshbdy->numElements(), std::distance( it, en ) );
        //saveGMSHMesh( _mesh=meshbdy, _filename=shape+"_sub.msh" );
        using namespace Feel::vf;
        double intm1 = integrate( elements( meshbdy ), cst( 1. ) ).evaluate()( 0,0 );
        double intm2 = integrate( boundaryelements( mesh ), cst( 1. ) ).evaluate()( 0,0 );
        BOOST_CHECK_CLOSE( intm1, intm2, 1e-12 );
        //double intm11 = integrate( boundaryfaces(meshbdy), cst(1.) ).evaluate()(0,0);
        //double intm21 = integrate( boundaryfaces(mesh), cst(1.) ).evaluate()(0,0);
        //BOOST_CHECK_CLOSE( intm11, intm21, 1e-12 );

        double intm12 = integrate( markedfaces( meshbdy,"Dirichlet" ), cst( 1. ) ).evaluate()( 0,0 );
        double intm22 = integrate( markedfaces( mesh,"Dirichlet" ), cst( 1. ) ).evaluate()( 0,0 );
        BOOST_CHECK_CLOSE( intm12, intm22, 1e-12 );
        double intm13 = integrate( markedfaces( meshbdy,"Neumann" ), cst( 1. ) ).evaluate()( 0,0 );
        double intm23 = integrate( markedfaces( mesh,"Neumann" ), cst( 1. ) ).evaluate()( 0,0 );
        BOOST_CHECK_CLOSE( intm13, intm23, 1e-12 );


        mesh_ptrtype meshint( new mesh_type );
        boost::tie( it,en ) = mesh->internalElements( 0 );
        mesh->createSubmesh( *meshint, it, en );
        BOOST_CHECK_EQUAL( meshint->numElements(), std::distance( it, en ) );
        //saveGMSHMesh( _mesh=meshbdy, _filename="meshbdy" );

        using namespace Feel::vf;
        double intm3 = integrate( elements( meshint ), cst( 1. ) ).evaluate()( 0,0 );
        double intm4 = integrate( internalelements( mesh ), cst( 1. ) ).evaluate()( 0,0 );
        BOOST_CHECK_CLOSE( intm3, intm4, 1e-12 );
        //saveGMSHMesh( _mesh=meshint, _filename="meshint" );
#if 0
        auto Xh = space_type::New( mesh );
        auto u = Xh->element();
        auto v = Xh->element();
        auto D = backend->newMatrix( Xh, Xh );
        form2( Xh, Xh, D, _init=true ) = integrate( internalelements( mesh ), idt( u )*id( v ) );
        D->close();

        auto Yh = space_type::New( meshint );
        auto ui = Yh->element();
        auto vi = Yh->element();
        auto Fi = backend->newVector( Yh );
        form1( _test=Yh, _vector=Fi, _init=true ) = integrate( elements( meshint ), id( vi ) );
        Fi->close();

        auto Di = backend->newMatrix( Yh, Yh );
        form2( _test=Yh, _trial=Yh, _matrix=Di, _init=true ) = integrate( elements( meshint ), gradt( ui )*trans( grad( vi ) ) );
        Di->close();
        form2( _test=Yh, _trial=Yh, _matrix=Di ) += on( boundaryfaces( meshint ), ui, Fi, cst( 0. ) );

        backend->solve( _matrix=Di, _rhs=Fi, _solution=ui );

        boost::shared_ptr<Exporter<mesh_type> > exporter( Exporter<mesh_type>::New( "gmsh", std::string( "submesh" )
                + "_"
                + mesh_type::shape_type::name() ) );

        exporter->step( 0 )->setMesh( meshint );
        exporter->step( 0 )->add( "u", ui );
        exporter->save();

        u = vf::project( Xh, elements( mesh ), cst( 1. ) );
        ui = vf::project( Yh, elements( meshint ), cst( 1. ) );

        BOOST_CHECK_CLOSE( Di->energy( ui, ui ), intm3, 5e-12 );
        BOOST_CHECK_CLOSE( D->energy( u, u ), Di->energy( ui, ui ), 5e-12 );
#endif
    }
    boost::shared_ptr<backend_type> backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh;

};

} // Feel
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description integrationoptions( "Test Submesh options" );
    integrationoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "h value" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (hypercube, simplex, ellipsoid)" )
    ;
    return integrationoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_submesh" ,
                           "test_submesh" ,
                           "0.2",
                           "1D/2D/3D submesh",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2010 Université Joseph Fourier (Grenoble I)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( submesh )

typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<1> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<3> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_submesh, T, dim_types )
{
    BOOST_TEST_MESSAGE( "Test submesh (" << T::value << "D)" );
    Feel::test_submesh<double,T::value> t;
    t();
    BOOST_TEST_MESSAGE( "Test submesh (" << T::value << "D) done." );
}

BOOST_AUTO_TEST_SUITE_END()

#if 0
int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    Feel::Assert::setLog( "test_submesh.assert" );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}
#endif

