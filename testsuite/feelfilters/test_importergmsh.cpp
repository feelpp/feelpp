/* -*- mode: c++: coding: utf-8 -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-16

  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_importergmsh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-06-16
 */
#include <sstream>

#include <boost/timer.hpp>
// Boost.Test
// make sure that the init_unit_test function is defined by UTF
#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE mesh filter testsuite
// disable the main function creation, use our own
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


using boost::unit_test::test_suite;


#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

template<int Dim>
void
checkCreateGmshMesh( std::string const& shape, std::string const& convex = "Simplex" )
{
    typedef Mesh<Simplex<Dim,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh;
    // simplex
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name=(boost::format( "%1%-%2%-%3%" )  % shape % convex % Dim).str() ,
                                         _convex=convex,
                                         _usenames=true,
                                         _addmidpoint=false,
                                         _shape=shape,
                                         _dim=Dim,
                                         _h=0.5 ) );

    auto neumann = markedfaces( mesh, "Neumann" );
    BOOST_CHECK_NE( std::distance( neumann.get<1>(), neumann.get<2>() ), 0 );
    auto dirichlet = markedfaces( mesh, "Dirichlet" );
    BOOST_CHECK_NE( std::distance( dirichlet.get<1>(), dirichlet.get<2>() ), 0 );
    BOOST_CHECK_EQUAL( std::distance( neumann.get<1>(), neumann.get<2>() )+
                       std::distance( dirichlet.get<1>(), dirichlet.get<2>() ),
                       std::distance( mesh->beginFaceOnBoundary(), mesh->endFaceOnBoundary() ) );
}
BOOST_AUTO_TEST_SUITE( gmshsuite )

typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( gmshsimplex, T, dim_types )
{
    checkCreateGmshMesh<T::value>( "simplex" );
}
BOOST_AUTO_TEST_CASE_TEMPLATE( gmshhypercube_simplex, T, dim_types )
{
    checkCreateGmshMesh<T::value>( "hypercube" );
}
BOOST_AUTO_TEST_CASE_TEMPLATE( gmshhypercube_hypercube, T, dim_types )
{
    checkCreateGmshMesh<T::value>( "hypercube", "Hypercube" );
}
BOOST_AUTO_TEST_CASE_TEMPLATE( gmshellipsoid, T, dim_types )
{
    checkCreateGmshMesh<T::value>( "ellipsoid" );
}

BOOST_AUTO_TEST_CASE( gmshgeo )
{
    typedef Mesh<Simplex<2,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh;
    // simplex
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=geo( _filename="feel.geo",
                                      _dim=2,
                                      _order=1,
                                      _h=0.2 ) );
    auto letters = markedfaces( mesh, "letters" );
    BOOST_CHECK_NE( std::distance( letters.get<1>(), letters.get<2>() ), 0 );
    auto wall = markedfaces( mesh, "wall" );
    BOOST_CHECK_NE( std::distance( wall.get<1>(), wall.get<2>() ), 0 );
    auto inlet = markedfaces( mesh, "inlet" );
    BOOST_CHECK_NE( std::distance( inlet.get<1>(), inlet.get<2>() ), 0 );
    auto outlet = markedfaces( mesh, "outlet" );
    BOOST_CHECK_NE( std::distance( outlet.get<1>(), outlet.get<2>() ), 0 );
    auto markedelts = markedelements( mesh, "feel" );
    BOOST_CHECK_EQUAL( std::distance( mesh->beginElement(), mesh->endElement() ),
                       std::distance( markedelts.get<1>(), markedelts.get<2>() ) );
}
BOOST_AUTO_TEST_CASE_TEMPLATE( gmshimportexport, T, dim_types )
{
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh,meshimp;
    // simplex
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name=(boost::format( "simplex-%1%" )  % T::value).str() ,
                                         _usenames=true,
                                         _addmidpoint=false,
                                         _shape="simplex",
                                         _dim=T::value,
                                         _h=0.5 ) );

    std::ostringstream estr;
    estr << "gmshexp-" << T::value;
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    export_ptrtype exporter( Exporter<mesh_type>::New( "gmsh", estr.str() ) );
    exporter->step(0)->setMesh( mesh );
    exporter->save();
    std::ostringstream fstr;
    fstr << "gmshexp-" << T::value << "-1_0.msh";
    meshimp = loadGMSHMesh( _mesh=new mesh_type,
                            _filename=fstr.str(),
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    auto neumann1 = markedfaces( mesh, "Neumann" );
    auto neumann2 = markedfaces( meshimp, mesh->markerName("Neumann") );
    BOOST_CHECK_EQUAL( std::distance( neumann1.get<1>(), neumann1.get<2>() ),
                       std::distance( neumann2.get<1>(), neumann2.get<2>() ) );
    auto dirichlet1 = markedfaces( mesh, "Dirichlet" );
    auto dirichlet2 = markedfaces( meshimp, mesh->markerName("Dirichlet") );
    BOOST_CHECK_EQUAL( std::distance( dirichlet1.get<1>(), dirichlet1.get<2>() ),
                       std::distance( dirichlet2.get<1>(), dirichlet2.get<2>() ) );
    BOOST_CHECK_EQUAL( std::distance( mesh->beginFaceOnBoundary(), mesh->endFaceOnBoundary() ),
                       std::distance( meshimp->beginFaceOnBoundary(), meshimp->endFaceOnBoundary() ) );
    BOOST_CHECK_EQUAL( std::distance( mesh->beginElement(), mesh->endElement() ),
                       std::distance( meshimp->beginElement(), meshimp->endElement() ) );
}

BOOST_AUTO_TEST_SUITE_END()

int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}
