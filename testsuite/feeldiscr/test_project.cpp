/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-08-09

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
   \file test_project.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-08-09
 */
#include <sstream>

#include <boost/timer.hpp>
// Boost.Test
// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE project discr testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


using boost::unit_test::test_suite;


#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

BOOST_AUTO_TEST_SUITE( projectsuite )
Feel::Environment env( boost::unit_test::framework::master_test_suite().argc,
                       boost::unit_test::framework::master_test_suite().argv );

typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( project1, T, dim_types )
{
    BOOST_TEST_MESSAGE( "check 1 vf::sum and vf::project for dim = " << T::value << "\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=(boost::format( "simplex-%1%" )  % T::value).str() ,
                                                      _usenames=true,
                                                      _addmidpoint=false,
                                                      _shape="simplex",
                                                      _h=2. ),
                                        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    typedef FunctionSpace<mesh_type,bases<Lagrange<1, Scalar> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    space_ptrtype P1h = space_type::New( mesh );
    auto p1meas = vf::sum( P1h, elements(mesh),  vf::meas() );
    for(int i = 0; i < P1h->nLocalDof(); ++i )
    {
        BOOST_CHECK_CLOSE( p1meas(i), mesh->beginElement()->measure(), 1e-13);
    }
    BOOST_CHECK_CLOSE( p1meas.sum(), mesh->numPoints()*mesh->beginElement()->measure(), 1e-13);
    BOOST_CHECK_EQUAL( mesh->beginElement()->numberOfPointElementNeighbors(), 1 );
    BOOST_CHECK_CLOSE( mesh->beginElement()->measurePointElementNeighbors(), mesh->beginElement()->measure(), 1e-13 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( project2, T, dim_types )
{
    BOOST_TEST_MESSAGE( "check 2 vf::sum and vf::project for dim = " << T::value << "\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=(boost::format( "simplex-%1%" )  % T::value).str() ,
                                                      _usenames=true,
                                                      _addmidpoint=false,
                                                      _shape="simplex",
                                                      _dim=T::value,
                                                      _h=(T::value==1)?0.49:0.5 ),
                                        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    typedef FunctionSpace<mesh_type,bases<Lagrange<1, Scalar> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    space_ptrtype P1h = space_type::New( mesh );
    auto elit = mesh->beginElement();
    auto elen = mesh->endElement();
    for( ; elit != elen; ++elit )
    {
        BOOST_TEST_MESSAGE( "  check element " << elit->id() << "\n" );
        BOOST_FOREACH( auto elid, elit->pointElementNeighborIds() )
        {
            BOOST_TEST_MESSAGE( "  element " << elit->id() << " -> element " << elid << "\n" );
        }
        switch ( T::value )
        {
        case 1:
            BOOST_CHECK_EQUAL( elit->numberOfPointElementNeighbors(), 2 );
            break;
        case 2:
            BOOST_CHECK_EQUAL( elit->numberOfPointElementNeighbors(), 2 );
            break;
        case 3:
            BOOST_CHECK_EQUAL( elit->numberOfPointElementNeighbors(), 4 );
            break;
        }
        BOOST_CHECK_CLOSE( elit->measurePointElementNeighbors(), mesh->measure(), 1e-13 );
    }
}

BOOST_AUTO_TEST_SUITE_END()

#if 0
int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}

#endif
