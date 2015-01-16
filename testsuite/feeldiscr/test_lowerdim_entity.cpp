/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-06-25

  Copyright (C) 2013 Université de Strasbourg

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
   \file test_lowerdim_entity.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-06-25
 */
#include <sstream>
#include <boost/timer.hpp>

#define BOOST_TEST_MODULE lower dimensional entities testsuite
#include <testsuite/testsuite.hpp>

#include <boost/mpl/list.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( lowerdimsuite )
//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
typedef boost::mpl::list<boost::mpl::int_<2> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( evaluate1, T, dim_types )
{
    BOOST_TEST_MESSAGE( "check 1 vf::sum and vf::evaluate for dim = " << T::value << "\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef Mesh<Simplex<1,1,T::value> > sub_mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    const std::string shape="hypercube";
    mesh_ptrtype mesh2d = createGMSHMesh( _mesh=new mesh_type,
                                          _desc=domain( _name=( boost::format( "simplex-%1%" )  % T::value ).str() ,
                                                        _usenames=true,
                                                        _addmidpoint=false,
                                                        _shape=shape,
                                                        _h=doption(_name="gmsh.hsize") ),
                                          _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    auto mesh = createSubmesh( mesh2d, boundaryfaces( mesh2d ), 0 );

    double exact_perim = 0.;
    if ( shape == "hypercube" )
        exact_perim = 4.;
    else
        exact_perim = 2.+math::sqrt(2.);

    typedef FunctionSpace<sub_mesh_type,bases<Lagrange<1, Scalar> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    space_ptrtype P1h = space_type::New( mesh );
    double perim = integrate( _range=elements(mesh), _expr=cst(1.) ).evaluate()( 0, 0 );
    BOOST_CHECK_CLOSE( perim, exact_perim, 1e-10 );
    BOOST_TEST_MESSAGE( "perimeter of square " << perim );
    auto u = P1h->element();
    u.setOnes();
    auto m = form2(_trial=P1h, _test=P1h );
    m = integrate( elements(mesh), idt(u)*id(u) );
    auto perim2 = m( u, u );
    BOOST_CHECK_CLOSE( perim2, perim, 1e-10 );

    auto s = form2(_trial=P1h, _test=P1h );
    s = integrate( elements(mesh), gradt(u)*trans(grad(u)) );
    auto zero = s( u, u );
    BOOST_CHECK_SMALL( zero, 1e-10 );

    u = project( _space=P1h, _range=elements(mesh), _expr=Px()+Py() );
    //u = project( _space=P1h, _range=elements(mesh), _expr=Px() );

    auto perim3 = s( u, u );
    BOOST_CHECK_CLOSE( perim3, perim, 1e-10 );

    double perim4 = integrate( _range=elements(mesh), _expr=gradv(u)*trans(gradv(u)) ).evaluate()( 0, 0 );
    BOOST_CHECK_CLOSE( perim4, exact_perim, 1e-10 );
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
