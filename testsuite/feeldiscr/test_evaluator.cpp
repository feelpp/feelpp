/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-05-28

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file test_evaluator.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-05-28
 */
#include <sstream>
#include <boost/timer.hpp>

#define BOOST_TEST_MODULE evaluate testsuite
#include <testsuite/testsuite.hpp>

#include <boost/mpl/list.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( evaluatesuite )
//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
typedef boost::mpl::list<boost::mpl::int_<2> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( evaluate1, T, dim_types )
{



    BOOST_TEST_MESSAGE( "check 1 vf::sum and vf::evaluate for dim = " << T::value << "\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "simplex-%1%" )  % T::value ).str() ,
                                                _usenames=true,
                                                _addmidpoint=false,
                                                _shape="simplex",
                                                _h=0.025 ),
                                        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    typedef FunctionSpace<mesh_type,bases<Lagrange<1, Scalar> > > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    space_ptrtype P1h = space_type::New( mesh );
    using namespace vf;
    auto e1 = normLinf( _range=elements(mesh), _pset=_Q<5>(), _expr=sin(2*constants::pi<double>()*Px()) );
    BOOST_CHECK_CLOSE( e1.value(), 1., 1e-1 );
    BOOST_TEST_MESSAGE( "maximum absolute value atteined at " << e1.arg() );

    BOOST_TEST_MESSAGE("Checking on faces...");
    auto e2 = normLinf( _range=boundaryfaces(mesh), _pset=_Q<5>(), _expr=sin(2*constants::pi<double>()*Px()) );
    BOOST_CHECK_CLOSE( e2.value(), 1., 1e-1 );
    BOOST_TEST_MESSAGE( "maximum over boundary absolute value atteined at " << e2.template get<1>() );

    BOOST_TEST_MESSAGE("Checking on faces...");
    auto e3 = minmax( _range=boundaryfaces(mesh), _pset=_Q<5>(), _expr=sin(2*constants::pi<double>()*Px()) );
    BOOST_CHECK_CLOSE( e3.min(), -1., 1e-1 );
    BOOST_CHECK_CLOSE( e3.max(), 1., 1e-1 );
    BOOST_TEST_MESSAGE( "minimum over boundary absolute value atteined at " << e3.argmin() );
    BOOST_TEST_MESSAGE( "maximum over boundary absolute value atteined at " << e3.argmax() );

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
