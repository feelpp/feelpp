/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-10-27

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
   \file test_mortar.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-10-27
 */
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE Mortar testsuite


#include <testsuite/testsuite.hpp>

#include <boost/timer.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <feel/feelcore/environment.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/moch.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/thch.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( mortar )

BOOST_AUTO_TEST_CASE( test_no_mortar_1 )
{
    BOOST_TEST_MESSAGE( "test_no_mortar_1" );
    using namespace Feel;
    auto mesh =loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>(mesh);
    BOOST_CHECK_MESSAGE( Xh->is_mortar == false, "Space should not be mortar" ) ;
}
BOOST_AUTO_TEST_CASE( test_no_mortar_2 )
{
    BOOST_TEST_MESSAGE( "test_no_mortar_2" );
    using namespace Feel;
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = THch<1>(mesh);
    BOOST_CHECK_MESSAGE( Xh->functionSpace<0>()->is_mortar == false, "Space 0 should not be mortar" ) ;
    BOOST_CHECK_MESSAGE( Xh->functionSpace<1>()->is_mortar == false, "Space 1 should not be mortar" ) ;

    BOOST_TEST_MESSAGE( "test_no_mortar_2 done" );
}

typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<4>  > order_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_mortar_1, T, order_types )
{
    using namespace Feel;
    Feel::Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/test_mortar_1/h_%2%/P%3%/" )
                                         % Feel::Environment::about().appName()
                                         % option(_name="gmsh.hsize").template as<double>()
                                         % T::value );
    BOOST_TEST_MESSAGE( "test_mortar_1 for order : " << T::value );
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<1,1,2>> );
    auto Xh = Moch<T::value>(mesh);
    BOOST_CHECK_MESSAGE(Xh->is_mortar == true, "Space should be mortar" ) ;
    BOOST_CHECK_MESSAGE(Xh->isMortar() == true, "Space should be mortar" ) ;


    BOOST_TEST_MESSAGE( "n elements : " << nelements( elements(mesh) ) );
    BOOST_TEST_MESSAGE( "n boundary elements : " << nelements( boundaryelements(mesh) ) );
    BOOST_TEST_MESSAGE( "n internal elements : " << nelements( internalelements(mesh) ) );

    for( auto const& dof : Xh->dof()->localDof() )
    {
        LOG(INFO) << "local dof element " << dof.first.elementId() << " id:" << dof.first.localDof()
                  << " global dof : " << dof.second.index();
    }

#if 1
    BOOST_CHECK_MESSAGE( nelements( internalelements(mesh) ) == nelements(elements(mesh))-2,
                         "invalid number of internal elements : " << nelements( internalelements(mesh) )
                         << " number of elements : " << nelements( elements(mesh) )  );
    for( auto const& e : internalelements(mesh) )
    {
        BOOST_CHECK_MESSAGE( e.isOnBoundary() == false, "element " << e.id() << " should be internal and not on boundary" );
        auto const& dofit = Xh->dof()->localDof( e.id() );
        int nldof = std::distance( dofit.first, dofit.second );
        BOOST_CHECK_MESSAGE(  nldof == T::value+1, "Invalid number of dof :  " << nldof );
    }
#endif
    BOOST_CHECK_MESSAGE( nelements( boundaryelements(mesh) ) == 2,
                         "invalid number of boundary elements : " << nelements( boundaryelements(mesh) ) );
    for( auto const& e : boundaryelements(mesh) )
    {
        BOOST_CHECK_MESSAGE( e.isOnBoundary() == true, "element " << e.id() << " should be on boundary" );
        auto const& dofit = Xh->dof()->localDof( e.id() );
        int nldof = std::distance( dofit.first, dofit.second );
        BOOST_CHECK_MESSAGE(  nldof == T::value, "Invalid number of dof :  " << nldof );
    }

    BOOST_TEST_MESSAGE( "test_mortar_1 for order : " << T::value << " done.");
}


BOOST_AUTO_TEST_CASE_TEMPLATE( test_mortar_integrate, T, order_types )
{
    using namespace Feel;
    Feel::Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/test_mortar_integrate/h_%2%/P%3%/" )
                                         % Feel::Environment::about().appName()
                                         % option(_name="gmsh.hsize").template as<double>()
                                         % T::value );


    BOOST_TEST_MESSAGE( "test_mortar_integrate for order : " << T::value );
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<1,1,2>> );
    //auto mesh = loadMesh( _mesh=new Mesh<Simplex<2,1,2> > );
    auto mesh2 = loadMesh( _mesh=new Mesh<Simplex<1,1,2> >, _h=option(_name="gmsh.hsize2").template as<double>() );
    // extract boundary faces marked by 1
    //auto mortarmesh = createSubmesh(mesh, markedfaces(mesh,1),Environment::worldComm() );
    auto Xh = Pch<T::value>(mesh);
    auto Vh = Pch<T::value>(mesh2);
    auto Mh = Moch<T::value>(mesh);

    BOOST_CHECK_MESSAGE(Mh->is_mortar == true, "Space should be mortar" ) ;
    BOOST_CHECK_MESSAGE(Mh->isMortar() == true, "Space should be mortar" ) ;


    BOOST_TEST_MESSAGE( "n elements : " << nelements( elements(mesh) ) );
    BOOST_TEST_MESSAGE( "n boundary elements : " << nelements( boundaryelements(mesh) ) );
    BOOST_TEST_MESSAGE( "n internal elements : " << nelements( internalelements(mesh) ) );

    for( auto const& dof : Mh->dof()->localDof() )
    {
        LOG(INFO) << "local dof element " << dof.first.elementId() << " id:" << dof.first.localDof()
                  << " global dof : " << dof.second.index();
    }

    BOOST_TEST_MESSAGE( "build Xh element" );
    auto u = Xh->element();
    auto v = Vh->element();
    BOOST_TEST_MESSAGE( "build Mh element" );
    auto l = Mh->element();
    BOOST_TEST_MESSAGE( "build bilinear form c_s(Mh,Xh)" );
    auto c_s = form2(_test=Mh, _trial=Xh), cs1= form2(_test=Mh, _trial=Xh);
    BOOST_TEST_MESSAGE( "integrate" );
    c_s = integrate(_range=internalelements(mesh),_expr=idt(u)*print(id(l),"test func internal"));
    c_s += integrate(_range=boundaryelements(mesh),_expr=idt(u)*print(id(l),"test func boundary"));
    cs1 = integrate(_range=internalelements(mesh),_expr=idt(u)*id(l));
    cs1 += integrate(_range=boundaryelements(mesh),_expr=idt(u));
    BOOST_TEST_MESSAGE( "printMatlab" );
    c_s.matrixPtr()->printMatlab( "C_s.m" );
    cs1.matrixPtr()->printMatlab( "C_s1.m" );
    BOOST_TEST_MESSAGE( "build bilinear form c_m(Mh,Vh)" );
    auto c_m = form2(_test=Mh, _trial=Vh);
    BOOST_TEST_MESSAGE( "integrate" );
    c_m = integrate(_range=elements(mesh),_expr=idt(v)*id(l));
    c_m.matrixPtr()->printMatlab( "C_m.m" );

}
BOOST_AUTO_TEST_SUITE_END()
