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
#if 0
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
//typedef boost::mpl::list<boost::mpl::int_<1>  > order_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_mortar_1, T, order_types )
{
    using namespace Feel;
    Feel::Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/test_mortar_1/h_%2%/P%3%/" )
                                         % Feel::Environment::about().appName()
                                         % doption(_name="gmsh.hsize")
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
                                         % doption(_name="gmsh.hsize")
                                         % T::value );


    BOOST_TEST_MESSAGE( "test_mortar_integrate for order : " << T::value );
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<1,1,2>> );
    //auto mesh = loadMesh( _mesh=new Mesh<Simplex<2,1,2> > );
    auto mesh2 = loadMesh( _mesh=new Mesh<Simplex<1,1,2> >, _h=doption(_name="gmsh.hsize2") );
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
    c_m = integrate(_range=internalelements(mesh),_expr=idt(v)*id(l));
    c_m += integrate(_range=boundaryelements(mesh),_expr=idt(v));
    c_m.matrixPtr()->printMatlab( "C_m.m" );

}

BOOST_AUTO_TEST_CASE_TEMPLATE( test_mortar_integrate_submesh, T, order_types )
{
    using namespace Feel;
    Feel::Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/test_mortar_integrate_submesh/h_%2%/P%3%/" )
                                         % Feel::Environment::about().appName()
                                         % doption(_name="gmsh.hsize")
                                         % T::value );


    BOOST_TEST_MESSAGE( "test_mortar_integrate_submesh for order : " << T::value );
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2,1,2> > );
    auto mesh2 = loadMesh( _mesh=new Mesh<Simplex<2,1,2> >, _h=doption(_name="gmsh.hsize2") );

    auto testmesh = createSubmesh(mesh, markedfaces(mesh,(boost::any)1),Environment::worldComm() );
    auto trialmesh = createSubmesh(mesh2, markedfaces(mesh2,(boost::any)1),Environment::worldComm() );
    auto Xh = Pch<T::value>(testmesh);
    auto Vh = Pch<T::value>(trialmesh);
    auto Mh = Moch<T::value>(testmesh);

    BOOST_CHECK_MESSAGE(Mh->is_mortar == true, "Space should be mortar" ) ;
    BOOST_CHECK_MESSAGE(Mh->isMortar() == true, "Space should be mortar" ) ;


    BOOST_TEST_MESSAGE( "n elements : " << nelements( elements(testmesh) ) );
    BOOST_TEST_MESSAGE( "n boundary elements : " << nelements( boundaryelements(testmesh) ) );
    BOOST_TEST_MESSAGE( "n internal elements : " << nelements( internalelements(testmesh) ) );

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
    c_s = integrate(_range=internalelements(testmesh),_expr=idt(u)*print(id(l),"test func internal"));
    c_s += integrate(_range=boundaryelements(testmesh),_expr=idt(u)*print(id(l),"test func boundary"));
    cs1 = integrate(_range=internalelements(testmesh),_expr=idt(u)*id(l));
    cs1 += integrate(_range=boundaryelements(testmesh),_expr=idt(u));
    BOOST_TEST_MESSAGE( "printMatlab" );
    c_s.matrixPtr()->printMatlab( "C_s.m" );
    cs1.matrixPtr()->printMatlab( "C_s1.m" );
    BOOST_TEST_MESSAGE( "build bilinear form c_m(Mh,Vh)" );
    auto c_m = form2(_test=Mh, _trial=Vh);
    BOOST_TEST_MESSAGE( "integrate" );
    c_m = integrate(_range=internalelements(testmesh),_expr=idt(v)*id(l));
    c_m += integrate(_range=boundaryelements(testmesh),_expr=idt(v));
    c_m.matrixPtr()->printMatlab( "C_m.m" );

}
#endif
typedef boost::mpl::list<boost::mpl::int_<1>, boost::mpl::int_<2>, boost::mpl::int_<3>,
                         boost::mpl::int_<4>, boost::mpl::int_<5>, boost::mpl::int_<10>  > order_types;
//typedef boost::mpl::list<boost::mpl::int_<2>> order_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_mortar_integrate_submesh2, T, order_types )
{
    LOG(INFO) << "test_mortar_integrate_submesh2: P" << T::value << " test case";
    using namespace Feel;
    Feel::Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/test_mortar_integrate_submesh2/h_%2%/P%3%/" )
                                         % Feel::Environment::about().appName()
                                         % doption(_name="gmsh.hsize2")
                                         % T::value );


    BOOST_TEST_MESSAGE( "test_mortar_integrate_submesh2 for order : " << T::value );
    LOG(INFO) << "[test_mortar_integrate_submesh2] for order : " << T::value;
    //auto mesh = loadMesh( _mesh=new Mesh<Simplex<2,1,2> >, _h=doption(_name="gmsh.hsize2") );
    //auto mesh2 = loadMesh( _mesh=new Mesh<Simplex<2,1,2> >, _h=doption(_name="gmsh.hsize2") );

#if 1
    auto mesh = createGMSHMesh( _mesh=new Mesh<Hypercube<2,1,2> >,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name="mesh", _addmidpoint=false, _usenames=false, _shape="hypercube",
                                              _dim=2, _h=doption(_name="gmsh.hsize"),
                                              _convex="Hypercube",_substructuring=true,
                                              _xmin=0., _xmax=1., _ymin=0., _ymax=1.
                                              ),
                                _structured=1
        );

    auto mesh2 = createGMSHMesh( _mesh=new Mesh<Hypercube<2,1,2> >,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                 _desc=domain( _name="mesh2", _addmidpoint=false, _usenames=false, _shape="hypercube",
                                               _dim=2, _h=doption(_name="gmsh.hsize2"),
                                               _convex="Hypercube",_substructuring=true,
                                               _xmin=0., _xmax=1., _ymin=1., _ymax=2.
                                               ),
                                 _structured=1
        );
#endif


    auto testmesh = createSubmesh(mesh, markedfaces(mesh,"NORTH"),Environment::worldComm() );
    auto trialmesh = createSubmesh(mesh2, markedfaces(mesh2,"SOUTH"),Environment::worldComm() );
    auto Xh = Pch<T::value,double,PointSetGaussLobatto>(testmesh);
    auto Vh = Pch<T::value,double,PointSetGaussLobatto>(trialmesh);
    auto Mh = Moch<T::value,PointSetGaussLobatto>(testmesh);

    BOOST_CHECK_MESSAGE(Mh->is_mortar == true, "Space should be mortar" ) ;
    BOOST_CHECK_MESSAGE(Mh->isMortar() == true, "Space should be mortar" ) ;


    BOOST_TEST_MESSAGE( "n elements : " << nelements( elements(testmesh) ) );
    BOOST_TEST_MESSAGE( "n boundary elements : " << nelements( boundaryelements(testmesh) ) );
    BOOST_TEST_MESSAGE( "n internal elements : " << nelements( internalelements(testmesh) ) );

    for( auto const& dof : Mh->dof()->localDof() )
    {
        LOG(INFO) << "test local dof element " << dof.first.elementId() << " id:" << dof.first.localDof()
                  << " global dof : " << dof.second.index() << " pts: " << Mh->dof()->dofPoint( dof.second.index() ).template get<0>();
    }

    for( auto const& dof : Xh->dof()->localDof() )
    {
        LOG(INFO) << "trial slave local dof element " << dof.first.elementId() << " id:" << dof.first.localDof()
                  << " global dof : " << dof.second.index() << " pts: " << Xh->dof()->dofPoint( dof.second.index() ).template get<0>();
    }
    for( auto const& dof : Vh->dof()->localDof() )
    {
        LOG(INFO) << "trial master local dof element " << dof.first.elementId() << " id:" << dof.first.localDof()
                  << " global dof : " << dof.second.index() << " pts: " << Vh->dof()->dofPoint( dof.second.index() ).template get<0>();
    }

    auto r = Xh->dof()->markerToDof(std::string("CrossPoints"));
    //for( auto const& dofid : Xh->dof()->markerToDof("CrossPoints") )
    for( auto it = r.first, en = r.second; it != en; ++ it )
    {
        auto dofid = *it;
        LOG(INFO) << "dof " << dofid.second << " is at a crosspoint pts: " << Xh->dof()->dofPoint( dofid.second ).template get<0>();
    }

    auto r1 = Vh->dof()->markerToDof(std::string("CrossPoints"));
    //for( auto const& dofid :  )
    for( auto it = r1.first, en = r1.second; it != en; ++ it )
    {
        auto dofid = *it;
        LOG(INFO) << "dof " << dofid.second << " is at a crosspoint pts: " << Vh->dof()->dofPoint( dofid.second ).template get<0>();
    }

    LOG_IF(WARNING, Xh->dof()->nDof() != Mh->dof()->nDof()+2 )
        << "Invalid Mortar space dimention, trial slave side: "
        << Xh->dof()->nDof()
        << " test side : " << Mh->dof()->nDof()
        << "  (+2 : " << Mh->dof()->nDof()+2;
    BOOST_CHECK( Xh->dof()->nDof() == Mh->dof()->nDof()+2 );


    BOOST_TEST_MESSAGE( "build Xh element" );
    auto u = Xh->element(),u1=Xh->element(),u2=Xh->element(),u3=Xh->element();
    u = vf::project(_space=Xh,_range=elements(testmesh),_expr=cst(1.0));
    u1 = vf::project(_space=Xh,_range=elements(testmesh),_expr=Px());
    u2 = vf::project(_space=Xh,_range=elements(testmesh),_expr=Px()*Px());
    u3 = vf::project(_space=Xh,_range=elements(testmesh),_expr=sin(pi*Px()) );

    auto v = Vh->element(), w=Vh->element(),z=Vh->element(),zz=Vh->element();
    v = vf::project(_space=Vh,_range=elements(trialmesh),_expr=cst(1.0));
    w = vf::project(_space=Vh,_range=elements(trialmesh),_expr=Px());
    z = vf::project(_space=Vh,_range=elements(trialmesh),_expr=Px()*Px());
    zz = vf::project(_space=Vh,_range=elements(trialmesh),_expr=sin(pi*Px()));


    BOOST_TEST_MESSAGE( "build Mh element" );
    auto l = Mh->element();
    //l = vf::project(_space=Mh,_range=elements(testmesh),_expr=cst(1.0));
    l.setOnes();
#if 0
    for ( testfunc : { "1:x:y", "x:x:y", "x+y:x:y" } )
    {
        BOOST_TEST_MESSAGE( "build bilinear form c_s(Mh,Xh)" );
        auto c_s = form2(_test=Mh, _trial=Xh);
        BOOST_TEST_MESSAGE( "integrate" );
        c_s = integrate(_range=internalelements(testmesh),_expr=idt(u)*id(l));
        c_s += integrate(_range=boundaryelements(testmesh),_expr=idt(u)*id(l));

        BOOST_TEST_MESSAGE( "printMatlab" );
        c_s.matrixPtr()->printMatlab( "C_s.m" );

        BOOST_CHECK_CLOSE( c_s( l, u ), 1, 1e-12 );
        BOOST_CHECK_CLOSE( c_s( l, u1 ), 0.5, 1e-12 );
        BOOST_CHECK_CLOSE( c_s( l, u2 ), 1./3., (T::value>=2)?1e-12:10 );

        BOOST_TEST_MESSAGE( "build bilinear form c_m(Mh,Vh)" );
        auto c_m = form2(_test=Mh, _trial=Vh);
        BOOST_TEST_MESSAGE( "integrate" );
        c_m = integrate(_range=elements(testmesh),_expr=idt(v)*id(l));
        c_m.matrixPtr()->printMatlab( "C_m.m" );
        BOOST_CHECK_CLOSE( c_m( l, v ), 1, 1e-12 );
        BOOST_CHECK_CLOSE( c_m( l, w ), 0.5, 1e-12 );
        BOOST_CHECK_CLOSE( c_m( l, z ), 1./3., (T::value>=2)?1e-12:10 );


        // build matrix C_m without mortar space
        BOOST_TEST_MESSAGE( "build bilinear form c_m(Xh,Vh)" );
        auto c_m1 = form2(_test=Xh, _trial=Vh);
        BOOST_TEST_MESSAGE( "integrate" );
        c_m1 = integrate(_range=elements(testmesh),_expr=idt(v)*id(u));
        c_m1.matrixPtr()->printMatlab( "C_m1.m" );

        BOOST_TEST_MESSAGE( "build bilinear form c_s2(Xh,Xh)" );
        auto c_s2 = form2(_test=Xh, _trial=Xh);
        BOOST_TEST_MESSAGE( "integrate" );
        c_s2 = integrate(_range=elements(testmesh),_expr=idt(u)*id(u));
        c_s2.matrixPtr()->printMatlab( "C_s2.m" );
    }
#endif
    BOOST_TEST_MESSAGE( "build bilinear form c_s(Mh,Xh)" );
    auto c_s = form2(_test=Mh, _trial=Xh), cs1= form2(_test=Mh, _trial=Xh);
    BOOST_TEST_MESSAGE( "integrate" );
    c_s = integrate(_range=elements(testmesh),_expr=idt(u)*id(l),_quad1=_Q<20>(),_quad=_Q<20>());

    BOOST_TEST_MESSAGE( "printMatlab" );
    c_s.matrixPtr()->printMatlab( "C_s.m" );

    double i1 = integrate(_range=elements(testmesh), _expr=cst(1.)).evaluate()(0,0);
    BOOST_TEST_MESSAGE("integrate(1)=" << i1 );
    double i2 = integrate(_range=markedfaces(mesh,"NORTH"), _expr=cst(1.)).evaluate()(0,0);
    BOOST_TEST_MESSAGE("integrate_2(1)=" << i2 );
    BOOST_CHECK_CLOSE( c_s( l, u ), i1, 1e-12 );
    BOOST_CHECK_CLOSE( c_s( l, u1 ), 0.5, 1e-12 );
    BOOST_CHECK_CLOSE( c_s( l, u2 ), 1./3., (T::value>=2)?1e-12:10 );
    BOOST_CHECK_CLOSE( c_s( l, u3 ), 2./pi, (T::value>=2)?1e-5:1e-02 );

    BOOST_TEST_MESSAGE( "build bilinear form c_m(Mh,Vh)" );
    auto c_m = form2(_test=Mh, _trial=Vh);
    BOOST_TEST_MESSAGE( "integrate" );
    c_m = integrate(_range=elements(testmesh),_expr=idt(v)*id(l),_quad1=_Q<20>(),_quad=_Q<20>());

    c_m.matrixPtr()->printMatlab( "C_m.m" );
    BOOST_CHECK_CLOSE( c_m( l, v ), 1, 1e-12 );
    BOOST_CHECK_CLOSE( c_m( l, w ), 0.5, 1e-12 );
    BOOST_CHECK_CLOSE( c_m( l, z ), 1./3., (T::value>=2)?1e-12:10 );
    BOOST_CHECK_CLOSE( c_m( l, zz ), 2./pi, (T::value>=2)?1e-5:1e-02 );

    // build matrix C_m without mortar space
    BOOST_TEST_MESSAGE( "build bilinear form c_m(Xh,Vh)" );
    auto c_m1 = form2(_test=Xh, _trial=Vh);
    BOOST_TEST_MESSAGE( "integrate" );
    c_m1 = integrate(_range=elements(testmesh),_expr=idt(v)*id(u));
    c_m1.matrixPtr()->printMatlab( "C_m1.m" );

    BOOST_TEST_MESSAGE( "build bilinear form c_s2(Xh,Xh)" );
    auto c_s2 = form2(_test=Xh, _trial=Xh);
    BOOST_TEST_MESSAGE( "integrate" );
    c_s2 = integrate(_range=elements(testmesh),_expr=idt(u)*id(u));
    c_s2.matrixPtr()->printMatlab( "C_s2.m" );

}


BOOST_AUTO_TEST_SUITE_END()
