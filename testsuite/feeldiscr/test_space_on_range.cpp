/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 7 june 2017

 Copyright (C) 2017 Feel++ Consortium

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

// give a name to the testsuite
#define BOOST_TEST_MODULE function space on range testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/concatenate.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feeldiscr/pdhv.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( space_on_range )

BOOST_AUTO_TEST_CASE( test_2d )
{
    using namespace Feel;
    typedef Mesh<Simplex<2,1>> mesh_type;
    auto mesh = loadMesh(_mesh=new mesh_type );

    typedef FunctionSpace<mesh_type,bases<Lagrange<2,Scalar> > > space_type;

    auto r1 = elements(mesh, pow(Px()-0.4,2)+pow(Py()-0.5,2) < pow(cst(0.23),2), _selector=select_elements_from_expression::with_value, _value=1 );
    auto r2 = elements(mesh, pow(Px()-0.7,2)+pow(Py()-0.5,2) < pow(cst(0.15),2), _selector=select_elements_from_expression::with_value, _value=1 );
    auto therange = concatenate(r1,r2);

    auto VhPS = space_type::New(_mesh=mesh,_range=therange);
    BOOST_TEST_MESSAGE( "VhPS->nDof() : " << VhPS->nDof() );
    auto VhFS = space_type::New(_mesh=mesh);
    BOOST_TEST_MESSAGE( "VhFS->nDof() : " << VhFS->nDof());

    // test element of space
    auto chiShapeFS = VhFS->element();
    chiShapeFS.setConstant(5.);
    chiShapeFS.on(_range=therange,_expr=cst(3.));
    auto chiShapePS = VhPS->element();
    chiShapePS.setConstant(3.);
    double err = normL2(_range=therange,_expr=idv(chiShapeFS)-idv(chiShapePS));
    BOOST_CHECK_SMALL( err,1e-12 );
    chiShapeFS.on(_range=therange,_expr=Px()*Py());
    chiShapePS.on(_range=therange,_expr=Px()*Py());
    err = normL2(_range=therange,_expr=idv(chiShapeFS)-idv(chiShapePS));
    BOOST_CHECK_SMALL( err,1e-12 );
    double chiFSnorm2 = chiShapeFS.l2Norm();
    double chiPSnorm2 = chiShapePS.l2Norm();
    BOOST_CHECK( chiFSnorm2 > chiPSnorm2 );

    // test a laplacian solve (strong dirichlet)
    auto u = VhPS->element("u");
    auto v = VhPS->element("v");
    auto g = VhPS->element();
    g.on(_range=therange,_expr=cst(1.));
    BOOST_CHECK( VhPS->dof()->hasMeshSupport() );
    VhPS->dof()->meshSupport()->updateBoundaryInternalFaces();
    auto myboundaryfaces = VhPS->dof()->meshSupport()->rangeBoundaryFaces();
    auto l = form1( _test=VhPS );
    l = integrate(_range=therange,
                  _expr=id(v));
    auto a = form2( _trial=VhPS, _test=VhPS);
    a = integrate(_range=therange,
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=myboundaryfaces, _rhs=l, _element=u, _expr=idv(g) );
    a.solve(_rhs=l,_solution=u);

    // export results
    auto e = exporter( _mesh=mesh,_name="test2d" );
    e->addRegions();
    e->add( "u", u );
    e->add( "chi-partial-support", chiShapeFS );
    e->add( "chi-full-support", chiShapePS );
    e->save();

    // test a laplacian solve (weak dirichlet)
    auto u2 = VhPS->element();
    auto v2 = VhPS->element();
    auto l2 = form1( _test=VhPS );
    l2 = integrate(_range=therange,
                   _expr=id(v2));
    auto a2 = form2( _trial=VhPS, _test=VhPS);
    a2 = integrate(_range=therange,
                   _expr=gradt(u2)*trans(grad(v2)) );
    double penaldir = 10.;
    a2 += integrate( _range=myboundaryfaces,
                     _expr=  ( -( gradt( u2 )*N() )*id( v2 )
                               -( grad( v2 )*N() )*idt( u2 )
                              +penaldir*id( v2 )*idt( u2 )/hFace() ) );
    l2 += integrate( _range=myboundaryfaces,
                    _expr=idv(g)*( -grad( v2 )*N()
                                   + penaldir*id( v2 )/hFace() ) );
    a2.solve(_rhs=l2,_solution=u2);
    double diffsol = normL2(_range=therange,_expr=idv(u)-idv(u2));
    BOOST_CHECK_SMALL( diffsol,1e-3 );
}

BOOST_AUTO_TEST_CASE( test_composite_2d )
{
    using namespace Feel;
    typedef Mesh<Simplex<2,1>> mesh_type;
    auto mesh = loadMesh(_mesh=new mesh_type );

    typedef FunctionSpace<mesh_type,bases<Lagrange<2,Vectorial>,Lagrange<2,Scalar> > > space_type;

    auto r1 = elements(mesh, pow(Px()-0.4,2)+pow(Py()-0.5,2) < pow(cst(0.23),2), _selector=select_elements_from_expression::with_value, _value=1 );
    auto r2 = elements( mesh, pow( Px() - 0.7, 2 ) + pow( Py() - 0.5, 2 ) < pow( cst( 0.15 ), 2 ), _selector = select_elements_from_expression::with_value, _value = 1 );
    auto therange = concatenate(r1,r2);
    auto VhPS1 = space_type::New(_mesh=mesh,_range=therange);
    BOOST_TEST_MESSAGE( "VhPS1->nDof() : " << VhPS1->nDof() );

    auto supp1 = std::make_shared<MeshSupport<mesh_type>>(mesh,r1);
    auto supp2 = std::make_shared<MeshSupport<mesh_type>>(mesh,r2);
    auto VhPS2 = space_type::New(_mesh=mesh,_range=fusion::make_vector(supp1,supp2));
    BOOST_TEST_MESSAGE( "VhPS2->nDof() : " << VhPS2->nDof() );

    auto VhFS = space_type::New(_mesh=mesh);
    BOOST_TEST_MESSAGE( "VhFS->nDof() : " << VhFS->nDof() );

    auto fieldFS_U = VhFS->element();
    auto fieldFS_u = fieldFS_U.element<0>();
    auto fieldFS_p = fieldFS_U.element<1>();
    fieldFS_u.on(_range=therange,_expr=3*one());
    fieldFS_p.on(_range=therange,_expr=cst(4.));
    auto fieldPS1_U = VhPS1->element();
    auto fieldPS1_u = fieldPS1_U.element<0>();
    auto fieldPS1_p = fieldPS1_U.element<1>();
    fieldPS1_u.setConstant(3.);
    fieldPS1_p.setConstant(4.);
    double err1_u = normL2(_range=therange,_expr=idv(fieldFS_u)-idv(fieldPS1_u));
    BOOST_CHECK_SMALL( err1_u,1e-12 );
    double err1_p = normL2(_range=therange,_expr=idv(fieldFS_p)-idv(fieldPS1_p));
    BOOST_CHECK_SMALL( err1_p,1e-12 );

    fieldFS_U.zero();
    fieldFS_u.on(_range=r1,_expr=3*one());
    fieldFS_p.on(_range=r2,_expr=cst(4.));
    auto fieldPS2_U = VhPS2->element();
    auto fieldPS2_u = fieldPS2_U.element<0>();
    auto fieldPS2_p = fieldPS2_U.element<1>();
    fieldPS2_u.setConstant(3.);
    fieldPS2_p.setConstant(4.);
    double err2_u = normL2(_range=r1,_expr=idv(fieldFS_u)-idv(fieldPS2_u));
    BOOST_CHECK_SMALL( err2_u,1e-12 );
    double err2_p = normL2(_range=r2,_expr=idv(fieldFS_p)-idv(fieldPS2_p));
    BOOST_CHECK_SMALL( err2_p,1e-12 );
}
BOOST_AUTO_TEST_CASE( test_composite_meshes_list_2d )
{
    using namespace Feel;
    static const uint16_type nDim = 2;
    typedef Mesh<Simplex<nDim,1>> mesh_type;
    auto mesh = loadMesh(_mesh=new mesh_type );

    typedef Mesh<Simplex<nDim-1,1,nDim>> submesh_type;
    auto submesh = createSubmesh(_mesh=mesh,_range=boundaryfaces(mesh) );

    typedef FunctionSpace< meshes<mesh_type,submesh_type>,bases<Lagrange<2,Scalar>,Lagrange<2,Scalar> > > space_type;

    auto r1 = elements(mesh, pow(Px()-0.4,2)+pow(Py()-0.5,2) < pow(cst(0.23),2), _selector=select_elements_from_expression::with_value, _value=1 );
    auto r2 = elements( submesh, Px() * Py() < cst( 0.25 ), _selector = select_elements_from_expression::with_value, _value = 1 );

    auto supp1 = std::make_shared<MeshSupport<mesh_type>>(mesh,r1);
    auto supp2 = std::make_shared<MeshSupport<submesh_type>>(submesh,r2);
    auto VhPS = space_type::New(_mesh=fusion::make_vector(mesh,submesh),_range=fusion::make_vector(supp1,supp2));
    BOOST_TEST_MESSAGE( "VhPS->nDof() : " << VhPS->nDof() );

    auto VhFS = space_type::New(_mesh=fusion::make_vector(mesh,submesh));
    BOOST_TEST_MESSAGE( "VhFS->nDof() : " << VhFS->nDof() );

    auto fieldFS_U = VhFS->element();
    auto fieldFS_u1 = fieldFS_U.element<0>();
    auto fieldFS_u2 = fieldFS_U.element<1>();
    fieldFS_u1.on(_range=r1,_expr=cst(3.));
    fieldFS_u2.on(_range=r2,_expr=cst(4.));

    auto fieldPS_U = VhPS->element();
    auto fieldPS_u1 = fieldPS_U.element<0>();
    auto fieldPS_u2 = fieldPS_U.element<1>();
    fieldPS_u1.setConstant(3.);
    fieldPS_u2.setConstant(4.);
    double err_u1 = normL2(_range=r1,_expr=idv(fieldFS_u1)-idv(fieldPS_u1));
    BOOST_CHECK_SMALL( err_u1,1e-12 );
    double err_u2 = normL2(_range=r2,_expr=idv(fieldFS_u2)-idv(fieldPS_u2));
    BOOST_CHECK_SMALL( err_u2,1e-12 );

    auto r = elements(mesh);
    auto rSize = r.size();
    auto r1Size = r1.size();
    auto r2Size = r2.size();
    BOOST_TEST_MESSAGE( "r size : " << rSize );
    BOOST_TEST_MESSAGE( "r1 size : " << r1Size );
    BOOST_TEST_MESSAGE( "r2 size : " << r2Size );
    auto VhPS1 = VhPS->template functionSpace<0>();
    auto fsSize = VhFS->rangeElements<0>().size(); // no meshsupport
    auto ps1Size = VhPS1->rangeElements<0>().size(); // no composite
    auto ps2Size = VhPS->rangeElements<1>().size(); // composite
    BOOST_CHECK_EQUAL( rSize,fsSize);
    BOOST_CHECK_EQUAL( r1Size,ps1Size);
    BOOST_CHECK_EQUAL( r2Size,ps2Size);
}
BOOST_AUTO_TEST_CASE( test_extended_2d )
{
    using namespace Feel;
    typedef Mesh<Simplex<2>> mesh_type;
    auto mesh = loadMesh( _mesh=new mesh_type );
    auto r1 = elements(mesh, pow(Px()-0.4,2)+pow(Py()-0.5,2)/*+pow(Pz()-0.5,2)*/ < pow(cst(0.23),2), _selector=select_elements_from_expression::with_value, _value=1 );
    auto r2 = elements(mesh, pow(Px()-0.7,2)+pow(Py()-0.5,2)/*+pow(Pz()-0.6,2)*/ < pow(cst(0.15),2), _selector=select_elements_from_expression::with_value, _value=1 );
    auto therange = concatenate(r1,r2);
    typedef FunctionSpace<mesh_type,bases<Lagrange<1,Scalar,Discontinuous> > > space_type;
    auto Vh = space_type::New(_mesh=mesh,_extended_doftable=true,_range=therange);

    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );
    l = integrate(_range=therange,
                  _expr=id(v));

    auto a = form2( _trial=Vh, _test=Vh,
                    _pattern=size_type(Pattern::EXTENDED) );
    a = integrate(_range=therange,
                  _expr=gradt(u)*trans(grad(v)) );
    Vh->dof()->meshSupport()->updateBoundaryInternalFaces();
    auto myinternalfaces = Vh->dof()->meshSupport()->rangeInternalFaces();
    a +=integrate( _range=myinternalfaces,
                   _expr=-averaget( gradt( u ) )*jump( id( v ) )
                   -average( grad( v ) )*jumpt( idt( u ) )
                   + 50* ( trans( jumpt( idt( u ) ) )*jump( id( v ) ) )/hFace() );
    auto myboundaryfaces = Vh->dof()->meshSupport()->rangeBoundaryFaces();
    a+=on(_range=myboundaryfaces, _rhs=l, _element=u, _expr=cst(0.) );
    a.solve(_rhs=l,_solution=u,_rebuild=true);
}
BOOST_AUTO_TEST_CASE( test_integrate_boundaryfaces )
{
    using namespace Feel;
    typedef Mesh<Simplex<2>> mesh_type;
    auto mesh = loadMesh( _mesh=new mesh_type );
    auto r1 = elements(mesh, pow(Px()-0.4,2)+pow(Py()-0.5,2) < pow(cst(0.23),2), _selector=select_elements_from_expression::with_value, _value=1 );
    auto r2 = elements(mesh, pow(Px()-0.7,2)+pow(Py()-0.5,2) < pow(cst(0.15),2), _selector=select_elements_from_expression::with_value, _value=1 );
    auto therange = concatenate(r1,r2);
    typedef FunctionSpace<mesh_type,bases<Lagrange<2,Scalar> > > space_type;
    auto Vh = space_type::New(_mesh=mesh,_range=therange);

    Vh->dof()->meshSupport()->updateBoundaryInternalFaces();
    auto myboundaryfaces = Vh->dof()->meshSupport()->rangeBoundaryFaces();

    double int1 = integrate(_range=myboundaryfaces,_expr=cst(1.)).evaluate()(0,0);
    auto uPS =  Vh->element(cst(1.));

    auto l = form1( _test=Vh );
    l = integrate(_range= myboundaryfaces,_expr=id(uPS));
    l.vectorPtr()->close();
    double int1PS = inner_product( *l.vectorPtr(),uPS);
    BOOST_CHECK_SMALL( std::abs(int1PS-int1),1e-12 );

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=myboundaryfaces,_expr=id(uPS)*idt(uPS));
    a.matrixPtr()->close();
    double int2PS = a.matrixPtr()->energy(uPS,uPS);
    BOOST_CHECK_SMALL( std::abs(int2PS-int1),1e-12 );

    auto submesh = createSubmesh(_mesh=mesh,_range=myboundaryfaces);
    double int1b = integrate(_range=elements(submesh),_expr=cst(1.)).evaluate()(0,0);
    BOOST_CHECK_SMALL( std::abs(int1PS-int1b),1e-12 );
    BOOST_CHECK_SMALL( std::abs(int2PS-int1b),1e-12 );
}

BOOST_AUTO_TEST_CASE( test_integrate_different_related_mesh )
{
    using namespace Feel;
    auto mesh = unitCube();
    auto Vh = Pdhv<1>(mesh,elements(mesh, Px() < cst(0.5), _selector=select_elements_from_expression::with_value, _value=1 ), true );
    //auto Vh = Pdhv<1>(mesh);
    auto submesh = createSubmesh(_mesh=mesh, _range=faces(support(Vh)),_update=0);
    //auto submesh = createSubmesh(_mesh=mesh, _range=boundaryfaces(support(Vh)));
    auto Xh = Pdh<1>(submesh,true);

    auto u = Vh->element();
    auto v = Xh->element();
    v.setConstant(2.);
    u.setConstant(3.);

    auto rangeIntegrateBoundary = intersect(boundaryfaces(support(Vh)),boundaryfaces(mesh)); // there is a problem
    auto a = form2(_test=Xh, _trial=Vh);
    //a = integrate(_range=boundaryfaces(mesh), _expr=idt(v)*normal(u));
    //a = integrate(_range=boundaryfaces(mesh), _expr=id(v)*inner(idt(u),N()));
    a = integrate(_range=rangeIntegrateBoundary, _expr=id(v)*inner(idt(u),one()));
    a.matrixPtr()->close();
    double theMatEnergyA = a(v,u);
    BOOST_TEST_MESSAGE( "theMatEnergyA=" << theMatEnergyA );
    auto submeshBoundarySupport = createSubmesh(_mesh=mesh, _range=rangeIntegrateBoundary/*,_update=0*/);
    double evalIntegrateA = integrate(_range=elements(submeshBoundarySupport),_expr=cst(2.)*inner(3*one(),one())).evaluate()(0,0);
    BOOST_TEST_MESSAGE( "evalIntegrateA=" << evalIntegrateA );

    BOOST_CHECK_CLOSE( theMatEnergyA, evalIntegrateA, 1e-9 );


    auto b = form2(_test=Xh, _trial=Vh);
    b = integrate(_range=internalfaces(support(Vh)/*mesh*/), _expr=id(v)*inner(idt(u),one()));
    b.matrixPtr()->close();
    double theMatEnergyB = b(v,u);
    BOOST_TEST_MESSAGE( "theMatEnergyB=" << theMatEnergyB );

    auto submeshInternalSupport = createSubmesh(_mesh=mesh, _range=internalfaces(support(Vh)),_update=0);
    double evalIntegrateB = integrate(_range=elements(submeshInternalSupport),_expr=cst(2.)*inner(3*one(),one())).evaluate()(0,0);
    BOOST_TEST_MESSAGE( "evalIntegrateB=" << evalIntegrateB );

    BOOST_CHECK_CLOSE( theMatEnergyB, 2*evalIntegrateB, 1e-9 );



    // TODO :  need to work on boundaryfaces(support(Vh)) : the support information is lost with range
    auto c = form1(_test=Xh);
    //c = integrate( _range=faces(support(Vh)), _expr=id(v) );
    //c = integrate( _range=boundaryfaces(support(Vh)), _expr=id(v) );
    c = integrate( _range=internalfaces(support(Vh)), _expr=id(v) );
    c.vectorPtr()->close();
    double evalIntLF = inner_product( *c.vectorPtr(),v);
    //double evalIntLF_check = integrate(_range=faces(support(Vh)),_expr=leftfacev(cst(2.))).evaluate()(0,0); // integrate only on one side
    // double evalIntLF_check = integrate(_range=boundaryfaces(support(Vh)),_expr=cst(2.)).evaluate()(0,0)
    //     +  2*integrate(_range=internalfaces(support(Vh)),_expr=cst(2.)).evaluate()(0,0);

    //auto submeshBoundarySupport2 = createSubmesh( _mesh=mesh, _range=internalfaces(support(Vh)), _update=0 );
    //double evalIntLF_check = integrate(_range=elements(submeshBoundarySupport2)/*boundaryfaces(support(Vh))*/,_expr=cst(2.)).evaluate()(0,0);
    //+  2*integrate(_range=internalfaces(support(Vh)),_expr=cst(2.)).evaluate()(0,0);a
    double evalIntLF_check = 2*integrate(_range=internalfaces(support(Vh)),_expr=cst(2.)).evaluate()(0,0);
    //double evalIntLF_check = 2*integrate(_range=elements(submeshBoundarySupport2),_expr=cst(2.)).evaluate()(0,0);
    BOOST_TEST_MESSAGE( "evalIntLF= " << evalIntLF );
    BOOST_CHECK_CLOSE( evalIntLF, evalIntLF_check, 1e-9 );

}

BOOST_AUTO_TEST_SUITE_END()

