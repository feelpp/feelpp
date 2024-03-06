/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Guillaume Dollé <gdolle@unistra.fr>
 Date: 26 Mar 2015

 Copyright (C) 2015 Feel++ Consortium

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

#define BOOST_TEST_MODULE test_on_dofs
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

template<int Dim, int PolyOrder=1>
void runTestAssign()
{
    // Check the on keyword for different topological entities.

    using mesh_type = Mesh< Simplex<Dim> >;
    auto mesh = loadMesh( _mesh=new mesh_type );
    auto Xh = Pch<PolyOrder>(mesh);
    auto p = Xh->element();
    auto l = Xh->element();
    auto s = Xh->element();
    auto v = Xh->element();

    size_type np = nelements( markedpoints( mesh,"P"),true );
    BOOST_CHECK( np > 0 );
    p.on( _range=markedpoints( mesh,"P"), _expr=cst(42.) );
    sync( p );
    double sump = p.sum();
    BOOST_CHECK_SMALL( sump - 42, 1e-12 );

    size_type ne = nelements( markededges( mesh,"L"),true );
    BOOST_CHECK( ne > 0 );
    l.on( _range=markededges( mesh,"L"), _expr=cst(42.) );
    sync( l );
    double maxl = l.max();
    BOOST_CHECK_SMALL( maxl - 42, 1e-12 );

    s.on( _range=markedfaces( mesh,"S"), _expr=cst(42.) );
    double ints = integrate(_range=markedfaces( mesh,"S"),_expr=idv(s) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( ints - 42, 1e-12 );

    v.on( _range=elements( mesh ), _expr=cst(42.) );
    double intv = integrate(_range=elements( mesh ),_expr=idv(v) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( intv - 42*27, 1e-10 );

#if 0
    auto e = exporter( _mesh=mesh );
    e->add( "p", p);
    e->add( "l", l);
    e->add( "s", s);
    e->add( "v", v);
    e->save();
#endif
}


template<int Dim, int PolyOrder=1>
void runTestElimination()
{
    using mesh_type = Mesh< Simplex<Dim> >;
    auto mesh = loadMesh( _mesh=new mesh_type );
    auto Xh = Pch<PolyOrder>( mesh );

    auto u = Xh->element();
    auto g = Xh->element();
    g.on(_range=elements(mesh),_expr=Px()+Py());

    auto mat = backend()->newMatrix(_test=Xh,_trial=Xh);
    auto rhs = backend()->newVector( Xh );
    auto l = form1( _test=Xh,_vector=rhs );
    //l = integrate(_range=elements(mesh),
    //              _expr=id(u));
    auto a = form2( _trial=Xh, _test=Xh,_matrix=mat );
    //a = integrate(_range=elements(mesh),
    //_expr=inner(gradt(u),grad(u)) );
    a+=on(_range=elements( mesh ), _rhs=rhs, _element=u, _expr=idv(g)*cst(2.) );
    a+=on(_range=markedfaces( mesh,"S"), _rhs=rhs, _element=u, _expr=idv(g)*cst(12.) );
    a+=on(_range=markededges( mesh,"L"), _rhs=rhs, _element=u, _expr=idv(g)*cst(22.) );
    a+=on(_range=markedpoints( mesh,"P"), _rhs=rhs, _element=u, _expr=idv(g)*cst(32.) );

    backend(_rebuild=true)->solve(_matrix=mat,_rhs=rhs,_solution=u );

    auto err = Xh->element();

    err.on( _range=elements( mesh ), _expr=abs(idv(u)-idv(g)*cst(2.)) );
    auto dofsOnElt = err.functionSpace()->dofs( elements( mesh), ComponentType::NO_COMPONENT, true );
    sync(err, "=", dofsOnElt); // useless here but just a check

    err.on( _range=markedfaces( mesh,"S"), _expr=abs(idv(u)-idv(g)*cst(12.)) );
    auto dofsOnFace = err.functionSpace()->dofs( markedfaces( mesh,"S"), ComponentType::NO_COMPONENT, true );
    sync(err, "=", dofsOnFace);

    err.on( _range=markededges( mesh,"L"), _expr=abs(idv(u)-idv(g)*cst(22.)) );
    auto dofsOnEdge = err.functionSpace()->dofs( markededges( mesh,"L"), ComponentType::NO_COMPONENT, true );
    sync(err, "=", dofsOnEdge);

    err.on( _range=markedpoints( mesh,"P"), _expr=abs(idv(u)-idv(g)*cst(32.)) );
    auto dofsOnPoint = err.functionSpace()->dofs( markedpoints( mesh,"P"), ComponentType::NO_COMPONENT, true );
    sync(err, "=", dofsOnPoint);

    BOOST_CHECK_SMALL( err.sum(), 1e-10 );
}


template<int Dim, int PolyOrder=1,template<class Convex, uint16_type Order, typename T> class PointSetT = PointSetFekete>
void runTestAssignMeshRelated()
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<Dim,1>>);
    auto therange = markedfaces(mesh,"Wall1");
    auto submesh = createSubmesh(_mesh=mesh,_range=therange);

    auto g = Px()+Py()+Pz();

    if constexpr ( PolyOrder>0 )
    {
        auto Vh = Pch<PolyOrder,double,PointSetT>( mesh );
        auto Qh = Pch<PolyOrder, double, PointSetT>( submesh );
        auto v = Vh->element();
        auto q1 = Qh->element();
        auto q2 = Qh->element();

        v.on(_range=therange,_expr=g);
        q1.on(_range=elements(submesh),_expr=g);
        q2.on(_range=therange,_expr=g);
        double int1 = integrate(_range=therange,_expr=idv(v)).evaluate()(0,0);
        double int2 = integrate(_range=elements(submesh),_expr=idv(q1)).evaluate()(0,0);
        double int3 = integrate(_range=therange,_expr=idv(q2)).evaluate()(0,0);
        BOOST_CHECK_CLOSE( int1, int2, 1e-12 );
        BOOST_CHECK_CLOSE( int1, int3, 1e-12 );
    }

    auto Wh = Pdhv<PolyOrder,PointSetT>( mesh );
    auto Rh = Pdhv<PolyOrder,PointSetT>( submesh );
    auto w = Wh->element();
    auto r = Rh->element();

    w.on(_range=therange,_expr=g*N());
    r.on(_range=therange,_expr=g*N());
    double intv1 = integrate(_range=therange,_expr=inner(idv(w))).evaluate()(0,0);
    double intv2 = integrate(_range=elements(submesh),_expr=inner(idv(r))).evaluate()(0,0);
    if constexpr ( PolyOrder>0 ) // currently intv1 is equal to 0 if P0dv because it seems that we can't apply w.on(...) from a range of faces
        BOOST_CHECK_CLOSE( intv1, intv2, 1e-12 );

#if 0
    auto e1 = exporter( _mesh=mesh,_name="e1");
    e1->addRegions();
    e1->add( "v", v );
    e1->add( "w", w );
    e1->save();
    auto e2 = exporter( _mesh=submesh,_name="e2");
    e2->addRegions();
    e2->add( "q1", q1 );
    e2->add( "q2", q2 );
    e2->add( "r", r );
    e2->save();
#endif
}

FEELPP_ENVIRONMENT_WITH_ABOUT_NO_OPTIONS(Feel::makeAboutDefault("test_on_dofs"))

BOOST_AUTO_TEST_SUITE( test_on_dofs )

BOOST_AUTO_TEST_CASE( assign_3d_assign )
{
    runTestAssign<3>();
}
BOOST_AUTO_TEST_CASE( assign_3d_meshrelated_order3 )
{
    runTestAssignMeshRelated<3, 3>();
}
BOOST_AUTO_TEST_CASE( assign_3d_meshrelated_order4 )
{
    // @warning this does not work with PointSetFekete
    runTestAssignMeshRelated<3,4,PointSetEquiSpaced>();
}
BOOST_AUTO_TEST_CASE( assign_3d_meshrelated_order1 )
{
    runTestAssignMeshRelated<3,1>();
}

BOOST_AUTO_TEST_CASE( elimination_3d )
{
    runTestElimination<3>();
}

BOOST_AUTO_TEST_SUITE_END()

