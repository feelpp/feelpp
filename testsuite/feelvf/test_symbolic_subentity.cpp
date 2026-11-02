// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-3.0-or-later

#define BOOST_TEST_MODULE symbolic_subentity
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

template <int dim>
void checkSubentitySymbols()
{
    auto mesh = createGMSHMesh(
        _mesh=new Mesh<Simplex<dim>>,
        _desc=domain( _name="symbolic-subentity-"+std::to_string(dim),
                      _shape="hypercube", _dim=dim, _h=0.3 ),
        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_ADD_ELEMENTS_INFO );
    auto Xh = Pch<1>( mesh );
    auto Vh = Pchv<1>( mesh );
    auto u = Xh->element();
    auto v = Vh->element();
    u.on( _range=elements( mesh ), _expr=Px()+2*Py() );
    v.on( _range=elements( mesh ), _expr=one()*( Px()+Py() ) );
    // Toolbox symbol tables contain normal derivatives even when the selected
    // expression only needs values. Points and 3D edges have no unique normal.
    auto symbols = symbolsExpr(
        symbolExpr( "u", idv( u ) ),
        symbolExpr( "du", dnv( u ) ),
        symbolExpr( "dv", dnv( v ), SymbolExprComponentSuffix( dim, 1 ) ) );
    auto expression = expr( expr( "u+1:u" ), symbols );
    auto result = Xh->element();
    auto expected = Xh->element();
    expected.on( _range=elements( mesh ), _expr=Px()+2*Py()+1 );

    auto pointRange = points( mesh, entity_process_t::ALL );
    BOOST_REQUIRE_GT( nelements( pointRange, true ), 0 );
    result.on( _range=pointRange, _expr=expression );
    sync( result );
    BOOST_CHECK_SMALL( normL2( _range=elements( mesh ), _expr=idv(result)-idv(expected) ), 1e-11 );
    if constexpr ( dim == 3 )
    {
        result.zero();
        auto edgeRange = edges( mesh, entity_process_t::ALL );
        BOOST_REQUIRE_GT( nelements( edgeRange, true ), 0 );
        result.on( _range=edgeRange, _expr=expression );
        sync( result );
        BOOST_CHECK_SMALL( normL2( _range=elements( mesh ), _expr=idv(result)-idv(expected) ), 1e-11 );
    }

    // Normal symbols and normal derivatives must still work on faces.
    auto normalExpression = expr( expr( "nx+2*ny:nx:ny" ), symbols );
    BOOST_CHECK_SMALL( normL2( _range=boundaryfaces( mesh ),
                              _expr=normalExpression-( Nx()+2*Ny() ) ), 1e-11 );
    auto derivativeExpression = expr( expr( "du:du" ), symbols );
    BOOST_CHECK_SMALL( normL2( _range=boundaryfaces( mesh ),
                              _expr=derivativeExpression-( Nx()+2*Ny() ) ), 1e-11 );
}

BOOST_AUTO_TEST_CASE( points_2d ) { checkSubentitySymbols<2>(); }
BOOST_AUTO_TEST_CASE( points_and_edges_3d ) { checkSubentitySymbols<3>(); }
