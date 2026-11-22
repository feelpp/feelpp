/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

// Static condensation must accept a right-hand side whose potential block
// (block 1) has no contribution on some or all elements, for example a zero
// source term with only boundary data. condense() used to call F1K.at(K)
// unconditionally and threw std::out_of_range on such elements.
//
// The test never assembles rhs(1_c). HDG of order 1 reproduces linear
// solutions exactly, so the computed solution is checked too.

#define BOOST_TEST_MODULE test_staticcondensation_rhs
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_WITH_ABOUT_NO_OPTIONS( Feel::makeAboutDefault( "test_staticcondensation_rhs" ) )

BOOST_AUTO_TEST_SUITE( staticcondensation_rhs_suite )

// Flux / potential / trace, Dirichlet on the whole boundary through the trace.
// Exact solution: p = x + y, u = -grad p = (-1,-1), div u = 0.
BOOST_AUTO_TEST_CASE( hdg_three_spaces_without_potential_rhs )
{
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto tau = cst( 1.0 );
    auto p_exact = Px() + Py();

    auto Vh = Pdhv<1>( mesh );
    auto Wh = Pdh<1>( mesh );
    auto face_mesh = createSubmesh( _mesh=mesh, _range=faces( mesh ), _update=0 );
    auto Mh = Pdh<1>( face_mesh );
    auto ps = product( Vh, Wh, Mh );

    backend( _rebuild=true );
    auto a = blockform2( ps, solve::strategy::static_condensation, backend() );
    auto rhs = blockform1( ps, solve::strategy::static_condensation, backend() );

    auto u = Vh->element(), v = Vh->element();
    auto p = Wh->element(), w = Wh->element();
    auto phat = Mh->element(), l = Mh->element();

    // No source term: rhs(1_c) is intentionally never assembled.
    rhs( 2_c ) += integrate( _range=boundaryfaces( mesh ), _expr=id( l ) * p_exact );

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ), _expr=trans( idt( u ) ) * id( v ) );
    a( 0_c, 1_c ) += integrate( _range=elements( mesh ), _expr=-idt( p ) * div( v ) );
    a( 0_c, 2_c ) += integrate( _range=internalfaces( mesh ),
                                _expr=idt( phat ) * ( leftface( normal( v ) ) + rightface( normal( v ) ) ) );
    a( 0_c, 2_c ) += integrate( _range=boundaryfaces( mesh ), _expr=idt( phat ) * normal( v ) );

    a( 1_c, 0_c ) += integrate( _range=elements( mesh ), _expr=id( w ) * divt( u ) );
    a( 1_c, 1_c ) += integrate( _range=internalfaces( mesh ),
                                _expr=tau * ( leftfacet( idt( p ) ) * leftface( id( w ) ) +
                                              rightfacet( idt( p ) ) * rightface( id( w ) ) ) );
    a( 1_c, 1_c ) += integrate( _range=boundaryfaces( mesh ), _expr=tau * id( w ) * idt( p ) );
    a( 1_c, 2_c ) += integrate( _range=internalfaces( mesh ),
                                _expr=-tau * idt( phat ) * ( leftface( id( w ) ) + rightface( id( w ) ) ) );
    a( 1_c, 2_c ) += integrate( _range=boundaryfaces( mesh ), _expr=-tau * idt( phat ) * id( w ) );

    a( 2_c, 0_c ) += integrate( _range=internalfaces( mesh ),
                                _expr=id( l ) * ( leftfacet( normalt( u ) ) + rightfacet( normalt( u ) ) ) );
    a( 2_c, 1_c ) += integrate( _range=internalfaces( mesh ),
                                _expr=tau * id( l ) * ( leftfacet( idt( p ) ) + rightfacet( idt( p ) ) ) );
    // 0.5: condensed assembly visits each internal face from both sides (as in qs_hdg_laplacian)
    a( 2_c, 2_c ) += integrate( _range=internalfaces( mesh ), _expr=-0.5 * tau * idt( phat ) * id( l ) );
    a( 2_c, 2_c ) += integrate( _range=boundaryfaces( mesh ), _expr=idt( phat ) * id( l ) );

    rhs.close();
    a.close();

    auto U = ps.element();
    BOOST_REQUIRE_NO_THROW( a.solve( _solution=U, _rhs=rhs, _condense=true ) );

    double const pError = normL2( _range=elements( mesh ), _expr=idv( U( 1_c ) ) - p_exact );
    double const uError = normL2( _range=elements( mesh ), _expr=idv( U( 0_c ) ) - vec( cst( -1. ), cst( -1. ) ) );
    BOOST_CHECK_SMALL( pError, 1e-7 );
    BOOST_CHECK_SMALL( uError, 1e-7 );
}

BOOST_AUTO_TEST_SUITE_END()
