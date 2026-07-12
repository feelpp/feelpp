/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_hdiv_runtime_order_parity.cpp
 * @brief Runtime-order constructor parity checks for H(div)/H(curl) wrappers.
 */
#define BOOST_TEST_MODULE test_hdiv_runtime_order_parity
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/neh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelvf/vf.hpp>

#include <stdexcept>
#include <string>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template<int Dim>
using simplex_mesh_type = Mesh<Simplex<Dim,1>>;

template<int Dim>
using simplex_mesh_ptrtype = std::shared_ptr<simplex_mesh_type<Dim>>;

template<int Dim>
simplex_mesh_ptrtype<Dim>
makeOneElementSimplexMesh( std::string const& name )
{
    auto desc = std::make_shared<GmshSimplexDomain>( Dim, 1, Gmsh::GMSH_REFERENCE_DOMAIN );
    desc->setPrefix( name );
    desc->usePhysicalNames( true );
    desc->setCharacteristicLength( 2.0 );
    return createGMSHMesh( _mesh = new simplex_mesh_type<Dim>, _desc = desc, _h = desc->h() );
}

template<typename StaticSpacePtrType, typename RuntimeSpacePtrType, typename MeshPtrType, typename ExprType>
void
checkOnParity( StaticSpacePtrType const& staticSpace,
               RuntimeSpacePtrType const& runtimeSpace,
               MeshPtrType const& mesh,
               ExprType const& expr,
               double tol,
               bool compareOrder = true )
{
    BOOST_REQUIRE_EQUAL( staticSpace->nDof(), runtimeSpace->nDof() );
    if ( compareOrder )
        BOOST_CHECK_EQUAL( staticSpace->runtimeOrder(), runtimeSpace->runtimeOrder() );

    auto uStatic = staticSpace->element();
    auto uRuntime = runtimeSpace->element();
    uStatic.on( _range = elements( mesh ), _expr = expr );
    uRuntime.on( _range = elements( mesh ), _expr = expr );

    auto const parityError = normL2( _range = elements( mesh ),
                                     _expr = idv( uStatic ) - idv( uRuntime ),
                                     _quad = _Q<12>() );
    BOOST_CHECK_SMALL( parityError, tol );
}

template<typename SpacePtrType, typename MeshPtrType, typename ExprType>
void
checkInterpolationExactness( SpacePtrType const& space,
                             MeshPtrType const& mesh,
                             ExprType const& expr,
                             double tol )
{
    auto u = space->element();
    u.on( _range = elements( mesh ), _expr = expr );

    auto const error = normL2( _range = elements( mesh ),
                               _expr = idv( u ) - expr,
                               _quad = _Q<14>() );
    BOOST_CHECK_SMALL( error, tol );
}

template<int Dim, int Order>
auto
rtExactExpr()
{
    if constexpr ( Dim == 2 && Order == 0 )
        return vec( cst( 1.25 ), cst( -0.75 ) );
    else if constexpr ( Dim == 3 && Order == 0 )
        return vec( cst( 1.25 ), cst( -0.75 ), cst( 0.5 ) );
    else if constexpr ( Dim == 2 && Order == 1 )
        return vec( cst( 1.0 ) + 2.0*Px() - 0.25*Py(),
                    cst( -0.5 ) + Px() + 3.0*Py() );
    else if constexpr ( Dim == 3 && Order == 1 )
        return vec( cst( 1.0 ) + 2.0*Px() - 0.25*Py(),
                    cst( -0.5 ) + Px() + 3.0*Py(),
                    cst( 0.25 ) - Py() + 0.5*Pz() );
    else if constexpr ( Dim == 2 && Order == 2 )
        return vec( cst( 1.0 ) + Px() + Py() + Px()*Py(),
                    cst( -0.5 ) + Px()*Px() + Py()*Py() );
    else if constexpr ( Dim == 3 && Order == 2 )
        return vec( cst( 1.0 ) + Px() + Py() + Px()*Py(),
                    cst( -0.5 ) + Px()*Px() + Py()*Py(),
                    cst( 0.25 ) + Pz() + Px()*Pz() );
    else if constexpr ( Dim == 2 )
        return vec( cst( 1.0 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py(),
                    cst( -0.5 ) + Px()*Px() + Py()*Py() + Px()*Py()*Py() );
    else
        return vec( cst( 1.0 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py(),
                    cst( -0.5 ) + Px()*Px() + Py()*Py() + Px()*Py()*Py(),
                    cst( 0.25 ) + Pz() + Px()*Pz() + Px()*Py()*Pz() );
}

template<int Dim, int Order>
auto
bdmExactExpr()
{
    if constexpr ( Dim == 2 && Order == 0 )
        return vec( cst( 1.0 ) + 0.5*Px() - Py(),
                    cst( -0.25 ) + 2.0*Py() + 0.125*Px() );
    else if constexpr ( Dim == 3 && Order == 0 )
        return vec( cst( 1.0 ) + 0.5*Px() - Py(),
                    cst( -0.25 ) + 2.0*Py() + 0.125*Px(),
                    cst( 0.5 ) - Px() + 1.5*Pz() );
    else if constexpr ( Dim == 2 && Order == 1 )
        return vec( cst( 0.75 ) + Px() + Px()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px() );
    else if constexpr ( Dim == 3 && Order == 1 )
        return vec( cst( 0.75 ) + Px() + Px()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px(),
                    cst( 0.5 ) + Pz() + Px()*Pz() );
    else if constexpr ( Dim == 2 && Order == 2 )
        return vec( cst( 0.75 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px() + Py()*Py() + Px()*Py()*Py() );
    else if constexpr ( Dim == 3 && Order == 2 )
        return vec( cst( 0.75 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px() + Py()*Py() + Px()*Py()*Py(),
                    cst( 0.5 ) + Pz() + Px()*Pz() + Py()*Pz() + Px()*Py()*Pz() );
    else if constexpr ( Dim == 2 )
        return vec( cst( 0.75 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py() + Px()*Px()*Py()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px() + Py()*Py() + Px()*Py()*Py() + Px()*Px()*Px()*Py() );
    else
        return vec( cst( 0.75 ) + Px() + Py() + Px()*Py() + Px()*Px()*Py() + Px()*Px()*Py()*Py(),
                    cst( -1.25 ) + Py() + Px()*Px() + Py()*Py() + Px()*Py()*Py() + Px()*Px()*Px()*Py(),
                    cst( 0.5 ) + Pz() + Px()*Pz() + Py()*Pz() + Px()*Py()*Pz() + Px()*Px()*Py()*Pz() );
}

template<int Dim, int Order>
void
checkRtStaticDynamicParityAndExactness( simplex_mesh_ptrtype<Dim> const& mesh )
{
    auto expr = rtExactExpr<Dim, Order>();
    auto rtStatic = RTh<Order>( mesh );
    auto rtDynamic = RTh<Dynamic>( mesh, RuntimeOrder{ Order } );

    BOOST_TEST_CONTEXT( "RT dim=" << Dim << " order=" << Order )
    {
        BOOST_CHECK_EQUAL( rtDynamic->basis()->runtimeOrder(), Order );
        BOOST_CHECK_EQUAL( rtDynamic->basis()->localDofPerComponent(), rtStatic->basis()->localDofPerComponent() );
        checkOnParity( rtStatic, rtDynamic, mesh, expr, 1e-11, false );
        checkInterpolationExactness( rtStatic, mesh, expr, 1e-10 );
        checkInterpolationExactness( rtDynamic, mesh, expr, 1e-10 );
    }
}

template<int Dim, int Order>
void
checkBdmStaticDynamicParityAndExactness( simplex_mesh_ptrtype<Dim> const& mesh )
{
    auto expr = bdmExactExpr<Dim, Order>();
    auto bdmStatic = BDMh<Order>( mesh );
    auto bdmDynamic = BDMh<Dynamic>( mesh, RuntimeOrder{ Order } );

    BOOST_TEST_CONTEXT( "BDM dim=" << Dim << " order=" << Order )
    {
        BOOST_CHECK_EQUAL( bdmDynamic->basis()->runtimeOrder(), Order );
        BOOST_CHECK_EQUAL( bdmDynamic->basis()->localDofPerComponent(), bdmStatic->basis()->localDofPerComponent() );
        checkOnParity( bdmStatic, bdmDynamic, mesh, expr, 1e-11, false );
        checkInterpolationExactness( bdmStatic, mesh, expr, 1e-10 );
        checkInterpolationExactness( bdmDynamic, mesh, expr, 1e-10 );
    }
}

} // namespace

BOOST_AUTO_TEST_SUITE( hdiv_runtime_order_parity_suite )

BOOST_AUTO_TEST_CASE( rt_runtime_order_constructor_parity )
{
    auto mesh = makeOneElementSimplexMesh<2>( "rt-runtime-parity" );
    auto expr = vec( cst( 1.0 ) + Px(), cst( -0.5 ) + 2.0 * Py() );

    auto rt0Static = RTh<0>( mesh );
    auto rt0Runtime = RTh<0>( mesh, RuntimeOrder{ 0 } );
    checkOnParity( rt0Static, rt0Runtime, mesh, expr, 1e-12 );
}

BOOST_AUTO_TEST_CASE( rt_dynamic_low_order_space_parity )
{
    auto mesh2d = makeOneElementSimplexMesh<2>( "rt-dynamic-low-order-parity-2d" );
    auto mesh3d = makeOneElementSimplexMesh<3>( "rt-dynamic-low-order-parity-3d" );
    checkRtStaticDynamicParityAndExactness<2, 0>( mesh2d );
    checkRtStaticDynamicParityAndExactness<2, 1>( mesh2d );
    checkRtStaticDynamicParityAndExactness<2, 2>( mesh2d );
    checkRtStaticDynamicParityAndExactness<2, 3>( mesh2d );
    checkRtStaticDynamicParityAndExactness<3, 0>( mesh3d );
    checkRtStaticDynamicParityAndExactness<3, 1>( mesh3d );
    checkRtStaticDynamicParityAndExactness<3, 2>( mesh3d );
    checkRtStaticDynamicParityAndExactness<3, 3>( mesh3d );
}

BOOST_AUTO_TEST_CASE( rt_static_higher_order_simplex_space_layout )
{
    auto mesh = makeOneElementSimplexMesh<2>( "rt-static-simplex-higher-order-layout" );

    auto rt1Static = RTh<1>( mesh );
    auto rt2Static = RTh<2>( mesh );

    BOOST_CHECK_EQUAL( rt1Static->basis()->localDofPerComponent(), 8 );
    BOOST_CHECK_EQUAL( rt2Static->basis()->localDofPerComponent(), 15 );
    BOOST_CHECK_EQUAL( rt1Static->nLocalDof(), 8 );
    BOOST_CHECK_EQUAL( rt2Static->nLocalDof(), 15 );
    BOOST_CHECK_EQUAL( rt1Static->basis()->nDofPerEdge, 2 );
    BOOST_CHECK_EQUAL( rt2Static->basis()->nDofPerEdge, 3 );
    BOOST_CHECK_EQUAL( rt1Static->basis()->nDofPerFace, 2 );
    BOOST_CHECK_EQUAL( rt2Static->basis()->nDofPerFace, 6 );

    auto rt1Dynamic = RTh<Dynamic>( mesh, RuntimeOrder{ 1 } );
    auto rt2Dynamic = RTh<Dynamic>( mesh, RuntimeOrder{ 2 } );

    BOOST_CHECK_EQUAL( rt1Dynamic->basis()->localDofPerComponent(), 8 );
    BOOST_CHECK_EQUAL( rt2Dynamic->basis()->localDofPerComponent(), 15 );
    BOOST_CHECK_EQUAL( rt1Dynamic->nLocalDof(), 8 );
    BOOST_CHECK_EQUAL( rt2Dynamic->nLocalDof(), 15 );
    BOOST_CHECK_EQUAL( rt1Dynamic->basis()->dofPerEdge(), 2 );
    BOOST_CHECK_EQUAL( rt2Dynamic->basis()->dofPerEdge(), 3 );
    BOOST_CHECK_EQUAL( rt1Dynamic->basis()->dofPerFace(), 2 );
    BOOST_CHECK_EQUAL( rt2Dynamic->basis()->dofPerFace(), 6 );
}

BOOST_AUTO_TEST_CASE( bdm_runtime_order_constructor_parity )
{
    auto mesh = makeOneElementSimplexMesh<2>( "bdm-runtime-parity" );
    auto mesh3d = makeOneElementSimplexMesh<3>( "bdm-runtime-parity-3d" );
    auto expr = bdmExactExpr<2, 0>();

    auto bdm0Static = BDMh<0>( mesh );
    auto bdm0Runtime = BDMh<0>( mesh, RuntimeOrder{ 0 } );
    checkOnParity( bdm0Static, bdm0Runtime, mesh, expr, 1e-12 );
    checkBdmStaticDynamicParityAndExactness<2, 0>( mesh );
    checkBdmStaticDynamicParityAndExactness<2, 1>( mesh );
    checkBdmStaticDynamicParityAndExactness<2, 2>( mesh );
    checkBdmStaticDynamicParityAndExactness<2, 3>( mesh );
    checkBdmStaticDynamicParityAndExactness<3, 0>( mesh3d );
    checkBdmStaticDynamicParityAndExactness<3, 1>( mesh3d );
    checkBdmStaticDynamicParityAndExactness<3, 2>( mesh3d );
    checkBdmStaticDynamicParityAndExactness<3, 3>( mesh3d );
}

BOOST_AUTO_TEST_CASE( neh_runtime_order_constructor_parity )
{
    auto mesh = makeOneElementSimplexMesh<2>( "neh-runtime-parity" );
    auto expr = vec( -Py(), Px() );

    auto neh0Static = Neh<0>( mesh );
    auto neh0Runtime = Neh<0>( mesh, RuntimeOrder{ 0 } );
    checkOnParity( neh0Static, neh0Runtime, mesh, expr, 1e-12 );
}

BOOST_AUTO_TEST_SUITE_END()
