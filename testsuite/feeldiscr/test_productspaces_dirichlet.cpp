/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Feel++ Consortium
 Date: 22 Apr 2026

 Copyright (C) 2026 Feel++ Consortium

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

#define BOOST_TEST_MODULE test_productspaces_dirichlet
#include <feel/feelcore/testsuite.hpp>

#include <algorithm>
#include <array>
#include <cmath>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

namespace
{
enum class MixedEntityCondition
{
    Elements,
    Faces,
    Edges,
    Points
};

struct DeferredDirichletStrategyCase
{
    char const* label;
    char const* type;
    double tolerance;
};

auto const deferred_dirichlet_strategies = std::array{
    DeferredDirichletStrategyCase{ .label = "elimination", .type = "elimination", .tolerance = 1e-10 },
    DeferredDirichletStrategyCase{ .label = "elimination_symmetric", .type = "elimination_symmetric", .tolerance = 1e-10 }
};

template<typename MeshPtrType, typename WElementType, typename TElementType, typename BlockFormType, typename BlockLinearFormType>
void assembleBlockformSystem( MeshPtrType const& mesh,
                              WElementType const& W,
                              TElementType const& T,
                              BlockFormType& a,
                              BlockLinearFormType& l )
{
    auto u = W( 0_c );
    auto p = W( 1_c );
    auto v = T( 0_c );
    auto q = T( 1_c );

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
    a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( p ), grad( q ) ) + idt( p ) * id( q ) );

    l( 0_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 1.0 ) * id( v ) );
    l( 1_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 2.0 ) * id( q ) );
}

template<typename MeshPtrType, typename WElementType, typename TElementType, typename BlockFormType, typename BlockLinearFormType>
void assembleStaticCondensationSystem( MeshPtrType const& mesh,
                                       WElementType const& W,
                                       TElementType const& T,
                                       BlockFormType& a,
                                       BlockLinearFormType& l )
{
    auto u = W( 0_c );
    auto alpha = W( 1_c );
    auto v = T( 0_c );
    auto beta = T( 1_c );

    a( 0_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
    a( 0_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=idt( alpha ) * id( v ) );
    a( 1_c, 0_c ) += integrate( _range=elements( mesh ),
                                _expr=idt( u ) * id( beta ) );
    a( 1_c, 1_c ) += integrate( _range=elements( mesh ),
                                _expr=cst( 2.0 ) * idt( alpha ) * id( beta ) );

    l( 0_c ) += integrate( _range=elements( mesh ),
                           _expr=cst( 1.0 ) * id( v ) );
}

template<typename BlockFormType, typename BlockLinearFormType, typename MeshType, typename ElementType, typename GElementType>
void applyMixedEntityConditions( BlockFormType& a,
                                 BlockLinearFormType& l,
                                 std::shared_ptr<MeshType> const& mesh,
                                 ElementType const& u,
                                 GElementType const& g,
                                 std::array<MixedEntityCondition, 4> const& order,
                                 char const* strategyType )
{
    for ( auto const condition : order )
    {
        switch ( condition )
        {
        case MixedEntityCondition::Elements:
            a.row( 0_c ) += on( _range=elements( mesh ),
                                _rhs=l( 0_c ),
                                _element=u,
                                _expr=idv( g ) * cst( 2.0 ),
                                _type=strategyType );
            break;
        case MixedEntityCondition::Faces:
            a.row( 0_c ) += on( _range=markedfaces( mesh, "S" ),
                                _rhs=l( 0_c ),
                                _element=u,
                                _expr=idv( g ) * cst( 12.0 ),
                                _type=strategyType );
            break;
        case MixedEntityCondition::Edges:
            a.row( 0_c ) += on( _range=markededges( mesh, "L" ),
                                _rhs=l( 0_c ),
                                _element=u,
                                _expr=idv( g ) * cst( 22.0 ),
                                _type=strategyType );
            break;
        case MixedEntityCondition::Points:
            a.row( 0_c ) += on( _range=markedpoints( mesh, "P" ),
                                _rhs=l( 0_c ),
                                _element=u,
                                _expr=idv( g ) * cst( 32.0 ),
                                _type=strategyType );
            break;
        }
    }
}

template<typename MeshType, typename FunctionSpaceType, typename SolutionElementType, typename GElementType>
void checkMixedEntityPrecedence( std::shared_ptr<MeshType> const& mesh,
                                 std::shared_ptr<FunctionSpaceType> const& Xh,
                                 SolutionElementType const& solution,
                                 GElementType const& g,
                                 double tolerance )
{
    auto err = Xh->element();

    err.on( _range=elements( mesh ), _expr=abs( idv( solution ) - idv( g ) * cst( 2.0 ) ) );
    auto dofsOnElement = err.functionSpace()->dofs( elements( mesh ), ComponentType::NO_COMPONENT, true );
    sync( err, "=", dofsOnElement );

    err.on( _range=markedfaces( mesh, "S" ), _expr=abs( idv( solution ) - idv( g ) * cst( 12.0 ) ) );
    auto dofsOnFace = err.functionSpace()->dofs( markedfaces( mesh, "S" ), ComponentType::NO_COMPONENT, true );
    sync( err, "=", dofsOnFace );

    err.on( _range=markededges( mesh, "L" ), _expr=abs( idv( solution ) - idv( g ) * cst( 22.0 ) ) );
    auto dofsOnEdge = err.functionSpace()->dofs( markededges( mesh, "L" ), ComponentType::NO_COMPONENT, true );
    sync( err, "=", dofsOnEdge );

    err.on( _range=markedpoints( mesh, "P" ), _expr=abs( idv( solution ) - idv( g ) * cst( 32.0 ) ) );
    auto dofsOnPoint = err.functionSpace()->dofs( markedpoints( mesh, "P" ), ComponentType::NO_COMPONENT, true );
    sync( err, "=", dofsOnPoint );

    BOOST_CHECK_SMALL( err.max(), tolerance );
}

void runBlockformMixedEntityOrderAllPermutations()
{
    using mesh_type = Mesh<Simplex<3>>;

    backend( _rebuild=true );

    auto mesh = loadMesh( _mesh=new mesh_type );
    BOOST_REQUIRE_GT( nelements( elements( mesh ), true ), 0 );
    BOOST_REQUIRE_GT( nelements( markedfaces( mesh, "S" ), true ), 0 );
    BOOST_REQUIRE_GT( nelements( markededges( mesh, "L" ), true ), 0 );
    BOOST_REQUIRE_GT( nelements( markedpoints( mesh, "P" ), true ), 0 );

    auto Xh = Pch<1>( mesh );
    auto Yh = Pch<1>( mesh );
    auto ps = product( Xh, Yh );
    auto g = Xh->element();
    g.on( _range=elements( mesh ), _expr=Px() + Py() + Pz() );

    auto order = std::array{
        MixedEntityCondition::Elements,
        MixedEntityCondition::Faces,
        MixedEntityCondition::Edges,
        MixedEntityCondition::Points
    };
    double const parallelToleranceFloor = 2e-9;

    for ( auto const& strategy : deferred_dirichlet_strategies )
    {
        double const tolerance = Environment::worldComm().globalSize() > 1 ?
            std::max( strategy.tolerance, parallelToleranceFloor ) :
            strategy.tolerance;
        auto permutation = order;
        int permutationCount = 0;
        do
        {
            BOOST_TEST_CONTEXT( "strategy=" << strategy.label << ", permutation=" << permutationCount )
            {
                auto W = ps.element();
                auto T = ps.element();
                auto u = W( 0_c );

                auto a = blockform2( ps, solve::strategy::monolithic, backend() );
                auto l = blockform1( ps, solve::strategy::monolithic, backend() );
                BOOST_CHECK( a.useDeferredDirichlet() );

                assembleBlockformSystem( mesh, W, T, a, l );
                l.close();
                a.close();

                applyMixedEntityConditions( a, l, mesh, u, g, permutation, strategy.type );

                BOOST_CHECK( a.hasPendingDirichletConstraints() );
                BOOST_CHECK( a.hasDirichletConstraints() );

                auto solution = ps.element();
                BOOST_REQUIRE_NO_THROW( a.solve( _solution=solution, _rhs=l ) );

                BOOST_CHECK( !a.hasPendingDirichletConstraints() );
                BOOST_CHECK( a.hasDirichletConstraints() );
                checkMixedEntityPrecedence( mesh, Xh, solution( 0_c ), g, tolerance );

                double const secondFieldNorm = normL2( _range=elements( mesh ),
                                                       _expr=idv( solution( 1_c ) ) );
                BOOST_CHECK( std::isfinite( secondFieldNorm ) );
            }
            ++permutationCount;
        }
        while ( std::next_permutation( permutation.begin(), permutation.end() ) );

        BOOST_CHECK_EQUAL( permutationCount, 24 );
    }
}

void runBlockformStaticCondensationMixedEntityOrderAllPermutations()
{
    using mesh_type = Mesh<Simplex<3>>;

    backend( _rebuild=true );

    auto mesh = loadMesh( _mesh=new mesh_type );
    BOOST_REQUIRE_GT( nelements( elements( mesh ), true ), 0 );
    BOOST_REQUIRE_GT( nelements( markedfaces( mesh, "S" ), true ), 0 );
    BOOST_REQUIRE_GT( nelements( markededges( mesh, "L" ), true ), 0 );
    BOOST_REQUIRE_GT( nelements( markedpoints( mesh, "P" ), true ), 0 );

    auto Uh = Pch<1>( mesh );
    auto Ah = Pdh<0>( mesh );
    auto ps = product( Uh, Ah );
    auto g = Uh->element();
    g.on( _range=elements( mesh ), _expr=Px() + Py() + Pz() );

    auto order = std::array{
        MixedEntityCondition::Elements,
        MixedEntityCondition::Faces,
        MixedEntityCondition::Edges,
        MixedEntityCondition::Points
    };
    double const parallelToleranceFloor = 2e-9;

    for ( auto const& strategy : deferred_dirichlet_strategies )
    {
        double const tolerance = Environment::worldComm().globalSize() > 1 ?
            std::max( strategy.tolerance, parallelToleranceFloor ) :
            strategy.tolerance;
        auto permutation = order;
        int permutationCount = 0;
        do
        {
            BOOST_TEST_CONTEXT( "strategy=" << strategy.label << ", permutation=" << permutationCount )
            {
                auto W = ps.element();
                auto T = ps.element();
                auto u = W( 0_c );

                auto a = blockform2( ps, solve::strategy::static_condensation, backend() );
                auto l = blockform1( ps, solve::strategy::static_condensation, backend() );
                BOOST_CHECK( a.useDeferredDirichlet() );

                assembleStaticCondensationSystem( mesh, W, T, a, l );
                l.close();
                a.close();

                applyMixedEntityConditions( a, l, mesh, u, g, permutation, strategy.type );

                BOOST_CHECK( a.hasPendingDirichletConstraints() );
                BOOST_CHECK( a.hasDirichletConstraints() );

                auto solution = ps.element();
                BOOST_REQUIRE_NO_THROW( a.solve( _solution=solution, _rhs=l,
                                                 _condense=true, _condenser=condenser_sb9() ) );

                BOOST_CHECK( !a.hasPendingDirichletConstraints() );
                BOOST_CHECK( a.hasDirichletConstraints() );
                checkMixedEntityPrecedence( mesh, Uh, solution( 0_c ), g, tolerance );

                double const secondFieldNorm = normL2( _range=elements( mesh ),
                                                       _expr=idv( solution( 1_c ) ) );
                BOOST_CHECK( std::isfinite( secondFieldNorm ) );
            }
            ++permutationCount;
        }
        while ( std::next_permutation( permutation.begin(), permutation.end() ) );

        BOOST_CHECK_EQUAL( permutationCount, 24 );
    }
}
} // namespace

FEELPP_ENVIRONMENT_WITH_ABOUT_NO_OPTIONS( Feel::makeAboutDefault( "test_productspaces_dirichlet" ) )

BOOST_AUTO_TEST_SUITE( productspaces_dirichlet_suite )

BOOST_AUTO_TEST_CASE( test_row_dirichlet_mixed_entity_order_all_permutations )
{
    runBlockformMixedEntityOrderAllPermutations();
}

BOOST_AUTO_TEST_CASE( test_row_dirichlet_static_condensation_mixed_entity_order_all_permutations )
{
    runBlockformStaticCondensationMixedEntityOrderAllPermutations();
}

BOOST_AUTO_TEST_SUITE_END()
