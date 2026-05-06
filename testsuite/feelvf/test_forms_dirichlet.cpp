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

#define BOOST_TEST_MODULE test_forms_dirichlet
#include <feel/feelcore/testsuite.hpp>

#include <algorithm>
#include <array>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

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

template<typename MeshPtrType, typename UElementType, typename VElementType, typename BilinearFormType, typename LinearFormType>
void assembleScalarForm2System( MeshPtrType const& mesh,
                                UElementType const& u,
                                VElementType const& v,
                                BilinearFormType& a,
                                LinearFormType& l )
{
    a += integrate( _range=elements( mesh ),
                    _expr=inner( gradt( u ), grad( v ) ) + idt( u ) * id( v ) );
    l += integrate( _range=elements( mesh ),
                    _expr=cst( 1.0 ) * id( v ) );
}

template<typename BilinearFormType, typename LinearFormType, typename MeshType, typename ElementType>
void applyMixedEntityConditions( BilinearFormType& a,
                                 LinearFormType& l,
                                 std::shared_ptr<MeshType> const& mesh,
                                 ElementType const& u,
                                 ElementType const& g,
                                 std::array<MixedEntityCondition, 4> const& order,
                                 char const* strategyType )
{
    for ( auto const condition : order )
    {
        switch ( condition )
        {
        case MixedEntityCondition::Elements:
            a += on( _range=elements( mesh ),
                     _rhs=l,
                     _element=u,
                     _expr=idv( g ) * cst( 2.0 ),
                     _type=strategyType );
            break;
        case MixedEntityCondition::Faces:
            a += on( _range=markedfaces( mesh, "S" ),
                     _rhs=l,
                     _element=u,
                     _expr=idv( g ) * cst( 12.0 ),
                     _type=strategyType );
            break;
        case MixedEntityCondition::Edges:
            a += on( _range=markededges( mesh, "L" ),
                     _rhs=l,
                     _element=u,
                     _expr=idv( g ) * cst( 22.0 ),
                     _type=strategyType );
            break;
        case MixedEntityCondition::Points:
            a += on( _range=markedpoints( mesh, "P" ),
                     _rhs=l,
                     _element=u,
                     _expr=idv( g ) * cst( 32.0 ),
                     _type=strategyType );
            break;
        }
    }
}

template<typename MeshType, typename FunctionSpaceType, typename ElementType>
void checkMixedEntityPrecedence( std::shared_ptr<MeshType> const& mesh,
                                 std::shared_ptr<FunctionSpaceType> const& Xh,
                                 ElementType const& solution,
                                 ElementType const& g,
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

void runDeferredDirichletMixedEntityOrderAllPermutations()
{
    using mesh_type = Mesh<Simplex<3>>;

    backend( _rebuild=true );

    auto mesh = loadMesh( _mesh=new mesh_type );
    BOOST_REQUIRE_GT( nelements( elements( mesh ), true ), 0 );
    BOOST_REQUIRE_GT( nelements( markedfaces( mesh, "S" ), true ), 0 );
    BOOST_REQUIRE_GT( nelements( markededges( mesh, "L" ), true ), 0 );
    BOOST_REQUIRE_GT( nelements( markedpoints( mesh, "P" ), true ), 0 );

    auto Xh = Pch<1>( mesh );
    auto g = Xh->element();
    g.on( _range=elements( mesh ), _expr=Px() + Py() + Pz() );

    auto order = std::array{
        MixedEntityCondition::Elements,
        MixedEntityCondition::Faces,
        MixedEntityCondition::Edges,
        MixedEntityCondition::Points
    };

    for ( auto const& strategy : deferred_dirichlet_strategies )
    {
        auto permutation = order;
        int permutationCount = 0;
        do
        {
            BOOST_TEST_CONTEXT( "strategy=" << strategy.label << ", permutation=" << permutationCount )
            {
                auto u = Xh->element();
                auto v = Xh->element();
                auto a = form2( _trial=Xh, _test=Xh );
                auto l = form1( _test=Xh );

                assembleScalarForm2System( mesh, u, v, a, l );
                BOOST_CHECK( a.useDeferredDirichlet() );

                applyMixedEntityConditions( a, l, mesh, u, g, permutation, strategy.type );

                BOOST_CHECK( a.hasPendingDirichletConstraints() );
                BOOST_CHECK( a.hasDirichletConstraints() );

                BOOST_REQUIRE_NO_THROW( a.solve( _solution=u, _rhs=l ) );

                BOOST_CHECK( !a.hasPendingDirichletConstraints() );
                BOOST_CHECK( a.hasDirichletConstraints() );
                checkMixedEntityPrecedence( mesh, Xh, u, g, strategy.tolerance );
            }
            ++permutationCount;
        }
        while ( std::next_permutation( permutation.begin(), permutation.end() ) );

        BOOST_CHECK_EQUAL( permutationCount, 24 );
    }
}
} // namespace

FEELPP_ENVIRONMENT_WITH_ABOUT_NO_OPTIONS( Feel::makeAboutDefault( "test_forms_dirichlet" ) )

BOOST_AUTO_TEST_SUITE( forms_dirichlet_suite )

BOOST_AUTO_TEST_CASE( test_form2_deferred_dirichlet_mixed_entity_order_all_permutations )
{
    runDeferredDirichletMixedEntityOrderAllPermutations();
}

BOOST_AUTO_TEST_SUITE_END()
