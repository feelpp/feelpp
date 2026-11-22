/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#define BOOST_TEST_MODULE functionspace composite fieldsplit
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/productfunctionspaces.hpp>
#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace Feel::Test
{
using mesh_type = Mesh<Simplex<2>>;
using mesh_ptrtype = std::shared_ptr<mesh_type>;

inline mesh_ptrtype
makeMesh()
{
    return loadMesh( _mesh = new mesh_type );
}

inline void
checkIndexSplitEqual( std::shared_ptr<IndexSplit> const& lhs,
                      std::shared_ptr<IndexSplit> const& rhs )
{
    BOOST_REQUIRE( lhs );
    BOOST_REQUIRE( rhs );
    BOOST_REQUIRE_EQUAL( lhs->size(), rhs->size() );

    for ( int k = 0; k < lhs->size(); ++k )
    {
        BOOST_CHECK_EQUAL( lhs->tag( k ), rhs->tag( k ) );
        BOOST_CHECK_EQUAL( lhs->firstIndex( k ), rhs->firstIndex( k ) );
        BOOST_CHECK_EQUAL( lhs->lastIndex( k ), rhs->lastIndex( k ) );
        BOOST_CHECK_EQUAL( lhs->nIndex( k ), rhs->nIndex( k ) );
        BOOST_CHECK_EQUAL( lhs->nIndexForSmallerRankId( k ), rhs->nIndexForSmallerRankId( k ) );
        BOOST_REQUIRE_EQUAL( lhs->split( k ).size(), rhs->split( k ).size() );

        for ( size_type i = 0; i < lhs->split( k ).size(); ++i )
            BOOST_CHECK_EQUAL( lhs->split( k )[i], rhs->split( k )[i] );
    }
}
}

BOOST_AUTO_TEST_SUITE( functionspace_composite_fieldsplit )

BOOST_AUTO_TEST_CASE( legacy_composite_owns_field_split )
{
    using mesh_type = Feel::Test::mesh_type;
    using legacy_space_type =
        FunctionSpace<mesh_type,
                      bases<Lagrange<2, Vectorial>,
                            Lagrange<1, Scalar>>>;

    auto mesh = Feel::Test::makeMesh();
    auto Xh = legacy_space_type::New( mesh );
    auto Vh = Xh->template functionSpace<0>();
    auto Qh = Xh->template functionSpace<1>();

    auto const& split = Xh->dof()->indexSplit();
    BOOST_REQUIRE( split );
    BOOST_REQUIRE_EQUAL( split->size(), 2 );

    BOOST_CHECK_EQUAL( split->nIndex( 0 ), Vh->dof()->nLocalDofWithoutGhost() );
    BOOST_CHECK_EQUAL( split->nIndex( 1 ), Qh->dof()->nLocalDofWithoutGhost() );
    BOOST_CHECK_EQUAL( split->split( 0 ).size(), Vh->dof()->nLocalDofWithoutGhost() );
    BOOST_CHECK_EQUAL( split->split( 1 ).size(), Qh->dof()->nLocalDofWithoutGhost() );
    BOOST_CHECK_EQUAL( split->firstIndex( 0 ), Xh->dof()->firstDofGlobalCluster() );
    BOOST_CHECK_EQUAL( split->firstIndex( 1 ),
                       Xh->dof()->firstDofGlobalCluster() + Vh->dof()->nLocalDofWithoutGhost() );

    Feel::Test::checkIndexSplitEqual( split, Xh->buildDofIndexSplit() );
}

BOOST_AUTO_TEST_CASE( legacy_composite_component_split_preserves_subspace_components )
{
    using mesh_type = Feel::Test::mesh_type;
    using legacy_space_type =
        FunctionSpace<mesh_type,
                      bases<Lagrange<2, Vectorial>,
                            Lagrange<1, Scalar>>>;

    auto mesh = Feel::Test::makeMesh();
    auto Xh = legacy_space_type::New( mesh );
    auto Vh = Xh->template functionSpace<0>();
    auto Qh = Xh->template functionSpace<1>();

    BOOST_REQUIRE( Xh->dof()->hasIndexSplitWithComponents() );

    auto const& splitWithComponents = Xh->dof()->indexSplitWithComponents();
    auto const& vectorComponents = Vh->dof()->indexSplitWithComponents();
    auto const& scalarSplit = Qh->dof()->indexSplit();

    BOOST_REQUIRE( splitWithComponents );
    BOOST_REQUIRE( vectorComponents );
    BOOST_REQUIRE( scalarSplit );
    BOOST_REQUIRE_EQUAL( vectorComponents->size(), mesh_type::nDim );
    BOOST_REQUIRE_EQUAL( scalarSplit->size(), 1 );
    BOOST_REQUIRE_EQUAL( splitWithComponents->size(), vectorComponents->size() + scalarSplit->size() );

    for ( int c = 0; c < vectorComponents->size(); ++c )
    {
        BOOST_CHECK_EQUAL( splitWithComponents->tag( c ), 0 );
        BOOST_CHECK_EQUAL( splitWithComponents->nIndex( c ), vectorComponents->nIndex( c ) );
        BOOST_CHECK_EQUAL( splitWithComponents->split( c ).size(), vectorComponents->split( c ).size() );
    }

    int const scalarSplitId = vectorComponents->size();
    BOOST_CHECK_EQUAL( splitWithComponents->tag( scalarSplitId ), 1 );
    BOOST_CHECK_EQUAL( splitWithComponents->nIndex( scalarSplitId ), scalarSplit->nIndex( 0 ) );
    BOOST_CHECK_EQUAL( splitWithComponents->split( scalarSplitId ).size(), scalarSplit->split( 0 ).size() );

    Feel::Test::checkIndexSplitEqual( splitWithComponents, Xh->buildDofIndexSplitWithComponents() );
}

BOOST_AUTO_TEST_CASE( scalar_legacy_composite_keeps_component_split_absent )
{
    using mesh_type = Feel::Test::mesh_type;
    using legacy_space_type =
        FunctionSpace<mesh_type,
                      bases<Lagrange<1, Scalar>,
                            Lagrange<0, Scalar, Discontinuous>>>;

    auto mesh = Feel::Test::makeMesh();
    auto Xh = legacy_space_type::New( mesh );

    BOOST_REQUIRE( Xh->dof()->indexSplit() );
    BOOST_CHECK( !Xh->dof()->hasIndexSplitWithComponents() );
    BOOST_CHECK_EQUAL( Xh->dof()->indexSplitWithComponents().get(), Xh->dof()->indexSplit().get() );
    BOOST_CHECK( !Xh->buildDofIndexSplitWithComponents() );
}

BOOST_AUTO_TEST_CASE( product_backed_facade_uses_same_composite_split_contract )
{
    auto mesh = Feel::Test::makeMesh();
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = productFunctionSpaces( Vh, Qh );

    auto const& split = Xh.dof()->indexSplit();
    BOOST_REQUIRE( split );
    BOOST_REQUIRE_EQUAL( split->size(), 2 );
    BOOST_CHECK_EQUAL( split->nIndex( 0 ), Vh->dof()->nLocalDofWithoutGhost() );
    BOOST_CHECK_EQUAL( split->nIndex( 1 ), Qh->dof()->nLocalDofWithoutGhost() );

    BOOST_REQUIRE( Xh.dof()->hasIndexSplitWithComponents() );
    auto const& splitWithComponents = Xh.dof()->indexSplitWithComponents();
    BOOST_REQUIRE( splitWithComponents );
    BOOST_CHECK_EQUAL( splitWithComponents->size(), Feel::Test::mesh_type::nDim + 1 );
}

BOOST_AUTO_TEST_SUITE_END()
