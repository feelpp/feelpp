/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#define BOOST_TEST_MODULE mixed space product helpers
#include <feel/feelcore/testsuite.hpp>

#include <type_traits>

#include <feel/feeldiscr/concepts.hpp>
#include <feel/feeldiscr/dhpdh.hpp>
#include <feel/feeldiscr/p2ch.hpp>
#include <feel/feeldiscr/thch.hpp>
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

template<typename ProductSpacePtr, typename LegacySpacePtr>
void checkProductHelperLayout( ProductSpacePtr const& productSpace, LegacySpacePtr const& legacySpace )
{
    BOOST_CHECK_EQUAL( productSpace->numberOfSpaces(), legacySpace->nSubFunctionSpace() );
    BOOST_CHECK_EQUAL( productSpace->nSubFunctionSpace(), legacySpace->nSubFunctionSpace() );
    BOOST_CHECK_EQUAL( productSpace->nDof(), legacySpace->nDof() );
    BOOST_CHECK_EQUAL( productSpace->nLocalDof(), legacySpace->nLocalDof() );
    BOOST_CHECK_EQUAL( productSpace->nLocalDofWithoutGhost(), legacySpace->nLocalDofWithoutGhost() );
    BOOST_CHECK_EQUAL( productSpace->nDofStart( 0 ), legacySpace->nDofStart( 0 ) );
    BOOST_CHECK_EQUAL( productSpace->nDofStart( 1 ), legacySpace->nDofStart( 1 ) );
    BOOST_CHECK_EQUAL( productSpace->nDofStart( 2 ), legacySpace->nDofStart( 2 ) );
    BOOST_CHECK_EQUAL( productSpace->nLocalDofStart( 1 ), legacySpace->nLocalDofStart( 1 ) );
    BOOST_CHECK_EQUAL( productSpace->template functionSpace<0>()->nDof(),
                       legacySpace->template functionSpace<0>()->nDof() );
    BOOST_CHECK_EQUAL( productSpace->template functionSpace<1>()->nDof(),
                       legacySpace->template functionSpace<1>()->nDof() );
    BOOST_CHECK_EQUAL( productSpace->blockMapPtr( 0 ), productSpace->template functionSpace<0>()->mapPtr() );
    BOOST_CHECK_EQUAL( productSpace->blockMapPtr( 1 ), productSpace->template functionSpace<1>()->mapPtr() );
    BOOST_CHECK_EQUAL( productSpace->mapPtr()->nDof(), productSpace->nDof() );

    auto U = productSpace->element( "U" );
    BOOST_CHECK_EQUAL( U.template functionSpace<0>(), productSpace->template functionSpace<0>() );
    BOOST_CHECK_EQUAL( U.template functionSpace<1>(), productSpace->template functionSpace<1>() );
    BOOST_CHECK_EQUAL( U.template element<0>().functionSpace(), productSpace->template functionSpace<0>() );
    BOOST_CHECK_EQUAL( U.template element<1>().functionSpace(), productSpace->template functionSpace<1>() );
}
}

BOOST_AUTO_TEST_SUITE( mixedspace_product_helpers )

BOOST_AUTO_TEST_CASE( taylor_hood_product_helper_matches_legacy_helper )
{
    using mesh_type = Feel::Test::mesh_type;
    using legacy_type = THch_type<1,mesh_type>;
    using product_type = THch_product_type<1,mesh_type>;

    static_assert( CompositeSpaceConcept<legacy_type> );
    static_assert( LegacyCompositeFunctionSpaceConcept<legacy_type> );
    static_assert( FunctionSpaceConcept<product_type> );
    static_assert( CompositeSpaceConcept<product_type> );
    static_assert( ProductBackedCompositeSpaceConcept<product_type> );
    static_assert( !LegacyCompositeFunctionSpaceConcept<product_type> );
    static_assert( ProductSpaceConcept<product_type> );
    static_assert( ProductSpacesConcept<product_type> );
    static_assert( !std::is_same_v<legacy_type, product_type> );

    auto mesh = Feel::Test::makeMesh();
    auto legacy = THch<1>( mesh );
    auto productSpace = THchProduct<1>( mesh );

    Feel::Test::checkProductHelperLayout( productSpace, legacy );
}

BOOST_AUTO_TEST_CASE( p2_product_helper_matches_legacy_helper )
{
    using mesh_type = Feel::Test::mesh_type;
    using velocity_basis_type = Lagrange<1, Vectorial>;
    using pressure_basis_type = Lagrange<0, Scalar, Discontinuous>;
    using legacy_type = P2ch_type<velocity_basis_type, pressure_basis_type, mesh_type>;
    using product_type = P2ch_product_type<velocity_basis_type, pressure_basis_type, mesh_type>;

    static_assert( CompositeSpaceConcept<legacy_type> );
    static_assert( LegacyCompositeFunctionSpaceConcept<legacy_type> );
    static_assert( FunctionSpaceConcept<product_type> );
    static_assert( CompositeSpaceConcept<product_type> );
    static_assert( ProductBackedCompositeSpaceConcept<product_type> );
    static_assert( !LegacyCompositeFunctionSpaceConcept<product_type> );
    static_assert( ProductSpaceConcept<product_type> );
    static_assert( ProductSpacesConcept<product_type> );
    static_assert( !std::is_same_v<legacy_type, product_type> );

    auto mesh = Feel::Test::makeMesh();
    auto legacy = P2ch<velocity_basis_type, pressure_basis_type>( mesh );
    auto productSpace = P2chProduct<velocity_basis_type, pressure_basis_type>( mesh );

    Feel::Test::checkProductHelperLayout( productSpace, legacy );
}

BOOST_AUTO_TEST_CASE( dhpdh_static_product_helper_matches_legacy_helper )
{
    using mesh_type = Feel::Test::mesh_type;
    using legacy_type = DhPdh_type<mesh_type,0>;
    using product_type = DhPdh_product_type<mesh_type,0>;

    static_assert( CompositeSpaceConcept<legacy_type> );
    static_assert( LegacyCompositeFunctionSpaceConcept<legacy_type> );
    static_assert( FunctionSpaceConcept<product_type> );
    static_assert( CompositeSpaceConcept<product_type> );
    static_assert( ProductBackedCompositeSpaceConcept<product_type> );
    static_assert( !LegacyCompositeFunctionSpaceConcept<product_type> );
    static_assert( ProductSpaceConcept<product_type> );
    static_assert( ProductSpacesConcept<product_type> );
    static_assert( !std::is_same_v<legacy_type, product_type> );

    auto mesh = Feel::Test::makeMesh();
    auto legacy = DhPdh<0>( mesh );
    auto productSpace = DhPdhProduct<0>( mesh );

    Feel::Test::checkProductHelperLayout( productSpace, legacy );
}

BOOST_AUTO_TEST_CASE( dhpdh_dynamic_product_helper_matches_legacy_helper )
{
    using mesh_type = Feel::Test::mesh_type;
    using legacy_type = DhPdh_type<mesh_type,Dynamic>;
    using product_type = DhPdh_product_type<mesh_type,Dynamic>;

    static_assert( CompositeSpaceConcept<legacy_type> );
    static_assert( LegacyCompositeFunctionSpaceConcept<legacy_type> );
    static_assert( FunctionSpaceConcept<product_type> );
    static_assert( CompositeSpaceConcept<product_type> );
    static_assert( ProductBackedCompositeSpaceConcept<product_type> );
    static_assert( !LegacyCompositeFunctionSpaceConcept<product_type> );
    static_assert( ProductSpaceConcept<product_type> );
    static_assert( ProductSpacesConcept<product_type> );
    static_assert( !std::is_same_v<legacy_type, product_type> );

    auto mesh = Feel::Test::makeMesh();
    RuntimeOrder order{ 0 };
    auto legacy = DhPdh<Dynamic>( mesh, order );
    auto productSpace = DhPdhProduct<Dynamic>( mesh, order );

    Feel::Test::checkProductHelperLayout( productSpace, legacy );
}

BOOST_AUTO_TEST_SUITE_END()
