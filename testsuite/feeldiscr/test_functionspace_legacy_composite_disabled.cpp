/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#define BOOST_TEST_MODULE functionspace legacy composite disabled
#include <feel/feelcore/testsuite.hpp>

#include <type_traits>

#include <feel/feeldiscr/concepts.hpp>
#include <feel/feeldiscr/dhpdh.hpp>
#include <feel/feeldiscr/p2ch.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/productfunctionspaces.hpp>
#include <feel/feeldiscr/thch.hpp>
#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace Feel::Test
{
template<typename ProductSpaceType>
void checkProductBackedCompositeType()
{
    static_assert( FunctionSpaceConcept<ProductSpaceType> );
    static_assert( CompositeSpaceConcept<ProductSpaceType> );
    static_assert( ProductBackedCompositeSpaceConcept<ProductSpaceType> );
    static_assert( ProductSpaceConcept<ProductSpaceType> );
    static_assert( ProductSpacesConcept<ProductSpaceType> );
    static_assert( !LegacyCompositeFunctionSpaceConcept<ProductSpaceType> );
    static_assert( !ProductSpaceType::legacy_composite_enabled );
    static_assert( !ProductSpaceType::uses_internal_composite );
    static_assert( ProductSpaceType::is_product_backed_composite );
}

template<typename ProductSpacePtr>
void checkProductHelperLayout( ProductSpacePtr const& productSpace )
{
    BOOST_REQUIRE( productSpace );
    auto const Xh0 = productSpace->template functionSpace<0>();
    auto const Xh1 = productSpace->template functionSpace<1>();
    BOOST_REQUIRE( Xh0 );
    BOOST_REQUIRE( Xh1 );

    BOOST_CHECK_EQUAL( productSpace->numberOfSpaces(), 2 );
    BOOST_CHECK_EQUAL( productSpace->nSubFunctionSpace(), 2 );
    BOOST_CHECK_EQUAL( productSpace->template space<0>(), Xh0 );
    BOOST_CHECK_EQUAL( productSpace->template space<1>(), Xh1 );
    BOOST_CHECK_EQUAL( productSpace->nDof(), Xh0->nDof() + Xh1->nDof() );
    BOOST_CHECK_EQUAL( productSpace->nDofStart( 0 ), 0 );
    BOOST_CHECK_EQUAL( productSpace->nDofStart( 1 ), Xh0->nDof() );
    BOOST_CHECK_EQUAL( productSpace->nDofStart( 2 ), productSpace->nDof() );
    BOOST_CHECK_EQUAL( productSpace->nLocalDof(), Xh0->nLocalDof() + Xh1->nLocalDof() );
    BOOST_CHECK_EQUAL( productSpace->nLocalDofStart( 1 ), Xh0->nLocalDof() );
    BOOST_CHECK_EQUAL( productSpace->blockDofStart( 1 ), productSpace->nDofStart( 1 ) );
    BOOST_CHECK_EQUAL( productSpace->blockLocalDofStart( 1 ), productSpace->nLocalDofStart( 1 ) );
    BOOST_CHECK_EQUAL( productSpace->blockMapPtr( 0 ), Xh0->mapPtr() );
    BOOST_CHECK_EQUAL( productSpace->blockMapPtr( 1 ), Xh1->mapPtr() );
    BOOST_CHECK_EQUAL( productSpace->mapPtr()->nDof(), productSpace->nDof() );

    auto U = productSpace->element( "U" );
    BOOST_CHECK_EQUAL( U.template functionSpace<0>(), Xh0 );
    BOOST_CHECK_EQUAL( U.template functionSpace<1>(), Xh1 );
    BOOST_CHECK_EQUAL( U.template element<0>().functionSpace(), Xh0 );
    BOOST_CHECK_EQUAL( U.template element<1>().functionSpace(), Xh1 );
}
}

BOOST_AUTO_TEST_SUITE( functionspace_legacy_composite_disabled )

BOOST_AUTO_TEST_CASE( product_backed_composite_is_available_when_raw_composite_is_disabled )
{
    static_assert( FEELPP_ENABLE_LEGACY_COMPOSITE_FUNCTIONSPACE == 0 );

    using mesh_type = Mesh<Simplex<2>>;
    using scalar_space_type = FunctionSpace<mesh_type, bases<Lagrange<1, Scalar>>>;
    using vectorial_space_type = FunctionSpace<mesh_type, bases<Lagrange<1, Vectorial>>>;

    static_assert( FunctionSpaceConcept<scalar_space_type> );
    static_assert( FunctionSpaceConcept<vectorial_space_type> );
    static_assert( NonCompositeSpaceConcept<scalar_space_type> );
    static_assert( NonCompositeSpaceConcept<vectorial_space_type> );
    static_assert( !scalar_space_type::uses_internal_composite );
    static_assert( !vectorial_space_type::uses_internal_composite );
    static_assert( !scalar_space_type::is_product_backed_composite );

    auto mesh = loadMesh( _mesh = new mesh_type );
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = productFunctionSpaces( Vh, Qh );

    using product_backed_type = std::remove_reference_t<decltype( Xh )>;
    static_assert( ProductBackedCompositeSpaceConcept<product_backed_type> );
    static_assert( !LegacyCompositeFunctionSpaceConcept<product_backed_type> );
    static_assert( !product_backed_type::legacy_composite_enabled );
    static_assert( !product_backed_type::uses_internal_composite );
    static_assert( product_backed_type::is_product_backed_composite );

    BOOST_CHECK_EQUAL( Xh.nSubFunctionSpace(), 2 );
    BOOST_CHECK_EQUAL( Xh.template functionSpace<0>(), Vh );
    BOOST_CHECK_EQUAL( Xh.template functionSpace<1>(), Qh );
    BOOST_CHECK_EQUAL( Xh.nDof(), Vh->nDof() + Qh->nDof() );
    BOOST_CHECK_EQUAL( Xh.nDofStart( 1 ), Vh->nDof() );

    auto U = Xh.element( "U" );
    BOOST_CHECK_EQUAL( U.template functionSpace<0>(), Vh );
    BOOST_CHECK_EQUAL( U.template functionSpace<1>(), Qh );
    BOOST_CHECK_EQUAL( U.template element<0>().functionSpace(), Vh );
    BOOST_CHECK_EQUAL( U.template element<1>().functionSpace(), Qh );
}

BOOST_AUTO_TEST_CASE( product_helper_factories_work_when_raw_composite_is_disabled )
{
    static_assert( FEELPP_ENABLE_LEGACY_COMPOSITE_FUNCTIONSPACE == 0 );

    using mesh_type = Mesh<Simplex<2>>;
    using velocity_basis_type = Lagrange<1, Vectorial>;
    using pressure_basis_type = Lagrange<0, Scalar, Discontinuous>;
    using taylor_hood_product_type = THchProduct_type<1, mesh_type>;
    using p2_product_type = P2chProduct_type<velocity_basis_type, pressure_basis_type, mesh_type>;
    using dhpdh_static_product_type = DhPdhProduct_type<mesh_type, 0>;
    using dhpdh_dynamic_product_type = DhPdhProduct_type<mesh_type, Dynamic>;

    Feel::Test::checkProductBackedCompositeType<taylor_hood_product_type>();
    Feel::Test::checkProductBackedCompositeType<p2_product_type>();
    Feel::Test::checkProductBackedCompositeType<dhpdh_static_product_type>();
    Feel::Test::checkProductBackedCompositeType<dhpdh_dynamic_product_type>();

    auto mesh = loadMesh( _mesh = new mesh_type );

    auto taylorHood = THchProduct<1>( mesh );
    Feel::Test::checkProductHelperLayout( taylorHood );

    auto p2 = P2chProduct<velocity_basis_type, pressure_basis_type>( mesh );
    Feel::Test::checkProductHelperLayout( p2 );

    auto dhpdhStatic = DhPdhProduct<0>( mesh );
    Feel::Test::checkProductHelperLayout( dhpdhStatic );

    RuntimeOrder order{ 0 };
    auto dhpdhDynamic = DhPdhProduct<Dynamic>( mesh, order );
    Feel::Test::checkProductHelperLayout( dhpdhDynamic );
}

BOOST_AUTO_TEST_SUITE_END()
