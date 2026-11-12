/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#define BOOST_TEST_MODULE product function spaces adapter
#include <feel/feelcore/testsuite.hpp>

#include <type_traits>

#include <boost/fusion/include/at_c.hpp>

#include <feel/feeldiscr/concepts.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/productfunctionspaces.hpp>
#include <feel/feeldiscr/thch.hpp>
#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;
using namespace boost::hana::literals;

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
}

BOOST_AUTO_TEST_SUITE( productfunctionspaces_adapter )

BOOST_AUTO_TEST_CASE( adapter_satisfies_product_and_composite_contracts )
{
    auto mesh = Feel::Test::makeMesh();
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = productFunctionSpaces( Vh, Qh );

    using adapter_type = std::remove_reference_t<decltype( Xh )>;
    static_assert( FunctionSpaceConcept<adapter_type> );
    static_assert( CompositeSpaceConcept<adapter_type> );
    static_assert( ProductBackedCompositeSpaceConcept<adapter_type> );
    static_assert( !LegacyCompositeFunctionSpaceConcept<adapter_type> );
    static_assert( ProductSpaceConcept<adapter_type> );
    static_assert( ProductSpacesConcept<adapter_type> );
    static_assert( adapter_type::is_composite );
    static_assert( !adapter_type::uses_internal_composite );
    static_assert( adapter_type::is_product_backed_composite );
    static_assert( adapter_type::nSpaces == 2 );
    static_assert( std::is_same_v<typename adapter_type::template sub_functionspace_type<0>,
                                  std::remove_reference_t<decltype( *Vh )>> );
    static_assert( std::is_same_v<typename adapter_type::template sub_functionspace_type<1>,
                                  std::remove_reference_t<decltype( *Qh )>> );

    BOOST_CHECK_EQUAL( Xh.numberOfSpaces(), 2 );
    BOOST_CHECK_EQUAL( Xh.nSubFunctionSpace(), 2 );
    BOOST_CHECK_EQUAL( Xh.template functionSpace<0>(), Vh );
    BOOST_CHECK_EQUAL( Xh.template functionSpace<1>(), Qh );
    BOOST_CHECK_EQUAL( Xh.template space<0>(), Vh );
    BOOST_CHECK_EQUAL( Xh.template space<1>(), Qh );
    BOOST_CHECK_EQUAL( Xh.space( 0_c ), Vh );
    BOOST_CHECK_EQUAL( Xh.space( 1_c ), Qh );
    BOOST_CHECK_EQUAL( boost::fusion::at_c<0>( Xh.functionSpaces() ), Vh );
    BOOST_CHECK_EQUAL( boost::fusion::at_c<1>( Xh.functionSpaces() ), Qh );
}

BOOST_AUTO_TEST_CASE( adapter_delegates_dof_accounting_to_product_spaces )
{
    auto mesh = Feel::Test::makeMesh();
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = productFunctionSpaces( product( Vh, Qh ) );

    BOOST_CHECK_EQUAL( Xh.nDof(), Vh->nDof() + Qh->nDof() );
    BOOST_CHECK_EQUAL( Xh.nLocalDof(), Vh->nLocalDof() + Qh->nLocalDof() );
    BOOST_CHECK_EQUAL( Xh.nLocalDofWithoutGhost(),
                       Vh->nLocalDofWithoutGhost() + Qh->nLocalDofWithoutGhost() );

    BOOST_CHECK_EQUAL( Xh.nDofStart( 0 ), 0 );
    BOOST_CHECK_EQUAL( Xh.nDofStart( 1 ), Vh->nDof() );
    BOOST_CHECK_EQUAL( Xh.nDofStart( 2 ), Xh.nDof() );
    BOOST_CHECK_EQUAL( Xh.nLocalDofStart( 1 ), Vh->nLocalDof() );
    BOOST_CHECK_EQUAL( Xh.nLocalDofWithoutGhostStart( 1 ), Vh->nLocalDofWithoutGhost() );
    BOOST_CHECK_EQUAL( Xh.blockDofStart( 1 ), Xh.nDofStart( 1 ) );
    BOOST_CHECK_EQUAL( Xh.blockLocalDofStart( 1 ), Xh.nLocalDofStart( 1 ) );
    BOOST_CHECK_EQUAL( Xh.blockMapPtr( 0 ), Vh->mapPtr() );
    BOOST_CHECK_EQUAL( Xh.blockMapPtr( 1 ), Qh->mapPtr() );

    BOOST_CHECK( Xh.mapPtr() );
    BOOST_CHECK_EQUAL( Xh.mapPtr()->nDof(), Xh.nDof() );
    BOOST_CHECK_EQUAL( Xh.mapPtr()->nLocalDofWithGhost(), Xh.nLocalDof() );
    BOOST_CHECK_EQUAL( Xh.dof()->nDof(), Xh.nDof() );
}

BOOST_AUTO_TEST_CASE( adapter_element_exposes_legacy_subaccess )
{
    auto mesh = Feel::Test::makeMesh();
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = productFunctionSpaces( Vh, Qh );
    auto U = Xh.element( "U" );

    using element_type = std::remove_reference_t<decltype( U )>;
    static_assert( FunctionSpaceElementConcept<element_type> );

    BOOST_CHECK_EQUAL( U.nDof(), Xh.nDof() );
    BOOST_CHECK_EQUAL( U.functionSpace().nDof(), Xh.nDof() );
    BOOST_CHECK_EQUAL( U.template functionSpace<0>(), Vh );
    BOOST_CHECK_EQUAL( U.template functionSpace<1>(), Qh );
    BOOST_CHECK_EQUAL( U.template element<0>().functionSpace(), Vh );
    BOOST_CHECK_EQUAL( U.template element<1>().functionSpace(), Qh );
    BOOST_CHECK_EQUAL( U.template element<0>().size(), Vh->nDof() );
    BOOST_CHECK_EQUAL( U.template element<1>().size(), Qh->nDof() );
    BOOST_CHECK_EQUAL( U[0_c].functionSpace(), Vh );
    BOOST_CHECK_EQUAL( U[1_c].functionSpace(), Qh );
}

BOOST_AUTO_TEST_CASE( adapter_matches_legacy_taylor_hood_layout )
{
    auto mesh = Feel::Test::makeMesh();
    auto legacy = THch<1>( mesh );
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = productFunctionSpaces( Vh, Qh );

    BOOST_CHECK_EQUAL( Xh.nDof(), legacy->nDof() );
    BOOST_CHECK_EQUAL( Xh.nLocalDof(), legacy->nLocalDof() );
    BOOST_CHECK_EQUAL( Xh.nLocalDofWithoutGhost(), legacy->nLocalDofWithoutGhost() );
    BOOST_CHECK_EQUAL( Xh.nDofStart( 1 ), legacy->nDofStart( 1 ) );
    BOOST_CHECK_EQUAL( Xh.nLocalDofStart( 1 ), legacy->nLocalDofStart( 1 ) );
    BOOST_CHECK_EQUAL( Xh.template functionSpace<0>()->nDof(),
                       legacy->template functionSpace<0>()->nDof() );
    BOOST_CHECK_EQUAL( Xh.template functionSpace<1>()->nDof(),
                       legacy->template functionSpace<1>()->nDof() );
}

BOOST_AUTO_TEST_SUITE_END()
