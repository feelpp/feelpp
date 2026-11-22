/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#define BOOST_TEST_MODULE functionspace composite guardrails
#include <feel/feelcore/testsuite.hpp>

#include <type_traits>

#include <feel/feeldiscr/concepts.hpp>
#include <feel/feeldiscr/dhpdh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/moch.hpp>
#include <feel/feeldiscr/p2ch.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
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

template<typename SpaceType>
concept HasLegacyFunctionSpacesAccessor = requires( SpaceType const& space )
{
    space.functionSpaces();
};

template<typename SpaceType>
concept HasLegacySubFunctionSpaceAccessor = requires( SpaceType& space )
{
    space.template functionSpace<0>();
};

template<typename SpaceType>
concept HasLegacyBuildDofIndexSplitAccessor = requires( SpaceType& space )
{
    space.buildDofIndexSplit();
};

template<typename ElementType>
concept HasLegacyElementSubFunctionSpaceAccessor = requires( ElementType const& element )
{
    element.template functionSpace<0>();
};
}

BOOST_AUTO_TEST_SUITE( functionspace_composite_guardrails )

BOOST_AUTO_TEST_CASE( single_space_contracts_use_existing_concepts )
{
    using mesh_type = Feel::Test::mesh_type;
    using scalar_space_type = FunctionSpace<mesh_type, bases<Lagrange<1, Scalar>>>;
    using vectorial_space_type = FunctionSpace<mesh_type, bases<Lagrange<1, Vectorial>>>;

    static_assert( FunctionSpaceConcept<scalar_space_type> );
    static_assert( NonCompositeSpaceConcept<scalar_space_type> );
    static_assert( !CompositeSpaceConcept<scalar_space_type> );
    static_assert( FunctionSpaceConcept<vectorial_space_type> );
    static_assert( NonCompositeSpaceConcept<vectorial_space_type> );
    static_assert( !vectorial_space_type::is_composite );
    static_assert( vectorial_space_type::is_product );
    static_assert( !scalar_space_type::has_type_level_periodicity );
    static_assert( !scalar_space_type::is_periodic );
    static_assert( !scalar_space_type::is_mortar );
    static_assert( std::is_empty_v<typename scalar_space_type::legacy_composite_storage_type> );
    static_assert( std::is_same_v<typename scalar_space_type::periodicity_0_type, NoPeriodicity> );
    static_assert( !Feel::Test::HasLegacyFunctionSpacesAccessor<scalar_space_type> );
    static_assert( !Feel::Test::HasLegacySubFunctionSpaceAccessor<scalar_space_type> );
    static_assert( !Feel::Test::HasLegacyBuildDofIndexSplitAccessor<scalar_space_type> );
    static_assert( !Feel::Test::HasLegacyElementSubFunctionSpaceAccessor<typename scalar_space_type::element_type> );

    auto mesh = Feel::Test::makeMesh();
    auto Xh = Pch<1>( mesh );
    auto Vh = Pchv<1>( mesh );

    static_assert( FunctionSpacePtrConcept<decltype( Xh )> );
    static_assert( FunctionSpacePtrConcept<decltype( Vh )> );

    BOOST_CHECK_EQUAL( Xh->nSubFunctionSpace(), 1 );
    BOOST_CHECK_EQUAL( Xh->nDof(), Xh->dof()->nDof() );
    BOOST_CHECK_EQUAL( Xh->nLocalDof(), Xh->dof()->nLocalDofWithGhost() );
    BOOST_CHECK_EQUAL( Xh->nLocalDofWithoutGhost(), Xh->dof()->nLocalDofWithoutGhost() );
    BOOST_CHECK_EQUAL( Xh->nDofStart( 0 ), 0 );
    BOOST_CHECK_EQUAL( Xh->nDofStart( 1 ), Xh->nDof() );
    BOOST_CHECK_EQUAL( Xh->nLocalDofStart( 0 ), 0 );
    BOOST_CHECK_EQUAL( Xh->nLocalDofStart( 1 ), Xh->nLocalDof() );
    BOOST_CHECK_EQUAL( Xh->nLocalDofWithGhostStart( 1 ), Xh->nLocalDofWithGhost() );
    BOOST_CHECK_EQUAL( Xh->nLocalDofWithoutGhostStart( 1 ), Xh->nLocalDofWithoutGhost() );
    BOOST_CHECK( !Xh->isMortar() );
    BOOST_CHECK( !Xh->meshHasPeriodicity() );
    BOOST_CHECK( !boost::fusion::at_c<0>( Xh->periodicity() ).isPeriodic() );

    BOOST_CHECK_EQUAL( Vh->nSubFunctionSpace(), 1 );
    BOOST_CHECK_EQUAL( Vh->nDof(), Vh->dof()->nDof() );
    BOOST_CHECK_EQUAL( Vh->nLocalDof(), Vh->dof()->nLocalDofWithGhost() );
    BOOST_CHECK_EQUAL( Vh->qDim(), mesh_type::nDim );
}

BOOST_AUTO_TEST_CASE( legacy_composite_space_offsets_and_dof_aggregation )
{
    using mesh_type = Feel::Test::mesh_type;
    using legacy_space_type =
        FunctionSpace<mesh_type,
                      bases<Lagrange<2, Vectorial>,
                            Lagrange<1, Scalar>,
                            Lagrange<0, Scalar, Discontinuous>>>;

    static_assert( FunctionSpaceConcept<legacy_space_type> );
    static_assert( CompositeSpaceConcept<legacy_space_type> );
    static_assert( LegacyCompositeFunctionSpaceConcept<legacy_space_type> );
    static_assert( !ProductBackedCompositeSpaceConcept<legacy_space_type> );
    static_assert( !NonCompositeSpaceConcept<legacy_space_type> );
    static_assert( legacy_space_type::is_composite );
    static_assert( legacy_space_type::is_legacy_composite );
    static_assert( legacy_space_type::uses_internal_composite );
    static_assert( !legacy_space_type::is_product_backed_composite );
    static_assert( legacy_space_type::legacy_composite_enabled );
    static_assert( legacy_space_type::nSpaces == 3 );
    static_assert( std::is_same_v<typename legacy_space_type::dof_type, DofComposite> );
    static_assert( Feel::Test::HasLegacyFunctionSpacesAccessor<legacy_space_type> );
    static_assert( Feel::Test::HasLegacySubFunctionSpaceAccessor<legacy_space_type> );
    static_assert( Feel::Test::HasLegacyBuildDofIndexSplitAccessor<legacy_space_type> );
    static_assert( Feel::Test::HasLegacyElementSubFunctionSpaceAccessor<typename legacy_space_type::element_type> );

    auto mesh = Feel::Test::makeMesh();
    auto Xh = legacy_space_type::New( mesh );
    auto U = Xh->element( "U" );

    auto const n0 = Xh->template functionSpace<0>()->nDof();
    auto const n1 = Xh->template functionSpace<1>()->nDof();
    auto const n2 = Xh->template functionSpace<2>()->nDof();

    BOOST_CHECK_EQUAL( Xh->nSubFunctionSpace(), 3 );
    BOOST_CHECK_EQUAL( Xh->nDof(), n0 + n1 + n2 );
    BOOST_CHECK_EQUAL( Xh->dof()->nDof(), Xh->nDof() );
    BOOST_CHECK_EQUAL( U.nDof(), Xh->nDof() );

    BOOST_CHECK_EQUAL( Xh->nLocalDof(),
                       Xh->template functionSpace<0>()->nLocalDof() +
                           Xh->template functionSpace<1>()->nLocalDof() +
                           Xh->template functionSpace<2>()->nLocalDof() );
    BOOST_CHECK_EQUAL( Xh->nLocalDofWithoutGhost(),
                       Xh->template functionSpace<0>()->nLocalDofWithoutGhost() +
                           Xh->template functionSpace<1>()->nLocalDofWithoutGhost() +
                           Xh->template functionSpace<2>()->nLocalDofWithoutGhost() );

    BOOST_CHECK_EQUAL( U.template functionSpace<0>()->nDof(), n0 );
    BOOST_CHECK_EQUAL( U.template functionSpace<1>()->nDof(), n1 );
    BOOST_CHECK_EQUAL( U.template functionSpace<2>()->nDof(), n2 );

    BOOST_CHECK_EQUAL( U.template element<0>().start(), 0 );
    BOOST_CHECK_EQUAL( U.template element<0>().size(), n0 );
    BOOST_CHECK_EQUAL( U.template element<1>().start(), n0 );
    BOOST_CHECK_EQUAL( U.template element<1>().size(), n1 );
    BOOST_CHECK_EQUAL( U.template element<2>().start(), n0 + n1 );
    BOOST_CHECK_EQUAL( U.template element<2>().size(), n2 );
}

BOOST_AUTO_TEST_CASE( legacy_helper_factories_remain_composite )
{
    auto mesh = Feel::Test::makeMesh();

    auto TH = THch<1>( mesh );
    using th_space_type = typename std::remove_reference_t<decltype( *TH )>;
    static_assert( CompositeSpaceConcept<th_space_type> );
    static_assert( th_space_type::nSpaces == 2 );
    BOOST_CHECK_EQUAL( TH->nDof(),
                       TH->template functionSpace<0>()->nDof() +
                           TH->template functionSpace<1>()->nDof() );

    auto P2 = P2ch<Lagrange<1, Vectorial>, Lagrange<0, Scalar, Discontinuous>>( mesh );
    using p2_space_type = typename std::remove_reference_t<decltype( *P2 )>;
    static_assert( CompositeSpaceConcept<p2_space_type> );
    static_assert( p2_space_type::nSpaces == 2 );
    BOOST_CHECK_EQUAL( P2->nDof(),
                       P2->template functionSpace<0>()->nDof() +
                           P2->template functionSpace<1>()->nDof() );

    auto Dh = DhPdh<0>( mesh );
    using dh_space_type = typename std::remove_reference_t<decltype( *Dh )>;
    static_assert( CompositeSpaceConcept<dh_space_type> );
    static_assert( dh_space_type::nSpaces == 2 );
    BOOST_CHECK_EQUAL( Dh->nDof(),
                       Dh->template functionSpace<0>()->nDof() +
                           Dh->template functionSpace<1>()->nDof() );
}

BOOST_AUTO_TEST_CASE( explicit_product_layer_matches_taylor_hood_dof_accounting )
{
    using namespace boost::hana::literals;
    using mesh_type = Feel::Test::mesh_type;

    auto mesh = Feel::Test::makeMesh();
    auto legacy = THch<1>( mesh );
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto ps = product( Vh, Qh );

    using product_type = decltype( ps );
    using runtime_product_type = ProductSpace<Pch_ptrtype<mesh_type, 1>, true>;

    static_assert( ProductSpaceConcept<product_type> );
    static_assert( ProductSpacesConcept<product_type> );
    static_assert( !CompositeSpaceConcept<product_type> );
    static_assert( !FunctionSpaceConcept<product_type> );
    static_assert( ProductSpaceConcept<runtime_product_type> );
    static_assert( !ProductSpacesConcept<runtime_product_type> );

    BOOST_CHECK_EQUAL( ps.numberOfSpaces(), 2 );
    BOOST_CHECK_EQUAL( ps[0_c]->nDof(), legacy->template functionSpace<0>()->nDof() );
    BOOST_CHECK_EQUAL( ps[1_c]->nDof(), legacy->template functionSpace<1>()->nDof() );
    BOOST_CHECK_EQUAL( ps.nDof(), legacy->nDof() );
    BOOST_CHECK_EQUAL( ps.nLocalDof(), legacy->nLocalDof() );
}

BOOST_AUTO_TEST_CASE( product_functionspaces_is_product_backed_composite_facade )
{
    using mesh_type = Feel::Test::mesh_type;

    auto mesh = Feel::Test::makeMesh();
    auto Vh = Pchv<2>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = productFunctionSpaces( Vh, Qh );

    using product_backed_type = std::remove_reference_t<decltype( Xh )>;

    static_assert( FunctionSpaceConcept<product_backed_type> );
    static_assert( CompositeSpaceConcept<product_backed_type> );
    static_assert( ProductBackedCompositeSpaceConcept<product_backed_type> );
    static_assert( !LegacyCompositeFunctionSpaceConcept<product_backed_type> );
    static_assert( !product_backed_type::uses_internal_composite );
    static_assert( product_backed_type::is_product_backed_composite );
    static_assert( std::is_same_v<typename product_backed_type::template sub_functionspace_type<0>,
                                  std::remove_reference_t<decltype( *Vh )>> );
    static_assert( std::is_same_v<typename product_backed_type::template sub_functionspace_type<1>,
                                  std::remove_reference_t<decltype( *Qh )>> );

    BOOST_CHECK_EQUAL( Xh.nSubFunctionSpace(), 2 );
    BOOST_CHECK_EQUAL( Xh.template functionSpace<0>(), Vh );
    BOOST_CHECK_EQUAL( Xh.template functionSpace<1>(), Qh );
    BOOST_CHECK_EQUAL( Xh.nDof(), Vh->nDof() + Qh->nDof() );
    BOOST_CHECK_EQUAL( Xh.blockMapPtr( 0 ), Vh->mapPtr() );
    BOOST_CHECK_EQUAL( Xh.blockMapPtr( 1 ), Qh->mapPtr() );
}

BOOST_AUTO_TEST_CASE( mortar_has_explicit_functionspace_abstraction )
{
    using mortar_mesh_type = Mesh<Simplex<1, 1, 2>>;
    using mortar_space_type = Moch_type<mortar_mesh_type, 1>;
    using raw_mortar_space_type =
        FunctionSpace<mortar_mesh_type,
                      bases<Lagrange<1, Scalar, Continuous, PointSetEquiSpaced>>,
                      double,
                      mortars<Mortar>>;

    static_assert( FunctionSpaceConcept<mortar_space_type> );
    static_assert( MortarFunctionSpaceConcept<mortar_space_type> );
    static_assert( NonCompositeSpaceConcept<mortar_space_type> );
    static_assert( mortar_space_type::is_mortar );
    static_assert( mortar_space_type::is_explicit_mortar_space );
    static_assert( !mortar_space_type::is_composite );
    static_assert( !mortar_space_type::has_type_level_periodicity );
    static_assert( std::is_same_v<typename mortar_space_type::periodicity_0_type, NoPeriodicity> );
    static_assert( raw_mortar_space_type::is_mortar );
    static_assert( !MortarFunctionSpaceConcept<raw_mortar_space_type> );
}

BOOST_AUTO_TEST_SUITE_END()
