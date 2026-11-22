/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#define BOOST_TEST_MODULE productspace contract
#include <feel/feelcore/testsuite.hpp>

#include <type_traits>

#include <feel/feeldiscr/concepts.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/product.hpp>
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

BOOST_AUTO_TEST_SUITE( productspace_contract )

BOOST_AUTO_TEST_CASE( compile_time_product_contract )
{
    using mesh_type = Feel::Test::mesh_type;

    auto mesh = Feel::Test::makeMesh();
    auto Xh = Pch<1>( mesh );
    auto Vh = Pchv<2>( mesh );
    auto ps = product( Xh, Vh );

    using scalar_space_type = std::remove_reference_t<decltype( *Xh )>;
    using vectorial_space_type = std::remove_reference_t<decltype( *Vh )>;
    using product_type = decltype( ps );

    static_assert( FunctionSpaceConcept<scalar_space_type> );
    static_assert( NonCompositeSpaceConcept<vectorial_space_type> );
    static_assert( !ProductSpaceConcept<scalar_space_type> );
    static_assert( ProductSpaceConcept<product_type> );
    static_assert( ProductSpacesConcept<product_type> );

    auto U = ps.element();
    static_assert( BlockElementConcept<decltype( U )> );

    BOOST_CHECK_EQUAL( ps.numberOfSpaces(), 2 );
    BOOST_CHECK_EQUAL( ps.template space<0>(), Xh );
    BOOST_CHECK_EQUAL( ps.template space<1>(), Vh );
    BOOST_CHECK_EQUAL( ps.space( 0_c ), Xh );
    BOOST_CHECK_EQUAL( ps.space( 1_c ), Vh );

    BOOST_CHECK_EQUAL( ps.nDof(), Xh->nDof() + Vh->nDof() );
    BOOST_CHECK_EQUAL( ps.nLocalDof(), Xh->nLocalDof() + Vh->nLocalDof() );
    BOOST_CHECK_EQUAL( ps.nDofStart( 0 ), 0 );
    BOOST_CHECK_EQUAL( ps.nDofStart( 1 ), Xh->nDof() );
    BOOST_CHECK_EQUAL( ps.nDofStart( 2 ), ps.nDof() );
    BOOST_CHECK_EQUAL( ps.blockDofStart( 1 ), ps.nDofStart( 1 ) );
    BOOST_CHECK_EQUAL( ps.nLocalDofStart( 1 ), Xh->nLocalDof() );
    BOOST_CHECK_EQUAL( ps.blockLocalDofStart( 1 ), ps.nLocalDofStart( 1 ) );
    BOOST_CHECK_EQUAL( ps.blockMapPtr( 0 ), Xh->mapPtr() );
    BOOST_CHECK_EQUAL( ps.blockMapPtr( 1 ), Vh->mapPtr() );
}

BOOST_AUTO_TEST_CASE( runtime_same_type_product_contract )
{
    using mesh_type = Feel::Test::mesh_type;
    using scalar_ptrtype = Pch_ptrtype<mesh_type, 1>;
    using runtime_product_type = ProductSpace<scalar_ptrtype, true>;

    auto mesh = Feel::Test::makeMesh();
    auto Xh = Pch<1>( mesh );
    runtime_product_type ps( 3, Xh );

    static_assert( ProductSpaceConcept<runtime_product_type> );
    static_assert( !ProductSpacesConcept<runtime_product_type> );

    auto U = ps.element();
    static_assert( BlockElementConcept<decltype( U )> );

    BOOST_CHECK_EQUAL( ps.numberOfSpaces(), 3 );
    BOOST_CHECK_EQUAL( ps.space( 0 ), Xh );
    BOOST_CHECK_EQUAL( ps.space( 1 ), Xh );
    BOOST_CHECK_EQUAL( ps.template space<2>(), Xh );
    BOOST_CHECK_EQUAL( ps.nDof(), 3*Xh->nDof() );
    BOOST_CHECK_EQUAL( ps.nLocalDof(), 3*Xh->nLocalDof() );
    BOOST_CHECK_EQUAL( ps.nDofStart( 0 ), 0 );
    BOOST_CHECK_EQUAL( ps.nDofStart( 2 ), 2*Xh->nDof() );
    BOOST_CHECK_EQUAL( ps.nDofStart( 3 ), ps.nDof() );
    BOOST_CHECK_EQUAL( ps.nLocalDofStart( 2 ), 2*Xh->nLocalDof() );
    BOOST_CHECK_EQUAL( ps.blockMapPtr( 2 ), Xh->mapPtr() );
}

BOOST_AUTO_TEST_CASE( mixed_static_runtime_product_contract )
{
    using mesh_type = Feel::Test::mesh_type;
    using scalar_ptrtype = Pch_ptrtype<mesh_type, 1>;
    using runtime_product_type = ProductSpace<scalar_ptrtype, true>;

    auto mesh = Feel::Test::makeMesh();
    auto Xh = Pch<1>( mesh );
    auto Vh = Pchv<2>( mesh );
    auto repeated = std::make_shared<runtime_product_type>( 3, Xh );
    auto ps = product2( repeated, Xh, Vh );

    using mixed_product_type = decltype( ps );
    static_assert( ProductSpaceConcept<mixed_product_type> );
    static_assert( ProductSpacesConcept<mixed_product_type> );

    auto U = ps.element();
    static_assert( BlockElementConcept<decltype( U )> );

    auto const staticDof = Xh->nDof() + Vh->nDof();
    auto const staticLocalDof = Xh->nLocalDof() + Vh->nLocalDof();

    BOOST_CHECK_EQUAL( ps.numberOfSpaces(), 5 );
    BOOST_CHECK_EQUAL( ps.template space<0>(), Xh );
    BOOST_CHECK_EQUAL( ps.template space<1>(), Vh );
    BOOST_CHECK_EQUAL( ps.template space<2>(), repeated );
    BOOST_CHECK_EQUAL( ps.nDof(), staticDof + repeated->nDof() );
    BOOST_CHECK_EQUAL( ps.nLocalDof(), staticLocalDof + repeated->nLocalDof() );
    BOOST_CHECK_EQUAL( ps.nDofStart( 0 ), 0 );
    BOOST_CHECK_EQUAL( ps.nDofStart( 1 ), Xh->nDof() );
    BOOST_CHECK_EQUAL( ps.nDofStart( 2 ), staticDof );
    BOOST_CHECK_EQUAL( ps.nDofStart( 3 ), staticDof + Xh->nDof() );
    BOOST_CHECK_EQUAL( ps.nDofStart( ps.numberOfSpaces() ), ps.nDof() );
    BOOST_CHECK_EQUAL( ps.nLocalDofStart( 2 ), staticLocalDof );
    BOOST_CHECK_EQUAL( ps.blockDofStart( 4 ), ps.nDofStart( 4 ) );
    BOOST_CHECK_EQUAL( ps.blockLocalDofStart( 4 ), ps.nLocalDofStart( 4 ) );
    BOOST_CHECK_EQUAL( ps.blockMapPtr( 0 ), Xh->mapPtr() );
    BOOST_CHECK_EQUAL( ps.blockMapPtr( 1 ), Vh->mapPtr() );
    BOOST_CHECK_EQUAL( ps.blockMapPtr( 2 ), Xh->mapPtr() );
    BOOST_CHECK_EQUAL( ps.blockMapPtr( 4 ), Xh->mapPtr() );
}

BOOST_AUTO_TEST_SUITE_END()
