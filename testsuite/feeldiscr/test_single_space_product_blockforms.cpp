/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#define BOOST_TEST_MODULE single_space_product_blockforms
#include <feel/feelcore/testsuite.hpp>

#include <memory>
#include <type_traits>

#include <feel/feeldiscr/concepts.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feeldiscr/productfunctionspaces.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

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

template<typename ProductT>
void
checkOneBlockGraphAndVector( ProductT&& ps )
{
    auto&& productSpace = remove_shared_ptr_f( std::forward<ProductT>( ps ) );

    auto graph = csrGraphBlocks( std::forward<ProductT>( ps ) );
    BOOST_CHECK_EQUAL( graph.nRow(), 1 );
    BOOST_CHECK_EQUAL( graph.nCol(), 1 );
    BOOST_CHECK( graph( 0, 0 ) );

    auto blocks = blockVector( productSpace );
    BOOST_CHECK_EQUAL( blocks.nRow(), 1 );
    BOOST_CHECK_EQUAL( blocks.nCol(), 1 );
    BOOST_CHECK( blocks( 0, 0 ) );
}
}

BOOST_AUTO_TEST_SUITE( single_space_product_blockforms )

BOOST_AUTO_TEST_CASE( one_block_product_contract )
{
    using mesh_type = Feel::Test::mesh_type;

    auto mesh = Feel::Test::makeMesh();
    auto Xh = Pch<1>( mesh );
    auto ps = product( Xh );

    using scalar_space_type = std::remove_reference_t<decltype( *Xh )>;
    using product_type = decltype( ps );

    static_assert( FunctionSpaceConcept<scalar_space_type> );
    static_assert( ProductSpaceConcept<product_type> );
    static_assert( ProductSpacesConcept<product_type> );
    static_assert( !FunctionSpaceConcept<product_type> );
    static_assert( !CompositeSpaceConcept<product_type> );

    BOOST_CHECK_EQUAL( ps.numberOfSpaces(), 1 );
    BOOST_CHECK_EQUAL( ps.template space<0>(), Xh );
    BOOST_CHECK_EQUAL( ps.space( 0_c ), Xh );
    BOOST_CHECK_EQUAL( ps.nDof(), Xh->nDof() );
    BOOST_CHECK_EQUAL( ps.nLocalDof(), Xh->nLocalDof() );
    BOOST_CHECK_EQUAL( ps.nDofStart( 0 ), 0 );
    BOOST_CHECK_EQUAL( ps.nDofStart( 1 ), ps.nDof() );
    BOOST_CHECK_EQUAL( ps.nLocalDofStart( 1 ), ps.nLocalDof() );
    BOOST_CHECK_EQUAL( ps.blockDofStart( 1 ), ps.nDofStart( 1 ) );
    BOOST_CHECK_EQUAL( ps.blockLocalDofStart( 1 ), ps.nLocalDofStart( 1 ) );
    BOOST_CHECK_EQUAL( ps.blockMapPtr( 0 ), Xh->mapPtr() );

    auto U = ps.element();
    static_assert( BlockElementConcept<decltype( U )> );
    BOOST_CHECK_EQUAL( U.functionSpace().template space<0>(), Xh );
    BOOST_CHECK_EQUAL( U( 0_c ).functionSpace(), Xh );

    Feel::Test::checkOneBlockGraphAndVector( ps );
    auto const& psConstRef = ps;
    Feel::Test::checkOneBlockGraphAndVector( psConstRef );
    auto psPtr = std::make_shared<product_type>( ps );
    Feel::Test::checkOneBlockGraphAndVector( psPtr );
    std::shared_ptr<product_type const> psConstPtr = psPtr;
    Feel::Test::checkOneBlockGraphAndVector( psConstPtr );
}

BOOST_AUTO_TEST_CASE( one_block_blockform_matches_ordinary_form )
{
    auto mesh = Feel::Test::makeMesh();
    auto Xh = Pch<1>( mesh );
    auto ps = product( Xh );
    auto b = backend( _rebuild = true );

    auto uOrd = Xh->element( "u_ordinary" );
    auto vOrd = Xh->element( "v_ordinary" );
    auto aOrd = form2( _trial = Xh, _test = Xh, _backend = b );
    auto lOrd = form1( _test = Xh, _backend = b );
    aOrd = integrate( _range = elements( mesh ), _expr = idt( uOrd )*id( vOrd ) );
    lOrd = integrate( _range = elements( mesh ), _expr = cst( 2.0 )*id( vOrd ) );
    aOrd.close();
    lOrd.close();
    aOrd.solve( _rhs = lOrd, _solution = uOrd, _rebuild = true );

    auto U = ps.element();
    auto& uBlock = U( 0_c );
    auto aBlock = blockform2( ps, solve::strategy::monolithic, b );
    auto lBlock = blockform1( ps, solve::strategy::monolithic, b );
    aBlock( 0_c, 0_c ) = integrate( _range = elements( mesh ), _expr = idt( uBlock )*id( uBlock ) );
    lBlock( 0_c ) = integrate( _range = elements( mesh ), _expr = cst( 2.0 )*id( uBlock ) );
    aBlock.close();
    lBlock.close();
    aBlock.solve( _rhs = lBlock, _solution = U, _rebuild = true );

    auto const diff = normL2( _range = elements( mesh ), _expr = idv( uOrd ) - idv( U( 0_c ) ) );
    auto const exact = normL2( _range = elements( mesh ), _expr = idv( U( 0_c ) ) - cst( 2.0 ) );
    BOOST_CHECK_SMALL( diff, 1e-8 );
    BOOST_CHECK_SMALL( exact, 1e-8 );
}

BOOST_AUTO_TEST_CASE( shared_pointer_inputs_assemble_like_value_inputs )
{
    auto mesh = Feel::Test::makeMesh();
    auto Xh = Pch<1>( mesh );
    auto ps = product( Xh );
    using product_type = decltype( ps );

    auto psPtr = std::make_shared<product_type>( ps );
    std::shared_ptr<product_type const> psConstPtr = psPtr;
    auto b = backend( _rebuild = true );

    auto U = ps.element();
    auto& u = U( 0_c );
    auto a = blockform2( psPtr, solve::strategy::monolithic, b );
    auto l = blockform1( psPtr, solve::strategy::monolithic, b );
    a( 0_c, 0_c ) += integrate( _range = elements( mesh ), _expr = idt( u )*id( u ) );
    l( 0_c ) += integrate( _range = elements( mesh ), _expr = cst( 3.0 )*id( u ) );
    a.close();
    l.close();
    a.solve( _rhs = l, _solution = U, _rebuild = true );

    auto const err = normL2( _range = elements( mesh ), _expr = idv( U( 0_c ) ) - cst( 3.0 ) );
    BOOST_CHECK_SMALL( err, 1e-8 );

    auto aConst = blockform2( psConstPtr, solve::strategy::monolithic, b );
    auto lConst = blockform1( psConstPtr, solve::strategy::monolithic, b );
    BOOST_CHECK( aConst.matrixPtr() );
    BOOST_CHECK( lConst.vectorPtr() );
}

BOOST_AUTO_TEST_SUITE_END()
