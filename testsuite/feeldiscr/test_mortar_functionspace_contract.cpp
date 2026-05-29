/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#define BOOST_TEST_MODULE mortar functionspace contract
#include <feel/feelcore/testsuite.hpp>

#include <type_traits>

#include <feel/feeldiscr/concepts.hpp>
#include <feel/feeldiscr/moch.hpp>
#include <feel/feeldiscr/mortarfunctionspace.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/geotool.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( mortar_functionspace_contract )

BOOST_AUTO_TEST_CASE( explicit_mortar_type_contract )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<1, 1, 2>>;
    using mortar_space_type = MortarLagrangeSpace<mesh_type, 1>;
    using moch_space_type = Moch_type<mesh_type, 1>;
    using ordinary_space_type = Pch_type<mesh_type, 1, double, PointSetEquiSpaced>;
    using raw_mortar_space_type =
        FunctionSpace<mesh_type,
                      bases<Lagrange<1, Scalar, Continuous, PointSetEquiSpaced>>,
                      double,
                      mortars<Mortar>>;

    static_assert( FunctionSpaceConcept<mortar_space_type> );
    static_assert( MortarFunctionSpaceConcept<mortar_space_type> );
    static_assert( NonCompositeSpaceConcept<mortar_space_type> );
    static_assert( !CompositeSpaceConcept<mortar_space_type> );
    static_assert( mortar_space_type::is_mortar );
    static_assert( mortar_space_type::is_explicit_mortar_space );
    static_assert( std::is_same_v<typename mortar_space_type::mortar_policy_type, Mortar> );
    static_assert( std::is_same_v<mortar_space_type, moch_space_type> );

    static_assert( FunctionSpaceConcept<ordinary_space_type> );
    static_assert( !ordinary_space_type::is_mortar );
    static_assert( !MortarFunctionSpaceConcept<ordinary_space_type> );

    static_assert( raw_mortar_space_type::is_mortar );
    static_assert( !MortarFunctionSpaceConcept<raw_mortar_space_type> );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( explicit_mortar_factory_matches_moch )
{
    using namespace Feel;
    using mesh_type = Mesh<Simplex<1, 1, 2>>;

    double meshSize = doption( _name = "gmsh.hsize" );
    GeoTool::Node x0( 0, 0 );
    GeoTool::Node x1( 1, 0 );
    GeoTool::Line line( meshSize, "MORTAR_CONTRACT_LINE", x0, x1 );
    auto mesh = line.createMesh( _mesh = new mesh_type, _name = "mortarContractLine" );

    auto Mh = mortarFunctionSpace<1>( mesh );
    auto MochSpace = Moch<1>( mesh );

    BOOST_CHECK( Mh->isMortar() );
    BOOST_CHECK( MochSpace->isMortar() );
    BOOST_CHECK( Mh->dof() );
    BOOST_CHECK_GT( Mh->nDof(), size_type( 0 ) );
    BOOST_CHECK_EQUAL( Mh->nDof(), Mh->dof()->nDof() );
    BOOST_CHECK_EQUAL( MochSpace->nDof(), Mh->nDof() );

    auto lambda = Mh->element( "lambda" );
    BOOST_CHECK_EQUAL( lambda.nDof(), Mh->nDof() );
}

BOOST_AUTO_TEST_SUITE_END()
