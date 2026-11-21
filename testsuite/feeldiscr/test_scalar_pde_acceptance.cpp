/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file test_scalar_pde_acceptance.cpp
 * @brief Static/runtime scalar PDE, submesh, and restart acceptance checks.
 */
#define BOOST_TEST_MODULE test_scalar_pde_acceptance
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelpoly/fecapabilities.hpp>
#include <feel/feelvf/vf.hpp>

#include "fe_parity_checks.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <fmt/format.h>

#include <filesystem>
#include <sstream>
#include <utility>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template<typename MeshPtrType>
void
checkPoissonParity( MeshPtrType const& mesh,
                    bool checkpoint = false,
                    bool exportSolution = false )
{
    fe::requireCapability<fe::Family::Lagrange>(
        *mesh, RuntimeOrder{ 2 }, fe::Operation::LinearSolve );

    auto staticSpace = Pch<2>( mesh );
    auto runtimeSpace = Pch<Dynamic>( mesh, RuntimeOrder{ 2 } );
    auto manufactured =
        cst( 0.75 ) + Px() - 0.5*Py() + 0.2*Px()*Py();

    auto staticTarget = staticSpace->element();
    auto runtimeTarget = runtimeSpace->element();
    staticTarget.on( _range = elements( mesh ), _expr = manufactured );
    runtimeTarget.on( _range = elements( mesh ), _expr = manufactured );

    auto assembleSystem = [&]( auto const& space, auto const& target )
    {
        auto trial = space->element();
        auto test = space->element();
        auto op = form2( _trial = space, _test = space );
        auto rhs = form1( _test = space );
        op = integrate(
            _range = elements( mesh ),
            _expr = inner( gradt( trial ), grad( test ) )+
                    idt( trial )*id( test ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        rhs = integrate(
            _range = elements( mesh ),
            _expr = inner( gradv( target ), grad( test ) )+
                    idv( target )*id( test ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        op += on( _range = boundaryfaces( mesh ),
                  _rhs = rhs,
                  _element = trial,
                  _expr = manufactured );
        return std::make_pair( std::move( op ), std::move( rhs ) );
    };
    auto [staticOperator, staticRhs] =
        assembleSystem( staticSpace, staticTarget );
    auto [runtimeOperator, runtimeRhs] =
        assembleSystem( runtimeSpace, runtimeTarget );

    Feel::test_detail::checkAssembledMatrixParity(
        staticOperator, runtimeOperator, 5e-11 );
    Feel::test_detail::checkAssembledVectorParity(
        staticRhs, runtimeRhs, 5e-11 );

    auto staticSolution = staticSpace->element();
    auto runtimeSolution = runtimeSpace->element();
    auto const staticSolve = staticOperator.solve(
        _rhs = staticRhs, _solution = staticSolution,
        _rebuild = true );
    auto const runtimeSolve = runtimeOperator.solve(
        _rhs = runtimeRhs, _solution = runtimeSolution,
        _rebuild = true );
    Feel::test_detail::checkSolveParity( staticSolve, runtimeSolve );

    BOOST_CHECK_SMALL(
        normL2( _range = elements( mesh ),
                _expr = idv( staticSolution )-idv( runtimeSolution ),
                _quad = _Q<14>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-9 );
    BOOST_CHECK_SMALL(
        normL2( _range = boundaryfaces( mesh ),
                _expr = idv( staticSolution )-idv( runtimeSolution ),
                _quad = _Q<14>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-9 );

    if ( exportSolution && Environment::worldComm().globalSize() == 1 )
    {
        fe::requireCapability<fe::Family::Lagrange>(
            *mesh, RuntimeOrder{ 2 }, fe::Operation::Export );
        auto solutionExporter = exporter(
            _mesh = mesh,
            _name = fmt::format(
                "scalar_pde_runtime_export_g{}", mesh->order() ) );
        BOOST_REQUIRE( solutionExporter );
        BOOST_REQUIRE( solutionExporter->doExport() );
        solutionExporter->add(
            "runtime_solution", runtimeSolution );
        solutionExporter->save();
        BOOST_REQUIRE( std::filesystem::exists(
            solutionExporter->path() ) );
        BOOST_CHECK( !std::filesystem::is_empty(
            solutionExporter->path() ) );
    }

    if ( checkpoint && Environment::worldComm().globalSize() == 1 )
    {
        fe::requireCapability<fe::Family::Lagrange>(
            *mesh, RuntimeOrder{ 2 }, fe::Operation::CheckpointRestart );
        auto const solutionNorm = normL2(
            _range = elements( mesh ), _expr = idv( runtimeSolution ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        std::stringstream meshBuffer;
        std::stringstream solutionBuffer;
        {
            boost::archive::text_oarchive meshArchive{ meshBuffer };
            meshArchive << *mesh;
        }
        {
            boost::archive::text_oarchive solutionArchive{ solutionBuffer };
            solutionArchive << runtimeSolution;
        }

        using mesh_type = typename MeshPtrType::element_type;
        auto restartedMesh = mesh_type::New( RuntimeOrder{ 1 } );
        {
            boost::archive::text_iarchive meshArchive{ meshBuffer };
            meshArchive >> *restartedMesh;
        }
        BOOST_REQUIRE_EQUAL( restartedMesh->order(), mesh->order() );
        auto restartedSpace = Pch<Dynamic>(
            restartedMesh, RuntimeOrder{ 2 } );
        auto restartedSolution = restartedSpace->element();
        {
            boost::archive::text_iarchive solutionArchive{ solutionBuffer };
            solutionArchive >> restartedSolution;
        }
        BOOST_CHECK_SMALL(
            normL2( _range = elements( restartedMesh ),
                    _expr = idv( restartedSolution ),
                    _quad = _Q<14>(),
                    _geomap = GeomapStrategyType::GEOMAP_HO )-
                solutionNorm,
            2e-9 );
    }
}
} // namespace

BOOST_AUTO_TEST_SUITE( scalar_pde_acceptance )

BOOST_AUTO_TEST_CASE( poisson_affine_quadratic_and_cubic_geometry )
{
    for ( uint16_type geometryOrder :
          { uint16_type{ 1 }, uint16_type{ 2 }, uint16_type{ 3 } } )
    {
        auto mesh = unitCircle(
            RuntimeOrder{ geometryOrder }, 0.8 );
        BOOST_REQUIRE( mesh );
        BOOST_REQUIRE_GT( mesh->numGlobalElements(), 0 );
        checkPoissonParity(
            mesh, geometryOrder == 2, geometryOrder <= 2 );
    }
}

BOOST_AUTO_TEST_CASE( poisson_on_runtime_order_submesh )
{
    auto mesh = unitCircle( RuntimeOrder{ 2 }, 0.8 );
    auto submesh = createSubmesh(
        _mesh = mesh, _range = elements( mesh ) );
    BOOST_REQUIRE( submesh );
    BOOST_CHECK_EQUAL( submesh->order(), mesh->order() );
    fe::requireCapability<fe::Family::Lagrange>(
        *submesh, RuntimeOrder{ 2 }, fe::Operation::Submesh );
    checkPoissonParity( submesh );
}

BOOST_AUTO_TEST_SUITE_END()
