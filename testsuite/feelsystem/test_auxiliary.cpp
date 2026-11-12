/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE feelsystem_auxiliary
#include <feel/feelcore/testsuite.hpp>

#include <boost/hana/integral_constant.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

#if defined( FEELPP_HAS_PETSC )
#include <feel/feelalg/petschpddm.hpp>
#include <feel/feelalg/preconditionerpetsc.hpp>
#include <feel/feelalg/topetsc.hpp>
#if defined( FEELPP_HAS_SLEPC )
#include <slepceps.h>
#endif
#endif

namespace
{

Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "FeelSystem auxiliary operator test options" );
    options.add( Feel::backend_options( "auxiliary-scalar" ) )
        .add( Feel::backend_options( "auxiliary-block" ) )
        .add( Feel::backend_options( "auxiliary-geneo" ) );
    return options;
}

template <typename AuxiliaryOperator>
double constantEnergy( AuxiliaryOperator const& op )
{
    auto ones = op.backend()->newVector( op.matrixPtr()->mapRowPtr() );
    ones->setConstant( 1.0 );
    ones->close();
    return op.matrixPtr()->energy( *ones, *ones );
}

#if defined( FEELPP_HAS_SLEPC ) && defined( FEELPP_HAS_PETSC_H )
template <typename AuxiliaryOperator>
std::vector<double> smallestEigenvalues( AuxiliaryOperator const& op, int nev )
{
    using namespace Feel;

    auto matrix = toPETSc( op.matrixPtr() );
    BOOST_REQUIRE( matrix );

    EPS eps = nullptr;
    int ierr = EPSCreate( PETSC_COMM_SELF, &eps );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    ierr = EPSSetOperators( eps, matrix->mat(), nullptr );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    ierr = EPSSetProblemType( eps, EPS_HEP );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    ierr = EPSSetType( eps, EPSLAPACK );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    ierr = EPSSetWhichEigenpairs( eps, EPS_SMALLEST_MAGNITUDE );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    ierr = EPSSetDimensions( eps, nev, PETSC_DEFAULT, PETSC_DEFAULT );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    ierr = EPSSolve( eps );
    BOOST_REQUIRE_EQUAL( ierr, 0 );

    PetscInt converged = 0;
    ierr = EPSGetConverged( eps, &converged );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    BOOST_REQUIRE_GE( converged, nev );

    std::vector<double> values;
    for ( PetscInt i = 0; i < nev; ++i )
    {
        PetscScalar real = 0, imaginary = 0;
        ierr = EPSGetEigenvalue( eps, i, &real, &imaginary );
        BOOST_REQUIRE_EQUAL( ierr, 0 );
        BOOST_CHECK_SMALL( PetscAbsScalar( imaginary ), 1e-12 );
        values.push_back( PetscRealPart( real ) );
    }
    ierr = EPSDestroy( &eps );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    std::sort( values.begin(), values.end(),
               []( double lhs, double rhs ) { return std::abs( lhs ) < std::abs( rhs ); } );
    return values;
}
#endif

} // namespace

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system_auxiliary" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( feelsystem_auxiliary_suite )

BOOST_AUTO_TEST_CASE( scalar_local_neumann_is_sequential_floating_and_mapped )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto stiffness = integrate( _range = elements( mesh ),
                                _expr = inner( grad( u ), grad( v ) ) );
    auto localStiffness = integrate( _range = elements( mesh, entity_process_t::ALL ),
                                     _expr = inner( grad( u ), grad( v ) ) );
    auto mass = integrate( _range = elements( mesh ),
                           _expr = u * v );

    auto problem =
        linearProblem( _space = Xh, _name = "auxiliary-scalar" )
            .lhs( stiffness )
            .lhs( mass )
            .rhs( integrate( _range = elements( mesh ), _expr = v ) )
            .dirichlet( _range = boundaryfaces( mesh ),
                        _element = uh,
                        _expr = cst( 0.0 ) )
            .localNeumann( "geneo", localStiffness )
            .localNeumannFromLhs( "lhs-copy" );

    BOOST_REQUIRE_EQUAL( problem.auxiliaryOperatorDescriptors().size(), 2 );
    auto system = problem.assemble();
    BOOST_REQUIRE( system.hasAuxiliaryOperator( "geneo" ) );
    BOOST_REQUIRE( system.hasAuxiliaryOperator( "lhs-copy" ) );

    auto const& geneo = system.auxiliaryOperator( "geneo" );
    BOOST_REQUIRE( geneo.matrixPtr() );
    BOOST_CHECK_EQUAL( geneo.matrixPtr()->worldComm().globalSize(), 1 );
    BOOST_CHECK_EQUAL( geneo.localSize(), Xh->dofOn()->nLocalDofWithGhost() );
    BOOST_CHECK_GT( geneo.matrixPtr()->l1Norm(), 0.0 );
    BOOST_CHECK_SMALL( constantEnergy( geneo ), 1e-9 );

#if defined( FEELPP_HAS_SLEPC ) && defined( FEELPP_HAS_PETSC_H )
    auto const eigenvalues = smallestEigenvalues( geneo, 2 );
    BOOST_REQUIRE_GE( eigenvalues.size(), 2 );
    BOOST_TEST_MESSAGE( "local Neumann eigenvalues: lambda0=" << eigenvalues[0]
                        << ", lambda1=" << eigenvalues[1] );
    BOOST_CHECK_SMALL( std::abs( eigenvalues[0] ), 1e-10 );
    BOOST_CHECK_GT( eigenvalues[1], 1e-8 );
#endif

    auto const& sourceIds = Xh->dofOn()->mapGlobalProcessToGlobalCluster();
    BOOST_REQUIRE_EQUAL( geneo.globalDofIds().size(), sourceIds.size() );
    BOOST_CHECK_EQUAL_COLLECTIONS( geneo.globalDofIds().begin(), geneo.globalDofIds().end(),
                                   sourceIds.begin(), sourceIds.end() );

    auto const& lhsCopy = system.auxiliaryOperator( "lhs-copy" );
    BOOST_CHECK_GT( lhsCopy.matrixPtr()->l1Norm(), 0.0 );
}

BOOST_AUTO_TEST_CASE( block_local_neumann_preserves_composite_overlap_layout )
{
    using namespace Feel;
    using namespace Feel::vf;
    using namespace boost::hana::literals;

    auto mesh = unitSquare();
    auto Vh = Pch<1>( mesh );
    auto Qh = Pch<1>( mesh );
    auto Xh = product( Vh, Qh );
    auto solution = Xh.element();
    auto u = trial( Vh );
    auto v = test( Vh );
    auto p = trial( Qh );
    auto q = test( Qh );

    auto problem =
        mixedProblem( _space = Xh, _name = "auxiliary-block" )
            .block( 0_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = inner( grad( u ), grad( v ) ) + u * v ) )
            .block( 1_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = inner( grad( p ), grad( q ) ) + p * q ) )
            .row( 0_c, integrate( _range = elements( mesh ), _expr = v ) )
            .row( 1_c, integrate( _range = elements( mesh ), _expr = q ) )
            .localNeumann(
                "geneo",
                [mesh, u, v, p, q]( auto& a )
                {
                    a( 0_c, 0_c ) += integrate( _range = elements( mesh, entity_process_t::ALL ),
                                                _expr = inner( grad( u ), grad( v ) ) );
                    a( 1_c, 1_c ) += integrate( _range = elements( mesh, entity_process_t::ALL ),
                                                _expr = inner( grad( p ), grad( q ) ) );
                } );

    auto system = problem.assemble();
    auto const& geneo = system.auxiliaryOperator( "geneo" );
    BOOST_REQUIRE( geneo.matrixPtr() );
    BOOST_CHECK_EQUAL( geneo.matrixPtr()->worldComm().globalSize(), 1 );
    BOOST_CHECK_EQUAL( geneo.localSize(), system.baseMatrixPtr()->mapRowPtr()->nLocalDofWithGhost() );
    BOOST_CHECK_GT( geneo.matrixPtr()->l1Norm(), 0.0 );
    BOOST_CHECK_SMALL( constantEnergy( geneo ), 1e-9 );

    auto const distributedMap = system.baseMatrixPtr()->mapRowPtr();
    bool idsInRange = true;
    for ( auto id : geneo.globalDofIds() )
        idsInRange = idsInRange && id < distributedMap->nDof();
    BOOST_CHECK( idsInRange );

    auto const result = system.solve( _solution = solution );
    BOOST_REQUIRE_MESSAGE( result.converged, result.message );
}

#if defined( FEELPP_HAS_PETSC )
BOOST_AUTO_TEST_CASE( scalar_problem_api_drives_hpddm_geneo_solve )
{
    using namespace Feel;
    using namespace Feel::vf;

    if ( !petscHasHpddmSymbol( "PCCreate_HPDDM" ) ||
         !petscHasHpddmSymbol( "PCHPDDMSetAuxiliaryMat" ) )
    {
        BOOST_TEST_MESSAGE( "Skipping FeelSystem GenEO test: PETSc HPDDM auxiliary API is unavailable" );
        return;
    }

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto solution = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto problem =
        linearProblem( _space = Xh, _name = "auxiliary-geneo" )
            .diffusion( _range = elements( mesh ),
                        _range_extended = elements( mesh, entity_process_t::ALL ),
                        _expr = inner( grad( u ), grad( v ) ) )
            .robin( _range = boundaryfaces( mesh ),
                    _range_extended = boundaryfaces( mesh, entity_process_t::ALL ),
                    _trial = u, _test = v,
                    _alpha = cst( 1.0 ), _expr = cst( 0.0 ) )
            .dirichlet( _range = boundaryfaces( mesh ),
                        _range_extended = boundaryfaces( mesh, entity_process_t::ALL ),
                        _element = solution, _expr = cst( 0.0 ) )
            .rhs( integrate( _range = elements( mesh ), _expr = v ) )
            ;

    auto system = problem.assemble();
    auto const result = system.solve( _solution = solution );
    BOOST_REQUIRE_MESSAGE( result.converged, result.message );
    BOOST_REQUIRE( system.hasAuxiliaryOperator( "geneo" ) );
    BOOST_CHECK_GT( constantEnergy( system.auxiliaryOperator( "geneo" ) ), 1e-10 );
#if defined( FEELPP_HAS_SLEPC ) && defined( FEELPP_HAS_PETSC_H )
    auto const inheritedEigenvalues = smallestEigenvalues(
        system.auxiliaryOperator( "geneo" ), 2 );
    BOOST_TEST_MESSAGE( "inherited Robin+Dirichlet eigenvalues: lambda0="
                        << inheritedEigenvalues[0]
                        << ", lambda1=" << inheritedEigenvalues[1] );
    BOOST_CHECK_GT( inheritedEigenvalues[1], 1e-8 );
    if ( Environment::worldComm().globalSize() == 1 )
        BOOST_CHECK_GT( inheritedEigenvalues[0], 1e-8 );
#endif

    auto petscPreconditioner = toPETSc( system.backend()->preconditioner() );
    BOOST_REQUIRE( petscPreconditioner );
    BOOST_REQUIRE( petscPreconditioner->hasHpddmAuxiliaryMatrix() );
    BOOST_REQUIRE( petscPreconditioner->hasHpddmAuxiliaryIS() );
    BOOST_CHECK_EQUAL( petscPreconditioner->hpddmAuxiliaryMatrix().get(),
                       system.auxiliaryOperator( "geneo" ).matrixPtr().get() );

    auto const& overlap = system.auxiliaryOperator( "geneo" ).globalDofIds();
    IS is = reinterpret_cast<IS>( petscPreconditioner->hpddmAuxiliaryIS().get() );
    PetscInt localSize = 0;
    int ierr = ISGetLocalSize( is, &localSize );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    BOOST_REQUIRE_EQUAL( static_cast<std::size_t>( localSize ), overlap.size() );
    PetscInt const* indices = nullptr;
    ierr = ISGetIndices( is, &indices );
    BOOST_REQUIRE_EQUAL( ierr, 0 );
    for ( PetscInt i = 0; i < localSize; ++i )
        BOOST_CHECK_EQUAL( indices[i], static_cast<PetscInt>( overlap[i] ) );
    ierr = ISRestoreIndices( is, &indices );
    BOOST_REQUIRE_EQUAL( ierr, 0 );

    if ( Environment::worldComm().globalSize() > 1 )
    {
        using get_complexities_type = PetscErrorCode (*)( PC, PetscReal*, PetscReal* );
        auto getComplexities = petscHpddmSymbol<get_complexities_type>(
            "PCHPDDMGetComplexities" );
        BOOST_REQUIRE( getComplexities );
        PetscReal gridComplexity = 0, operatorComplexity = 0;
        ierr = getComplexities( petscPreconditioner->pc(),
                                &gridComplexity, &operatorComplexity );
        BOOST_REQUIRE_EQUAL( ierr, 0 );
        BOOST_CHECK_GT( gridComplexity, 1.0 );
        BOOST_CHECK_GT( operatorComplexity, 1.0 );
    }
}
#endif

BOOST_AUTO_TEST_SUITE_END()
