/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#define BOOST_TEST_MODULE tao testsuite

#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/program_options/parsers.hpp>

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelopt/solveroptimizationpetsc.hpp>
#include <feel/feelalg/vectorublas.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace Feel
{

namespace
{

/**
 * Build a contiguous, ghost-free distributed map for optimization vectors.
 *
 * @param localSize number of entries owned by each rank
 * @param worldComm communicator defining the rank layout
 * @return fully initialized distributed data map
 */
std::shared_ptr<DataMap<>>
makeContiguousMap( uint32_type localSize, worldcomm_ptr_t const& worldComm )
{
    auto map = std::make_shared<DataMap<>>( worldComm );
    auto const globalSize = localSize * worldComm->globalSize();
    map->setNDof( globalSize );

    for ( rank_type rank = 0; rank < worldComm->globalSize(); ++rank )
    {
        auto const first = localSize * rank;
        map->setNLocalDofWithoutGhost( rank, localSize );
        map->setNLocalDofWithGhost( rank, localSize );
        map->setFirstDofGlobalCluster( rank, first );
        map->setLastDofGlobalCluster( rank, first + localSize - 1 );
    }

    auto const first = localSize * worldComm->globalRank();
    map->resizeMapGlobalProcessToGlobalCluster( localSize );
    for ( uint32_type index = 0; index < localSize; ++index )
        map->setMapGlobalProcessToGlobalCluster( index, first + index );

    map->initNumberOfDofIdToContainerId( 1 );
    map->initDofIdToContainerIdIdentity( 0, localSize );
    map->updateWorldIndexForUse();
    map->buildIndexSplit();
    return map;
}

/**
 * Create a distributed PETSc vector with a constant initial value.
 *
 * @param localSize number of entries owned by each rank
 * @param worldComm communicator defining the rank layout
 * @param value initial value for every component
 * @return initialized distributed vector
 */
std::shared_ptr<VectorPetscMPI<double>>
makeDistributedVector( uint32_type localSize, worldcomm_ptr_t const& worldComm,
                       double value = 0.0 )
{
    auto vector = std::make_shared<VectorPetscMPI<double>>(
        makeContiguousMap( localSize, worldComm ) );
    vector->setConstant( value );
    return vector;
}

/**
 * Evaluate the distributed quadratic objective around a constant minimizer.
 *
 * @param state distributed optimization state
 * @param minimizer constant minimizer for every vector entry
 * @return one half of the squared Euclidean distance to the minimizer
 */
double
quadraticObjective( Vector<double> const& state, double minimizer )
{
    return 0.5 * ( state.dot( state ) - 2.0 * minimizer * state.sum() +
                   minimizer * minimizer * state.size() );
}

/**
 * Fill the gradient of the distributed quadratic objective.
 *
 * @param state distributed optimization state
 * @param minimizer constant minimizer for every vector entry
 * @param gradient output gradient vector
 */
void
quadraticGradient( Vector<double> const& state, double minimizer,
                   Vector<double>& gradient )
{
    gradient = state;
    gradient.add( -minimizer );
}

/**
 * Throw a descriptive exception when a PETSc test setup operation fails.
 *
 * @param error PETSc error code to inspect
 * @param operation name of the failed PETSc operation
 */
void
checkPetscSetup( PetscErrorCode error, char const* operation )
{
    if ( !error )
        return;

    char const* message = nullptr;
    PetscErrorMessage( error, &message, nullptr );
    throw std::runtime_error(
        std::string( operation ) + " failed: " + ( message ? message : "unknown PETSc error" ) );
}

/** RAII registration of one process-local PETSc runtime option for a test. */
class ScopedPetscOption
{
public:
    /**
     * Set a PETSc option until this object is destroyed.
     *
     * @param name complete option key including its leading dash
     * @param value option value
     */
    ScopedPetscOption( std::string name, std::string value )
        :
        M_name( std::move( name ) )
    {
        checkPetscSetup( PetscOptionsSetValue( nullptr, M_name.c_str(), value.c_str() ),
                         "PetscOptionsSetValue" );
    }

    /** Clear the registered PETSc option without throwing. */
    ~ScopedPetscOption()
    {
        if ( auto const error = PetscOptionsClearValue( nullptr, M_name.c_str() ) )
            LOG( ERROR ) << "PetscOptionsClearValue failed with PETSc error " << error;
    }

    /** Option guards are unique owners of their cleanup action. */
    ScopedPetscOption( ScopedPetscOption const& ) = delete;

    /** Option guards cannot be copy-assigned. */
    ScopedPetscOption& operator=( ScopedPetscOption const& ) = delete;

private:
    std::string M_name; ///< Complete PETSc option key removed at destruction.
};

/**
 * Fill an already allocated square sparse matrix with the identity operator.
 *
 * Existing nonzero storage is retained so repeated Hessian evaluations reuse
 * the same matrix structure.
 *
 * @param matrix matrix whose locally owned diagonal entries are set to one
 */
void
fillIdentity( MatrixSparse<double>& matrix )
{
    matrix.zero();
    for ( auto row = matrix.rowStart(); row < matrix.rowStop(); ++row )
        matrix.set( row, row, 1.0 );
    matrix.close();
}

/**
 * Create an assembled distributed PETSc matrix with the supplied Feel++ map.
 *
 * The returned Feel++ wrapper owns the native matrix, while TAO and Hessian
 * callback views only borrow it.
 *
 * @param map row and column distribution map
 * @param worldComm communicator used to allocate the native matrix
 * @return owning Feel++ matrix initialized with an identity sparsity pattern
 */
std::shared_ptr<MatrixPetsc<double>>
makeDistributedMatrix( std::shared_ptr<DataMap<>> const& map,
                       worldcomm_ptr_t const& worldComm )
{
    Mat native = nullptr;
    checkPetscSetup(
        MatCreateAIJ( worldComm->globalComm(),
                      static_cast<PetscInt>( map->nLocalDofWithoutGhost() ),
                      static_cast<PetscInt>( map->nLocalDofWithoutGhost() ),
                      static_cast<PetscInt>( map->nDof() ),
                      static_cast<PetscInt>( map->nDof() ),
                      1, nullptr, 0, nullptr, &native ),
        "MatCreateAIJ" );
    auto matrix = std::make_shared<MatrixPetsc<double>>( native, map, map, true );

    PetscInt firstRow = 0;
    PetscInt lastRow = 0;
    checkPetscSetup( MatGetOwnershipRange( native, &firstRow, &lastRow ),
                     "MatGetOwnershipRange" );
    for ( PetscInt row = firstRow; row < lastRow; ++row )
    {
        PetscScalar const one = 1.0;
        checkPetscSetup(
            MatSetValues( native, 1, &row, 1, &row, &one, INSERT_VALUES ),
            "MatSetValues" );
    }
    checkPetscSetup( MatAssemblyBegin( native, MAT_FINAL_ASSEMBLY ),
                     "MatAssemblyBegin" );
    checkPetscSetup( MatAssemblyEnd( native, MAT_FINAL_ASSEMBLY ),
                     "MatAssemblyEnd" );

    return matrix;
}

/**
 * Check that a solve rethrows the original callback exception.
 *
 * @tparam SolveCallable nullary callable invoking the solve
 * @param solve callable expected to throw
 * @param expectedMessage exact callback exception message
 */
template<typename SolveCallable>
void
checkCallbackException( SolveCallable&& solve, std::string const& expectedMessage )
{
    bool caught = false;
    try
    {
        std::invoke( std::forward<SolveCallable>( solve ) );
    }
    catch ( std::runtime_error const& error )
    {
        caught = true;
        BOOST_CHECK_EQUAL( error.what(), expectedMessage );
    }
    BOOST_CHECK( caught );
}

} // namespace

BOOST_AUTO_TEST_SUITE( tao_suite )

/** Verify that prefixed optimization options are registered only on construction. */
BOOST_AUTO_TEST_CASE( feelpp_options_are_registered_lazily )
{
    std::string const prefix = "tao_lazy_options";
    std::string const typeOption = prefixvm( prefix, "tao-type" );
    std::string const maximumIterationsOption = prefixvm( prefix, "tao-maxit" );

    BOOST_CHECK_EQUAL( Environment::vm().count( typeOption ), 0 );
    BOOST_CHECK(
        Environment::optionsDescription().find_nothrow( typeOption, false ) == nullptr );

    auto options = solveroptimization_options( prefix );
    po::variables_map values;
    std::vector<std::string> arguments{
        "--" + typeOption + "=cg",
        "--" + maximumIterationsOption + "=17",
        "--" + prefixvm( prefix, "tao-maxfcn" ) + "=31",
        "--" + prefixvm( prefix, "tao-gatol" ) + "=1e-11",
        "--" + prefixvm( prefix, "tao-grtol" ) + "=2e-11",
        "--" + prefixvm( prefix, "tao-gttol" ) + "=3e-11",
        "--" + prefixvm( prefix, "tao-steptol" ) + "=4e-12",
        "--" + prefixvm( prefix, "tao-monitor" ) + "=true",
        "--" + prefixvm( prefix, "tao-converged-reason" ) + "=true",
        "--" + prefixvm( prefix, "tao-view" ) + "=true"
    };
    po::store( po::command_line_parser( arguments ).options( options ).run(), values );
    po::notify( values );

    auto solver = optimizationSolver<double>(
        _backend = "petsc", _name = prefix,
        _worldcomm = Environment::worldCommPtr(), _vm = values );

    BOOST_REQUIRE( solver );
    BOOST_CHECK_EQUAL( Environment::vm().count( typeOption ), 1 );
    BOOST_CHECK(
        Environment::optionsDescription().find_nothrow( typeOption, false ) != nullptr );
    BOOST_CHECK_EQUAL( solver->type(), "cg" );
    BOOST_CHECK_EQUAL( values.at( maximumIterationsOption ).as<int>(), 17 );
    BOOST_CHECK_EQUAL(
        values.at( prefixvm( prefix, "tao-maxfcn" ) ).as<int>(), 31 );
    BOOST_CHECK_CLOSE(
        values.at( prefixvm( prefix, "tao-gatol" ) ).as<double>(),
        1e-11, 1e-10 );
    BOOST_CHECK_CLOSE(
        values.at( prefixvm( prefix, "tao-steptol" ) ).as<double>(),
        4e-12, 1e-10 );
    BOOST_CHECK( solver->monitorEnabled() );
    BOOST_CHECK( solver->convergedReasonEnabled() );
    BOOST_CHECK( solver->viewEnabled() );
}

/** Verify that a combined callback is preferred over separately registered callbacks. */
BOOST_AUTO_TEST_CASE( combined_objective_gradient_precedence )
{
    auto const worldComm = Environment::worldCommPtr();
    uint32_type const localSize = 4;
    double constexpr minimizer = 2.5;

    auto variable = std::make_shared<VectorPetscMPI<double>>(
        makeContiguousMap( localSize, worldComm ) );
    variable->zero();

    SolverOptimizationPetsc<double> solver( "tao_quadratic", worldComm );
    solver.setType( "lmvm" );
    solver.setMaximumIterations( 50 );
    solver.setGradientTolerances( 1e-10, 1e-10, 1e-10 );

    std::size_t objectiveEvaluations = 0;
    std::size_t gradientEvaluations = 0;
    std::size_t combinedEvaluations = 0;
    solver.setObjective(
        [&]( Vector<double> const& state )
        {
            ++objectiveEvaluations;
            return quadraticObjective( state, minimizer );
        } );
    solver.setGradient(
        [&]( Vector<double> const& state, Vector<double>& gradient )
        {
            ++gradientEvaluations;
            quadraticGradient( state, minimizer, gradient );
        } );
    solver.setObjectiveGradient(
        [&]( Vector<double> const& state, Vector<double>& gradient )
        {
            ++combinedEvaluations;
            quadraticGradient( state, minimizer, gradient );
            return 0.5 * gradient.dot( gradient );
        } );

    auto const result = solver.solve( variable );

    BOOST_CHECK( result.converged );
    BOOST_CHECK( result.status == OptimizationStatus::Converged );
    BOOST_CHECK_GT( result.rawReason, 0 );
    BOOST_CHECK_EQUAL( objectiveEvaluations, 0 );
    BOOST_CHECK_EQUAL( gradientEvaluations, 0 );
    BOOST_CHECK_GT( combinedEvaluations, 0 );
    BOOST_CHECK_GT( result.functionEvaluations, 0 );
    BOOST_CHECK_EQUAL( result.solverType, "lmvm" );
    BOOST_CHECK_EQUAL( result.optionsPrefix, "tao_quadratic_" );
    BOOST_CHECK_SMALL( result.objective, 1e-16 );
    BOOST_CHECK_SMALL( result.gradientNorm, 1e-8 );

    variable->add( -minimizer );
    BOOST_CHECK_SMALL( variable->l2Norm(), 1e-8 );

    auto const combinedEvaluationsBeforeReset = combinedEvaluations;
    solver.resetObjectiveGradient();
    auto const separateResult = solver.solve( variable );

    BOOST_CHECK( separateResult.converged );
    BOOST_CHECK_GT( objectiveEvaluations, 0 );
    BOOST_CHECK_GT( gradientEvaluations, 0 );
    BOOST_CHECK_EQUAL( combinedEvaluations, combinedEvaluationsBeforeReset );

    variable->add( -minimizer );
    BOOST_CHECK_SMALL( variable->l2Norm(), 1e-8 );
}

/** Verify separate objective and gradient callbacks with distributed borrowed views. */
BOOST_AUTO_TEST_CASE( separate_objective_and_gradient )
{
    auto const worldComm = Environment::worldCommPtr();
    uint32_type const localSize = 4;
    double constexpr minimizer = -1.75;

    auto variable = std::make_shared<VectorPetscMPI<double>>(
        makeContiguousMap( localSize, worldComm ) );
    variable->zero();

    SolverOptimizationPetsc<double> solver( "tao_separate", worldComm );
    solver.setType( "lmvm" );
    solver.setMaximumIterations( 50 );
    solver.setGradientTolerances( 1e-10, 1e-10, 1e-10 );

    std::size_t objectiveEvaluations = 0;
    std::size_t gradientEvaluations = 0;
    solver.setObjective(
        [&]( Vector<double> const& state )
        {
            ++objectiveEvaluations;
            return quadraticObjective( state, minimizer );
        } );
    solver.setGradient(
        [&]( Vector<double> const& state, Vector<double>& gradient )
        {
            ++gradientEvaluations;
            quadraticGradient( state, minimizer, gradient );
        } );

    auto const result = solver.solve( variable );

    BOOST_CHECK( result.converged );
    BOOST_CHECK( result.status == OptimizationStatus::Converged );
    BOOST_CHECK_GT( objectiveEvaluations, 0 );
    BOOST_CHECK_GT( gradientEvaluations, 0 );
    BOOST_CHECK_GT( result.functionEvaluations, 0 );
    BOOST_CHECK_SMALL( result.objective, 1e-12 );
    BOOST_CHECK_SMALL( result.gradientNorm, 1e-8 );

    variable->add( -minimizer );
    BOOST_CHECK_SMALL( variable->l2Norm(), 1e-8 );
}

/** Verify exact Hessian assembly, matrix reuse, and shared/separate preconditioners. */
BOOST_AUTO_TEST_CASE( distributed_exact_hessian )
{
    using solver_type = SolverOptimizationPetsc<double>;
    using objective_gradient_data_type = solver_type::objective_gradient_data_type;
    using hessian_data_type = solver_type::hessian_data_type;
    auto const worldComm = Environment::worldCommPtr();
    uint32_type const localSize = 4;
    double constexpr minimizer = 1.5;

    auto map = makeContiguousMap( localSize, worldComm );
    auto variable = std::make_shared<VectorPetscMPI<double>>( map );
    variable->zero();
    auto hessian = makeDistributedMatrix( map, worldComm );
    auto preconditioner = makeDistributedMatrix( map, worldComm );
    Mat const nativeHessian = hessian->mat();
    Mat const nativePreconditioner = preconditioner->mat();

    solver_type solver( "tao_exact_hessian", worldComm );
    solver.setType( "nls" );
    solver.setMaximumIterations( 20 );
    solver.setGradientTolerances( 1e-10, 1e-10, 1e-10 );
    solver.setObjectiveGradient(
        []( objective_gradient_data_type& data )
        {
            quadraticGradient( data.state, minimizer, data.gradient );
            data.objective = 0.5 * data.gradient.dot( data.gradient );
        } );

    std::size_t separateHessianEvaluations = 0;
    bool reusedSeparateMatrices = true;
    bool receivedDistinctMatrices = true;
    solver.setHessian(
        hessian, preconditioner,
        [&]( hessian_data_type& data )
        {
            ++separateHessianEvaluations;
            auto* petscHessian = dynamic_cast<MatrixPetsc<double>*>( &data.hessian );
            auto* petscPreconditioner =
                dynamic_cast<MatrixPetsc<double>*>( &data.preconditioner );
            reusedSeparateMatrices =
                reusedSeparateMatrices && petscHessian && petscPreconditioner &&
                petscHessian->mat() == nativeHessian &&
                petscPreconditioner->mat() == nativePreconditioner;
            receivedDistinctMatrices = receivedDistinctMatrices && !data.sharesMatrix();
            fillIdentity( data.hessian );
            fillIdentity( data.preconditioner );
        } );

    auto const separateResult = solver.solve( variable );
    BOOST_CHECK( separateResult.converged );
    BOOST_CHECK_EQUAL( separateResult.solverType, "nls" );
    BOOST_CHECK( solver.nativeTao() != nullptr );
    BOOST_CHECK( solver.nativeKsp() != nullptr );
    BOOST_CHECK_GT( separateHessianEvaluations, 0 );
    BOOST_CHECK( reusedSeparateMatrices );
    BOOST_CHECK( receivedDistinctMatrices );
    BOOST_CHECK_EQUAL( hessian->mat(), nativeHessian );
    BOOST_CHECK_EQUAL( preconditioner->mat(), nativePreconditioner );
    variable->add( -minimizer );
    BOOST_CHECK_SMALL( variable->l2Norm(), 1e-8 );

    variable->setConstant( 0.0 );
    std::size_t sharedHessianEvaluations = 0;
    bool receivedAliasedMatrix = true;
    solver.setHessian(
        hessian,
        [&]( hessian_data_type& data )
        {
            ++sharedHessianEvaluations;
            receivedAliasedMatrix = receivedAliasedMatrix && data.sharesMatrix();
            fillIdentity( data.hessian );
        } );

    auto const sharedResult = solver.solve( variable );
    BOOST_CHECK( sharedResult.converged );
    BOOST_CHECK_GT( sharedHessianEvaluations, 0 );
    BOOST_CHECK( receivedAliasedMatrix );
    BOOST_CHECK_EQUAL( hessian->mat(), nativeHessian );
    variable->add( -minimizer );
    BOOST_CHECK_SMALL( variable->l2Norm(), 1e-8 );

    auto const evaluationsBeforeReset = sharedHessianEvaluations;
    variable->setConstant( 0.0 );
    solver.resetHessian();
    solver.setType( "lmvm" );
    auto const resultWithoutHessian = solver.solve( variable );
    BOOST_CHECK( resultWithoutHessian.converged );
    BOOST_CHECK_EQUAL( resultWithoutHessian.solverType, "lmvm" );
    BOOST_CHECK_EQUAL( sharedHessianEvaluations, evaluationsBeforeReset );
}

/** Verify C++ and option-driven PETSc monitors, registration control, and native TAO access. */
BOOST_AUTO_TEST_CASE( cpp_monitors_and_native_escape_hatch )
{
    using solver_type = SolverOptimizationPetsc<double>;
    using monitor_record_type = solver_type::monitor_record_type;
    auto const worldComm = Environment::worldCommPtr();
    uint32_type const localSize = 4;
    double constexpr minimizer = 2.0;
    ScopedPetscOption petscMonitor( "-tao_cpp_monitor_tao_monitor", "" );

    auto variable = makeDistributedVector( localSize, worldComm );
    solver_type solver( "tao_cpp_monitor", worldComm );
    solver.setType( "lmvm" );
    solver.setMaximumIterations( 50 );
    solver.setGradientTolerances( 1e-10, 1e-10, 1e-10 );
    solver.setObjectiveGradient(
        []( Vector<double> const& state, Vector<double>& gradient )
        {
            quadraticGradient( state, minimizer, gradient );
            return 0.5 * gradient.dot( gradient );
        } );

    std::vector<monitor_record_type> records;
    std::size_t removedMonitorCalls = 0;
    auto const retainedMonitor = solver.addMonitor(
        [&]( monitor_record_type const& record ) { records.push_back( record ); } );
    auto const removedMonitor = solver.addMonitor(
        [&]( monitor_record_type const& ) { ++removedMonitorCalls; } );

    BOOST_CHECK_NE( retainedMonitor, removedMonitor );
    BOOST_CHECK_EQUAL( solver.monitorCount(), 2 );
    BOOST_CHECK( solver.removeMonitor( removedMonitor ) );
    BOOST_CHECK( !solver.removeMonitor( removedMonitor ) );
    BOOST_CHECK_EQUAL( solver.monitorCount(), 1 );
    BOOST_CHECK( solver.nativeTao() != nullptr );

    auto const result = solver.solve( variable );
    BOOST_CHECK( result.converged );
    BOOST_CHECK( !records.empty() );
    BOOST_CHECK_EQUAL( removedMonitorCalls, 0 );
    for ( std::size_t index = 0; index < records.size(); ++index )
    {
        BOOST_CHECK( std::isfinite( records[index].objective ) );
        BOOST_CHECK( std::isfinite( records[index].gradientNorm ) );
        if ( index > 0 )
            BOOST_CHECK_GE( records[index].iteration, records[index - 1].iteration );
    }

    auto const recordsBeforeReset = records.size();
    solver.resetMonitors();
    BOOST_CHECK( !solver.hasMonitors() );
    variable->setConstant( 0.0 );
    auto const resultWithoutMonitors = solver.solve( variable );
    BOOST_CHECK( resultWithoutMonitors.converged );
    BOOST_CHECK_EQUAL( records.size(), recordsBeforeReset );

    solver.clear();
    BOOST_CHECK( solver.nativeTao() == nullptr );
}

/** Verify PETSc runtime options remain isolated by each solver prefix. */
BOOST_AUTO_TEST_CASE( options_prefix_isolation )
{
    auto const worldComm = Environment::worldCommPtr();
    uint32_type const localSize = 4;
    double constexpr minimizer = -1.25;
    ScopedPetscOption firstType( "-tao_prefix_first_tao_type", "cg" );

    auto firstVariable = makeDistributedVector( localSize, worldComm );
    auto secondVariable = makeDistributedVector( localSize, worldComm );
    SolverOptimizationPetsc<double> firstSolver( "tao_prefix_first", worldComm );
    SolverOptimizationPetsc<double> secondSolver( "tao_prefix_second", worldComm );
    firstSolver.setType( "lmvm" );
    secondSolver.setType( "lmvm" );

    auto objectiveGradient = []( Vector<double> const& state, Vector<double>& gradient )
    {
        quadraticGradient( state, minimizer, gradient );
        return 0.5 * gradient.dot( gradient );
    };
    firstSolver.setObjectiveGradient( objectiveGradient );
    secondSolver.setObjectiveGradient( objectiveGradient );

    auto const firstResult = firstSolver.solve( firstVariable );
    auto const secondResult = secondSolver.solve( secondVariable );
    BOOST_CHECK( firstResult.converged );
    BOOST_CHECK( secondResult.converged );
    BOOST_CHECK_EQUAL( firstResult.solverType, "cg" );
    BOOST_CHECK_EQUAL( secondResult.solverType, "lmvm" );
    BOOST_CHECK_EQUAL( firstResult.optionsPrefix, "tao_prefix_first_" );
    BOOST_CHECK_EQUAL( secondResult.optionsPrefix, "tao_prefix_second_" );
}

/** Verify divergence reasons are normalized without treating TaoSolve success as convergence. */
BOOST_AUTO_TEST_CASE( result_reason_normalization )
{
    auto const worldComm = Environment::worldCommPtr();
    uint32_type constexpr localSize = 4;
    auto variable = makeDistributedVector( localSize, worldComm );
    auto minimizer = makeDistributedVector( localSize, worldComm, 3.0 );
    auto weights = makeDistributedVector( localSize, worldComm );
    for ( uint32_type index = 0; index < localSize; ++index )
        weights->set( index, static_cast<double>( index + 1 ) );
    weights->close();

    SolverOptimizationPetsc<double> solver( "tao_reason", worldComm );
    solver.setType( "lmvm" );
    solver.setMaximumIterations( 1 );
    solver.setObjectiveGradient(
        [minimizer, weights]( Vector<double> const& state, Vector<double>& gradient )
        {
            auto displacement = state.clone();
            displacement->add( -1.0, *minimizer );
            gradient.pointwiseMult( *weights, *displacement );
            return 0.5 * gradient.dot( *displacement );
        } );

    auto const result = solver.solve( variable );
    BOOST_CHECK( !result.converged );
    BOOST_CHECK( result.status == OptimizationStatus::DivergedLineSearch );
    BOOST_CHECK_EQUAL( result.rawReason, static_cast<int>( TAO_DIVERGED_LS_FAILURE ) );
    BOOST_CHECK( !result.reason.empty() );
    BOOST_CHECK( !result.diagnostic.empty() );
}

/** Verify named construction, generic-base use, callback replacement, and solver reuse. */
BOOST_AUTO_TEST_CASE( generic_interface_factory_and_repeated_solves )
{
    using generic_solver_type = SolverOptimization<double>;
    static_assert(
        std::derived_from<SolverOptimizationPetsc<double>, generic_solver_type> );

    auto const worldComm = Environment::worldCommPtr();
    uint32_type constexpr localSize = 4;
    auto variable = makeDistributedVector( localSize, worldComm, 4.0 );
    auto solver = optimizationSolver<double>(
        _backend = "petsc", _name = "tao_generic", _worldcomm = worldComm );

    BOOST_REQUIRE( solver );
    BOOST_CHECK_EQUAL( solver->prefix(), "tao_generic" );
    solver->setType( "lmvm" );
    solver->setMaximumIterations( 50 );
    solver->setGradientTolerances( 1e-10, 1e-10, 1e-10 );
    BOOST_CHECK_THROW( solver->setStepTolerance( -1.0 ), std::invalid_argument );
    solver->setStepTolerance( 1e-14 );

    std::size_t firstCallbackCalls = 0;
    double constexpr firstMinimizer = 1.5;
    solver->setObjectiveGradient(
        [&]( Vector<double> const& state, Vector<double>& gradient )
        {
            ++firstCallbackCalls;
            quadraticGradient( state, firstMinimizer, gradient );
            return 0.5 * gradient.dot( gradient );
        } );

    auto const firstResult = solver->solve( variable );
    BOOST_CHECK( firstResult.converged );
    BOOST_CHECK_GT( firstCallbackCalls, 0 );
    variable->add( -firstMinimizer );
    BOOST_CHECK_SMALL( variable->l2Norm(), 1e-8 );

    auto const firstCallsAfterSolve = firstCallbackCalls;
    std::size_t secondCallbackCalls = 0;
    double constexpr secondMinimizer = -0.75;
    variable->setConstant( 3.0 );
    solver->setObjectiveGradient(
        [&]( Vector<double> const& state, Vector<double>& gradient )
        {
            ++secondCallbackCalls;
            quadraticGradient( state, secondMinimizer, gradient );
            return 0.5 * gradient.dot( gradient );
        } );

    auto const secondResult = solver->solve( variable );
    BOOST_CHECK( secondResult.converged );
    BOOST_CHECK_GT( secondCallbackCalls, 0 );
    BOOST_CHECK_EQUAL( firstCallbackCalls, firstCallsAfterSolve );
    BOOST_CHECK_EQUAL( secondResult.optionsPrefix, "tao_generic_" );
    variable->add( -secondMinimizer );
    BOOST_CHECK_SMALL( variable->l2Norm(), 1e-8 );

    BOOST_CHECK_THROW(
        static_cast<void>( generic_solver_type::build(
            "unsupported", "tao_invalid", worldComm ) ),
        std::invalid_argument );
}

/** Verify a PETSc solver rejects a non-PETSc Feel++ vector with a clear diagnostic. */
BOOST_AUTO_TEST_CASE( non_petsc_vector_rejection )
{
    auto const worldComm = Environment::worldCommPtr();
    auto solver = optimizationSolver<double>(
        _backend = "petsc", _name = "tao_non_petsc_vector",
        _worldcomm = worldComm );
    solver->setObjectiveGradient(
        []( Vector<double> const& state, Vector<double>& gradient )
        {
            quadraticGradient( state, 0.0, gradient );
            return 0.5 * gradient.dot( gradient );
        } );

    VectorUblas<double> incompatible( 4 );
    try
    {
        static_cast<void>( solver->solve( incompatible ) );
        BOOST_FAIL( "A PETSc TAO solver accepted a non-PETSc Feel++ vector" );
    }
    catch ( std::invalid_argument const& error )
    {
        BOOST_CHECK_NE(
            std::string( error.what() ).find( "must be VectorPetsc-compatible" ),
            std::string::npos );
    }
}

/** Verify two-sided, one-sided, constant, and component-wise variable bounds. */
BOOST_AUTO_TEST_CASE( distributed_variable_bounds )
{
    using solver_type = SolverOptimizationPetsc<double>;
    auto const worldComm = Environment::worldCommPtr();
    uint32_type const localSize = 4;

    {
        double constexpr minimizer = 2.5;
        double constexpr expected = 1.25;
        auto variable = makeDistributedVector( localSize, worldComm );
        auto lower = makeDistributedVector( localSize, worldComm, -1.0 );
        auto upper = makeDistributedVector( localSize, worldComm, expected );

        solver_type solver( "tao_box_bounds", worldComm );
        solver.setType( "blmvm" );
        solver.setMaximumIterations( 100 );
        solver.setGradientTolerances( 1e-10, 1e-10, 1e-10 );
        solver.setVariableBounds( lower, upper );
        solver.setObjectiveGradient(
            []( Vector<double> const& state, Vector<double>& gradient )
            {
                quadraticGradient( state, minimizer, gradient );
                return 0.5 * gradient.dot( gradient );
            } );

        auto const result = solver.solve( variable );
        BOOST_CHECK( result.converged );
        BOOST_CHECK_EQUAL( result.solverType, "blmvm" );
        variable->add( -expected );
        BOOST_CHECK_SMALL( variable->l2Norm(), 1e-8 );
    }

    {
        double constexpr minimizer = -2.0;
        double constexpr expected = -0.5;
        auto variable = makeDistributedVector( localSize, worldComm );

        solver_type solver( "tao_lower_bound", worldComm );
        solver.setType( "blmvm" );
        solver.setMaximumIterations( 100 );
        solver.setGradientTolerances( 1e-10, 1e-10, 1e-10 );
        solver.setConstantLowerBound( expected );
        solver.setObjectiveGradient(
            []( Vector<double> const& state, Vector<double>& gradient )
            {
                quadraticGradient( state, minimizer, gradient );
                return 0.5 * gradient.dot( gradient );
            } );

        auto const result = solver.solve( variable );
        BOOST_CHECK( result.converged );
        variable->add( -expected );
        BOOST_CHECK_SMALL( variable->l2Norm(), 1e-8 );
    }

    {
        double constexpr minimizer = 2.0;
        double constexpr finiteUpper = 0.75;
        auto variable = makeDistributedVector( localSize, worldComm );
        auto upper = makeDistributedVector( localSize, worldComm, finiteUpper );
        upper->set( 0, solver_type::positiveInfinity() );
        upper->close();

        solver_type solver( "tao_upper_bound", worldComm );
        solver.setType( "blmvm" );
        solver.setMaximumIterations( 100 );
        solver.setGradientTolerances( 1e-10, 1e-10, 1e-10 );
        solver.setUpperBound( upper );
        solver.setObjectiveGradient(
            []( Vector<double> const& state, Vector<double>& gradient )
            {
                quadraticGradient( state, minimizer, gradient );
                return 0.5 * gradient.dot( gradient );
            } );

        auto const result = solver.solve( variable );
        BOOST_CHECK( result.converged );

        auto expected = makeDistributedVector( localSize, worldComm, finiteUpper );
        expected->set( 0, minimizer );
        expected->close();
        variable->add( -1.0, *expected );
        BOOST_CHECK_SMALL( variable->l2Norm(), 1e-8 );
    }
}

/** Verify exception transport for all objective, derivative, and monitor callbacks. */
BOOST_AUTO_TEST_CASE( callback_exceptions_are_rethrown )
{
    auto const worldComm = Environment::worldCommPtr();
    uint32_type const localSize = 4;

    {
        auto variable = std::make_shared<VectorPetscMPI<double>>(
            makeContiguousMap( localSize, worldComm ) );
        SolverOptimizationPetsc<double> solver( "tao_objective_exception", worldComm );
        solver.setObjective(
            []( Vector<double> const& ) -> double
            {
                throw std::runtime_error( "objective callback failure" );
            } );
        solver.setGradient(
            []( Vector<double> const& state, Vector<double>& gradient )
            {
                quadraticGradient( state, 1.0, gradient );
            } );
        checkCallbackException(
            [&]() { static_cast<void>( solver.solve( variable ) ); },
            "objective callback failure" );
    }

    {
        auto variable = std::make_shared<VectorPetscMPI<double>>(
            makeContiguousMap( localSize, worldComm ) );
        SolverOptimizationPetsc<double> solver( "tao_gradient_exception", worldComm );
        solver.setObjective(
            []( Vector<double> const& state )
            {
                return quadraticObjective( state, 1.0 );
            } );
        solver.setGradient(
            []( Vector<double> const&, Vector<double>& )
            {
                throw std::runtime_error( "gradient callback failure" );
            } );
        checkCallbackException(
            [&]() { static_cast<void>( solver.solve( variable ) ); },
            "gradient callback failure" );
    }

    {
        auto variable = std::make_shared<VectorPetscMPI<double>>(
            makeContiguousMap( localSize, worldComm ) );
        SolverOptimizationPetsc<double> solver( "tao_combined_exception", worldComm );
        solver.setObjectiveGradient(
            []( Vector<double> const&, Vector<double>& ) -> double
            {
                throw std::runtime_error( "combined callback failure" );
            } );
        checkCallbackException(
            [&]() { static_cast<void>( solver.solve( variable ) ); },
            "combined callback failure" );
    }

    {
        auto map = makeContiguousMap( localSize, worldComm );
        auto variable = std::make_shared<VectorPetscMPI<double>>( map );
        auto hessian = makeDistributedMatrix( map, worldComm );
        SolverOptimizationPetsc<double> solver( "tao_hessian_exception", worldComm );
        solver.setType( "nls" );
        solver.setObjectiveGradient(
            []( Vector<double> const& state, Vector<double>& gradient )
            {
                quadraticGradient( state, 1.0, gradient );
                return 0.5 * gradient.dot( gradient );
            } );
        solver.setHessian(
            hessian,
            []( Vector<double> const&, MatrixSparse<double>&,
                MatrixSparse<double>& )
            {
                throw std::runtime_error( "Hessian callback failure" );
            } );
        checkCallbackException(
            [&]() { static_cast<void>( solver.solve( variable ) ); },
            "Hessian callback failure" );
    }
    {
        auto variable = makeDistributedVector( localSize, worldComm );
        SolverOptimizationPetsc<double> solver( "tao_monitor_exception", worldComm );
        solver.setObjectiveGradient(
            []( Vector<double> const& state, Vector<double>& gradient )
            {
                quadraticGradient( state, 1.0, gradient );
                return 0.5 * gradient.dot( gradient );
            } );
        static_cast<void>( solver.addMonitor(
            []( SolverOptimizationPetsc<double>::monitor_record_type const& )
            {
                throw std::runtime_error( "monitor callback failure" );
            } ) );
        checkCallbackException(
            [&]() { static_cast<void>( solver.solve( variable ) ); },
            "monitor callback failure" );
    }
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace Feel
