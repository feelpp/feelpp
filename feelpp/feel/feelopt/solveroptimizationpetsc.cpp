/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <sstream>
#include <stdexcept>
#include <type_traits>

#include <feel/feelopt/solveroptimizationpetsc.hpp>

#if defined( FEELPP_HAS_PETSC_TAO )

namespace Feel
{

template<typename T, typename SizeT>
SolverOptimizationPetsc<T, SizeT>::SolverOptimizationPetsc(
    std::string const& prefix, worldcomm_ptr_t const& worldComm,
    po::variables_map const& vm )
    :
    super( prefix, worldComm, vm ),
    M_petscPrefix( normalizePrefix( prefix ) )
{
    static_assert( std::is_same_v<T, PetscScalar>,
                   "SolverOptimizationPetsc requires T to match PetscScalar" );
    static_assert( std::is_same_v<real_type, PetscReal>,
                   "SolverOptimizationPetsc requires a real PETSc scalar configuration" );

    if ( !worldComm )
        throw std::invalid_argument( "SolverOptimizationPetsc requires a valid WorldComm" );

    M_context.solver = this;
    try
    {
        checkPetsc( TaoCreate( this->worldComm().globalComm(), &M_tao ), "TaoCreate" );
        checkPetsc( TaoSetOptionsPrefix( M_tao, M_petscPrefix.c_str() ),
                    "TaoSetOptionsPrefix" );
    }
    catch ( ... )
    {
        if ( M_tao )
            TaoDestroy( &M_tao );
        throw;
    }
}

template<typename T, typename SizeT>
SolverOptimizationPetsc<T, SizeT>::~SolverOptimizationPetsc()
{
    this->clear();
}

template<typename T, typename SizeT>
void
SolverOptimizationPetsc<T, SizeT>::clear() noexcept
{
    M_context.worldComm.reset();
    M_context.map.reset();
    M_context.variable = nullptr;
    M_context.hessian.reset();
    M_context.preconditioner.reset();
    M_context.exception = {};
    M_hessianConfigured = false;

    if ( !M_tao )
    {
        M_materializedLowerBound.reset();
        M_materializedUpperBound.reset();
        return;
    }

    PetscErrorCode const error = TaoDestroy( &M_tao );
    if ( error )
        LOG( ERROR ) << "TaoDestroy failed with PETSc error " << error;
    M_tao = nullptr;
    M_materializedLowerBound.reset();
    M_materializedUpperBound.reset();
}

template<typename T, typename SizeT>
void
SolverOptimizationPetsc<T, SizeT>::setFromOptions()
{
    this->configureMonitors();
    this->applyProgrammaticConfiguration();

    std::vector<std::string> injectedOptions;
    auto injectDefaultOption = [&]( std::string const& suffix,
                                    std::string const& value )
    {
        std::string const option = "-" + M_petscPrefix + suffix;
        PetscBool runtimeOptionExists = PETSC_FALSE;
        checkPetsc( PetscOptionsHasName( nullptr, nullptr, option.c_str(),
                                         &runtimeOptionExists ),
                    "PetscOptionsHasName(optimization option)" );
        if ( runtimeOptionExists )
            return;

        checkPetsc( PetscOptionsSetValue( nullptr, option.c_str(), value.c_str() ),
                    "PetscOptionsSetValue(optimization option)" );
        injectedOptions.push_back( option );
    };

    if ( auto const tolerance = this->stepTolerance() )
    {
        std::ostringstream value;
        value.precision( std::numeric_limits<real_type>::max_digits10 );
        value << *tolerance;
        injectDefaultOption( "tao_steptol", value.str() );
    }
    if ( this->monitorEnabled() )
        injectDefaultOption( "tao_monitor", "" );
    if ( this->convergedReasonEnabled() )
        injectDefaultOption( "tao_converged_reason", "" );
    if ( this->viewEnabled() )
        injectDefaultOption( "tao_view", "" );

    try
    {
        checkPetsc( TaoSetFromOptions( M_tao ), "TaoSetFromOptions" );
    }
    catch ( ... )
    {
        for ( auto const& option : injectedOptions )
        {
            if ( auto const error = PetscOptionsClearValue( nullptr, option.c_str() ) )
                LOG( ERROR ) << "PetscOptionsClearValue failed with PETSc error " << error;
        }
        throw;
    }

    for ( auto const& option : injectedOptions )
        checkPetsc( PetscOptionsClearValue( nullptr, option.c_str() ),
                    "PetscOptionsClearValue(optimization option)" );
}

template<typename T, typename SizeT>
typename SolverOptimizationPetsc<T, SizeT>::result_type
SolverOptimizationPetsc<T, SizeT>::solve( vector_type& variable )
{
    if ( !M_tao )
        throw std::logic_error( "SolverOptimizationPetsc has been cleared" );
    if ( !this->hasObjectiveGradient() && !this->hasObjective() )
        throw std::invalid_argument(
            "SolverOptimizationPetsc requires an objective callback when no combined "
            "objective-gradient callback is registered" );
    if ( !this->hasObjectiveGradient() && !this->hasGradient() )
        throw std::invalid_argument(
            "SolverOptimizationPetsc requires a gradient callback when no combined "
            "objective-gradient callback is registered" );

    Vec const solution = this->nativeVector( variable, "optimization variable" );

    M_context.worldComm = variable.worldCommPtr();
    M_context.map = variable.mapPtr();
    M_context.variable = &variable;
    M_context.hessian = this->hessianMatrix();
    M_context.preconditioner = this->hessianPreconditionerMatrix();
    M_context.exception = {};
    /** Reset borrowed per-solve state on every exit path. */
    struct ContextReset
    {
        CallbackContext& context; ///< Callback state whose retained solve data is released.

        /** Release all retained algebraic solve data after the solve. */
        ~ContextReset()
        {
            context.worldComm.reset();
            context.map.reset();
            context.variable = nullptr;
            context.hessian.reset();
            context.preconditioner.reset();
        }
    } reset{ M_context };

    checkPetsc( TaoSetSolution( M_tao, solution ), "TaoSetSolution" );
    this->configureVariableBounds( variable );
    this->configureCallbacks();
    this->configureHessian( variable );

    this->setFromOptions();

    checkPetsc( PetscPushErrorHandler( PetscReturnErrorHandler, nullptr ),
                "PetscPushErrorHandler" );
    PetscErrorCode const solveError = TaoSolve( M_tao );
    PetscErrorCode const popError = PetscPopErrorHandler();

    if ( M_context.exception )
        std::rethrow_exception( M_context.exception );
    checkPetsc( popError, "PetscPopErrorHandler" );
    checkPetsc( solveError, "TaoSolve" );

    PetscInt iterations = 0;
    PetscInt functionEvaluations = 0;
    PetscReal objective = 0;
    PetscReal gradientNorm = 0;
    PetscReal constraintNorm = 0;
    PetscReal stepNorm = 0;
    TaoConvergedReason reason = TAO_CONTINUE_ITERATING;
    TaoType taoType = nullptr;

    checkPetsc( TaoGetSolutionStatus( M_tao, &iterations, &objective, &gradientNorm,
                                      &constraintNorm, &stepNorm, &reason ),
                "TaoGetSolutionStatus" );
    checkPetsc( TaoGetCurrentFunctionEvaluations( M_tao, &functionEvaluations ),
                "TaoGetCurrentFunctionEvaluations" );
    checkPetsc( TaoGetType( M_tao, &taoType ), "TaoGetType" );

    result_type result;
    result.converged = reason > 0;
    result.status = normalizeReason( reason );
    result.rawReason = static_cast<int>( reason );
    result.reason = reasonString( reason );
    result.iterations = static_cast<size_type>( iterations );
    result.functionEvaluations = static_cast<size_type>( functionEvaluations );
    result.objective = static_cast<real_type>( objective );
    result.gradientNorm = static_cast<real_type>( gradientNorm );
    result.constraintNorm = static_cast<real_type>( constraintNorm );
    result.stepNorm = static_cast<real_type>( stepNorm );
    result.solverType = taoType ? taoType : "";
    result.optionsPrefix = M_petscPrefix;
    if ( !result.converged )
    {
        std::ostringstream diagnostic;
        diagnostic << "TAO solver '" << result.solverType << "' with prefix '"
                   << result.optionsPrefix << "' stopped after " << result.iterations
                   << " iterations: " << result.reason << " (raw reason "
                   << result.rawReason << ")";
        result.diagnostic = diagnostic.str();
    }
    return result;
}

template<typename T, typename SizeT>
typename SolverOptimizationPetsc<T, SizeT>::vector_ptrtype
SolverOptimizationPetsc<T, SizeT>::makeVectorView(
    Vec vector, CallbackContext const& context )
{
    if ( !context.worldComm || !context.map )
        throw std::logic_error(
            "SolverOptimizationPetsc callback requires a WorldComm and DataMap" );

    if ( context.worldComm->globalSize() > 1 )
        return std::make_shared<VectorPetscMPI<T>>( vector, context.map, false );
    return std::make_shared<VectorPetsc<T>>( vector, context.map, false );
}

template<typename T, typename SizeT>
typename SolverOptimizationPetsc<T, SizeT>::sparse_matrix_ptrtype
SolverOptimizationPetsc<T, SizeT>::makeMatrixView(
    Mat matrix, datamap_ptrtype const& rowMap,
    datamap_ptrtype const& columnMap, CallbackContext const& context )
{
    if ( !context.worldComm || !rowMap || !columnMap )
        throw std::logic_error(
            "SolverOptimizationPetsc matrix callback requires a WorldComm and DataMaps" );

    if ( context.worldComm->globalSize() > 1 )
        return std::make_shared<MatrixPetscMPI<T>>(
            matrix, rowMap, columnMap, false, false );
    return std::make_shared<MatrixPetsc<T>>(
        matrix, rowMap, columnMap, false );
}

template<typename T, typename SizeT>
PetscErrorCode
SolverOptimizationPetsc<T, SizeT>::objectiveBridge(
    Tao, Vec state, PetscReal* objective, void* context ) noexcept
{
    auto& callbackContext = *static_cast<CallbackContext*>( context );
    try
    {
        auto stateView = makeVectorView( state, callbackContext );
        *objective = static_cast<PetscReal>(
            callbackContext.solver->objective()( *stateView ) );
        return PETSC_SUCCESS;
    }
    catch ( ... )
    {
        if ( !callbackContext.exception )
            callbackContext.exception = std::current_exception();
        return PETSC_ERR_USER;
    }
}

template<typename T, typename SizeT>
PetscErrorCode
SolverOptimizationPetsc<T, SizeT>::gradientBridge(
    Tao, Vec state, Vec gradient, void* context ) noexcept
{
    auto& callbackContext = *static_cast<CallbackContext*>( context );
    try
    {
        auto stateView = makeVectorView( state, callbackContext );
        auto gradientView = makeVectorView( gradient, callbackContext );
        callbackContext.solver->gradient()( *stateView, *gradientView );
        if ( !gradientView->closed() )
            gradientView->close();
        return PETSC_SUCCESS;
    }
    catch ( ... )
    {
        if ( !callbackContext.exception )
            callbackContext.exception = std::current_exception();
        return PETSC_ERR_USER;
    }
}

template<typename T, typename SizeT>
PetscErrorCode
SolverOptimizationPetsc<T, SizeT>::objectiveGradientBridge(
    Tao, Vec state, PetscReal* objective, Vec gradient, void* context ) noexcept
{
    auto& callbackContext = *static_cast<CallbackContext*>( context );
    try
    {
        auto stateView = makeVectorView( state, callbackContext );
        auto gradientView = makeVectorView( gradient, callbackContext );
        objective_gradient_data_type data( *stateView, *gradientView );
        callbackContext.solver->objectiveGradient()( data );
        *objective = static_cast<PetscReal>( data.objective );
        if ( !gradientView->closed() )
            gradientView->close();
        return PETSC_SUCCESS;
    }
    catch ( ... )
    {
        if ( !callbackContext.exception )
            callbackContext.exception = std::current_exception();
        return PETSC_ERR_USER;
    }
}

template<typename T, typename SizeT>
PetscErrorCode
SolverOptimizationPetsc<T, SizeT>::hessianBridge(
    Tao, Vec state, Mat hessian, Mat preconditioner, void* context ) noexcept
{
    auto& callbackContext = *static_cast<CallbackContext*>( context );
    try
    {
        auto* solver = callbackContext.solver;
        auto const& hessianStorage = callbackContext.hessian;
        auto const& preconditionerStorage = callbackContext.preconditioner;
        if ( !hessianStorage || !preconditionerStorage )
            throw std::logic_error(
                "SolverOptimizationPetsc Hessian storage disappeared during a solve" );

        auto stateView = makeVectorView( state, callbackContext );
        auto hessianView = makeMatrixView(
            hessian, hessianStorage->mapRowPtr(), hessianStorage->mapColPtr(),
            callbackContext );

        if ( hessian == preconditioner )
        {
            hessian_data_type data( *stateView, *hessianView, *hessianView );
            solver->hessianFunction()( data );
            hessianView->closeIfNeeded();
        }
        else
        {
            auto preconditionerView = makeMatrixView(
                preconditioner, preconditionerStorage->mapRowPtr(),
                preconditionerStorage->mapColPtr(), callbackContext );
            hessian_data_type data(
                *stateView, *hessianView, *preconditionerView );
            solver->hessianFunction()( data );
            hessianView->closeIfNeeded();
            preconditionerView->closeIfNeeded();
        }
        return PETSC_SUCCESS;
    }
    catch ( ... )
    {
        if ( !callbackContext.exception )
            callbackContext.exception = std::current_exception();
        return PETSC_ERR_USER;
    }
}

template<typename T, typename SizeT>
PetscErrorCode
SolverOptimizationPetsc<T, SizeT>::monitorBridge( Tao tao, void* context ) noexcept
{
    auto& callbackContext = *static_cast<CallbackContext*>( context );
    try
    {
        auto* solver = callbackContext.solver;
        if ( !solver->hasMonitors() )
            return PETSC_SUCCESS;

        PetscInt iteration = 0;
        PetscReal objective = 0;
        PetscReal gradientNorm = 0;
        PetscReal constraintNorm = 0;
        PetscReal stepNorm = 0;
        TaoConvergedReason reason = TAO_CONTINUE_ITERATING;
        checkPetsc( TaoGetSolutionStatus( tao, &iteration, &objective, &gradientNorm,
                                          &constraintNorm, &stepNorm, &reason ),
                    "TaoGetSolutionStatus(monitor)" );

        monitor_record_type record;
        record.iteration = static_cast<size_type>( iteration );
        record.objective = static_cast<real_type>( objective );
        record.gradientNorm = static_cast<real_type>( gradientNorm );
        record.constraintNorm = static_cast<real_type>( constraintNorm );
        record.stepNorm = static_cast<real_type>( stepNorm );
        record.rawReason = static_cast<int>( reason );
        record.status = normalizeReason( reason );

        for ( auto const& registration : solver->monitors() )
            registration.second( record );
        return PETSC_SUCCESS;
    }
    catch ( ... )
    {
        if ( !callbackContext.exception )
            callbackContext.exception = std::current_exception();
        return PETSC_ERR_USER;
    }
}

template<typename T, typename SizeT>
std::string
SolverOptimizationPetsc<T, SizeT>::normalizePrefix( std::string prefix )
{
    if ( !prefix.empty() && prefix.back() != '_' )
        prefix.push_back( '_' );
    return prefix;
}

template<typename T, typename SizeT>
OptimizationStatus
SolverOptimizationPetsc<T, SizeT>::normalizeReason( TaoConvergedReason reason )
{
    switch ( reason )
    {
    case TAO_CONVERGED_GATOL:
    case TAO_CONVERGED_GRTOL:
    case TAO_CONVERGED_GTTOL:
    case TAO_CONVERGED_STEPTOL:
    case TAO_CONVERGED_MINF:
        return OptimizationStatus::Converged;
    case TAO_CONVERGED_USER:
    case TAO_DIVERGED_USER:
        return OptimizationStatus::UserStopped;
    case TAO_DIVERGED_MAXITS:
        return OptimizationStatus::MaximumIterations;
    case TAO_DIVERGED_MAXFCN:
        return OptimizationStatus::MaximumFunctionEvaluations;
    case TAO_DIVERGED_LS_FAILURE:
        return OptimizationStatus::DivergedLineSearch;
    case TAO_DIVERGED_TR_REDUCTION:
        return OptimizationStatus::DivergedTrustRegion;
    case TAO_DIVERGED_NAN:
        return OptimizationStatus::InvalidNumber;
    default:
        return OptimizationStatus::Unknown;
    }
}

template<typename T, typename SizeT>
std::string
SolverOptimizationPetsc<T, SizeT>::reasonString( TaoConvergedReason reason )
{
    switch ( reason )
    {
    case TAO_CONTINUE_ITERATING: return "continue iterating";
    case TAO_CONVERGED_GATOL: return "absolute gradient tolerance reached";
    case TAO_CONVERGED_GRTOL: return "relative gradient tolerance reached";
    case TAO_CONVERGED_GTTOL: return "initial-gradient tolerance reached";
    case TAO_CONVERGED_STEPTOL: return "step tolerance reached";
    case TAO_CONVERGED_MINF: return "objective lower bound reached";
    case TAO_CONVERGED_USER: return "user convergence test succeeded";
    case TAO_DIVERGED_MAXITS: return "maximum iterations reached";
    case TAO_DIVERGED_NAN: return "invalid numeric value encountered";
    case TAO_DIVERGED_MAXFCN: return "maximum function evaluations reached";
    case TAO_DIVERGED_LS_FAILURE: return "line search failed";
    case TAO_DIVERGED_TR_REDUCTION: return "trust-region reduction failed";
    case TAO_DIVERGED_USER: return "user convergence test failed";
    default: return "unknown TAO convergence reason";
    }
}

template<typename T, typename SizeT>
KSP
SolverOptimizationPetsc<T, SizeT>::nativeKsp() const
{
    if ( !M_tao )
        throw std::logic_error( "SolverOptimizationPetsc has been cleared" );

    KSP ksp = nullptr;
    checkPetsc( TaoGetKSP( M_tao, &ksp ), "TaoGetKSP" );
    return ksp;
}

template<typename T, typename SizeT>
void
SolverOptimizationPetsc<T, SizeT>::applyProgrammaticConfiguration()
{
    if ( !this->type().empty() )
        checkPetsc( TaoSetType( M_tao, this->type().c_str() ), "TaoSetType" );
    if ( auto const maximum = this->maximumIterations() )
        checkPetsc( TaoSetMaximumIterations( M_tao, static_cast<PetscInt>( *maximum ) ),
                    "TaoSetMaximumIterations" );
    if ( auto const maximum = this->maximumFunctionEvaluations() )
        checkPetsc( TaoSetMaximumFunctionEvaluations( M_tao, static_cast<PetscInt>( *maximum ) ),
                    "TaoSetMaximumFunctionEvaluations" );
    if ( auto const tolerances = this->gradientTolerances() )
        checkPetsc( TaoSetTolerances( M_tao, tolerances->absolute, tolerances->relative,
                                      tolerances->initial ),
                    "TaoSetTolerances" );
}

template<typename T, typename SizeT>
void
SolverOptimizationPetsc<T, SizeT>::configureCallbacks()
{
    if ( this->hasObjectiveGradient() )
    {
        checkPetsc( TaoSetObjective( M_tao, nullptr, nullptr ), "TaoSetObjective(clear)" );
        checkPetsc( TaoSetGradient( M_tao, nullptr, nullptr, nullptr ),
                    "TaoSetGradient(clear)" );
        checkPetsc( TaoSetObjectiveAndGradient(
                        M_tao, nullptr, &self_type::objectiveGradientBridge, &M_context ),
                    "TaoSetObjectiveAndGradient" );
        return;
    }

    checkPetsc( TaoSetObjectiveAndGradient( M_tao, nullptr, nullptr, nullptr ),
                "TaoSetObjectiveAndGradient(clear)" );
    checkPetsc( TaoSetObjective( M_tao, &self_type::objectiveBridge, &M_context ),
                "TaoSetObjective" );
    checkPetsc( TaoSetGradient(
                    M_tao, nullptr, &self_type::gradientBridge, &M_context ),
                "TaoSetGradient" );
}

template<typename T, typename SizeT>
void
SolverOptimizationPetsc<T, SizeT>::configureMonitors()
{
    checkPetsc( TaoMonitorCancel( M_tao ), "TaoMonitorCancel" );
    checkPetsc( TaoMonitorSet( M_tao, &self_type::monitorBridge, &M_context, nullptr ),
                "TaoMonitorSet" );
}

template<typename T, typename SizeT>
void
SolverOptimizationPetsc<T, SizeT>::configureHessian( vector_type const& variable )
{
    if ( !this->hasHessian() )
    {
        if ( M_hessianConfigured )
        {
            checkPetsc( TaoSetHessian( M_tao, nullptr, nullptr, nullptr, nullptr ),
                        "TaoSetHessian(clear)" );
            M_hessianConfigured = false;
        }
        return;
    }

    auto const& hessian = this->hessianMatrix();
    auto const& preconditioner = this->hessianPreconditionerMatrix();
    if ( !hessian || !preconditioner )
        throw std::invalid_argument(
            "SolverOptimizationPetsc requires Hessian and preconditioning matrices" );

    Mat const nativeHessian = this->nativeMatrix( *hessian, "Hessian matrix" );
    Mat const nativePreconditioner =
        this->nativeMatrix( *preconditioner, "Hessian preconditioning matrix" );

    auto validateLayout = [&]( sparse_matrix_type const& matrix, Mat native,
                               char const* role )
    {
        if ( !matrix.mapRowPtr() || !matrix.mapColPtr() )
            throw std::invalid_argument(
                std::string( "SolverOptimizationPetsc " ) + role +
                " requires row and column distribution maps" );
        if ( !haveSameLayout( variable.map(), *matrix.mapRowPtr() ) ||
             !haveSameLayout( variable.map(), *matrix.mapColPtr() ) )
            throw std::invalid_argument(
                std::string( "SolverOptimizationPetsc " ) + role +
                " must be square with the optimization variable's distributed layout" );

        PetscInt localRows = 0;
        PetscInt localColumns = 0;
        checkPetsc( MatGetLocalSize( native, &localRows, &localColumns ),
                    "MatGetLocalSize" );
        if ( localRows != static_cast<PetscInt>( variable.map().nLocalDofWithoutGhost() ) ||
             localColumns != static_cast<PetscInt>( variable.map().nLocalDofWithoutGhost() ) )
            throw std::invalid_argument(
                std::string( "SolverOptimizationPetsc " ) + role +
                " native local dimensions do not match the optimization variable" );
    };
    validateLayout( *hessian, nativeHessian, "Hessian matrix" );
    validateLayout( *preconditioner, nativePreconditioner,
                    "Hessian preconditioning matrix" );

    checkPetsc( TaoSetHessian( M_tao, nativeHessian, nativePreconditioner,
                               &self_type::hessianBridge, &M_context ),
                "TaoSetHessian" );
    M_hessianConfigured = true;
}

template<typename T, typename SizeT>
void
SolverOptimizationPetsc<T, SizeT>::configureVariableBounds( vector_type const& variable )
{
    if ( !this->hasVariableBounds() )
    {
        checkPetsc( TaoSetVariableBounds( M_tao, nullptr, nullptr ),
                    "TaoSetVariableBounds(clear)" );
        M_materializedLowerBound.reset();
        M_materializedUpperBound.reset();
        return;
    }

    Vec lower = nullptr;
    Vec upper = nullptr;
    bool materializedLower = false;
    bool materializedUpper = false;

    if ( this->lowerBound() )
    {
        this->validateBoundLayout( variable, *this->lowerBound(), "lower bound" );
        lower = this->nativeVector( *this->lowerBound(), "lower bound" );
    }
    else
    {
        real_type const value = this->constantLowerBound().value_or( super::negativeInfinity() );
        lower = this->materializeBound(
            variable, M_materializedLowerBound, value, "materialized lower bound" );
        materializedLower = true;
    }

    if ( this->upperBound() )
    {
        this->validateBoundLayout( variable, *this->upperBound(), "upper bound" );
        upper = this->nativeVector( *this->upperBound(), "upper bound" );
    }
    else
    {
        real_type const value = this->constantUpperBound().value_or( super::positiveInfinity() );
        upper = this->materializeBound(
            variable, M_materializedUpperBound, value, "materialized upper bound" );
        materializedUpper = true;
    }

    checkPetsc( TaoSetVariableBounds( M_tao, lower, upper ), "TaoSetVariableBounds" );
    if ( !materializedLower )
        M_materializedLowerBound.reset();
    if ( !materializedUpper )
        M_materializedUpperBound.reset();
}

template<typename T, typename SizeT>
void
SolverOptimizationPetsc<T, SizeT>::validateBoundLayout(
    vector_type const& variable, vector_type const& bound, char const* role ) const
{
    if ( !haveSameLayout( variable, bound ) )
        throw std::invalid_argument(
            std::string( "SolverOptimizationPetsc " ) + role +
            " must have the same distributed layout as the optimization variable" );
}

template<typename T, typename SizeT>
Vec
SolverOptimizationPetsc<T, SizeT>::materializeBound(
    vector_type const& variable, vector_ptrtype& storage, real_type value, char const* role )
{
    if ( !storage || !haveSameLayout( variable, *storage ) )
        storage = variable.clone();
    storage->setConstant( value );
    this->validateBoundLayout( variable, *storage, role );
    return this->nativeVector( *storage, role );
}

template<typename T, typename SizeT>
bool
SolverOptimizationPetsc<T, SizeT>::haveSameLayout(
    vector_type const& left, vector_type const& right )
{
    return haveSameLayout( left.map(), right.map() );
}

template<typename T, typename SizeT>
bool
SolverOptimizationPetsc<T, SizeT>::haveSameLayout(
    datamap_type const& left, datamap_type const& right )
{
    return left.nDof() == right.nDof() &&
           left.nLocalDofWithoutGhost() == right.nLocalDofWithoutGhost() &&
           left.nLocalDofWithGhost() == right.nLocalDofWithGhost() &&
           left.firstDofGlobalCluster() == right.firstDofGlobalCluster() &&
           left.lastDofGlobalCluster() == right.lastDofGlobalCluster() &&
           left.mapGlobalProcessToGlobalCluster() ==
               right.mapGlobalProcessToGlobalCluster();
}

template<typename T, typename SizeT>
void
SolverOptimizationPetsc<T, SizeT>::checkCommunicator(
    PetscObject object, char const* role ) const
{
    MPI_Comm objectComm = MPI_COMM_NULL;
    checkPetsc( PetscObjectGetComm( object, &objectComm ),
                "PetscObjectGetComm" );

    MPI_Comm solverComm = this->worldComm().globalComm();
    int comparison = MPI_UNEQUAL;
    int const error = MPI_Comm_compare( solverComm, objectComm, &comparison );
    if ( error != MPI_SUCCESS || ( comparison != MPI_IDENT && comparison != MPI_CONGRUENT ) )
        throw std::invalid_argument(
            std::string( "SolverOptimizationPetsc " ) + role +
            " uses an incompatible communicator" );
}

template<typename T, typename SizeT>
Vec
SolverOptimizationPetsc<T, SizeT>::nativeVector(
    vector_type const& vector, char const* role ) const
{
    auto const* petscVector = dynamic_cast<VectorPetsc<T> const*>( &vector );
    if ( !petscVector )
        throw std::invalid_argument(
            std::string( "SolverOptimizationPetsc " ) + role +
            " must be VectorPetsc-compatible" );
    if ( !petscVector->isInitialized() || !petscVector->vec() )
        throw std::invalid_argument(
            std::string( "SolverOptimizationPetsc " ) + role +
            " is an uninitialized PETSc vector" );

    this->checkCommunicator(
        reinterpret_cast<PetscObject>( petscVector->vec() ), role );
    return petscVector->vec();
}

template<typename T, typename SizeT>
Mat
SolverOptimizationPetsc<T, SizeT>::nativeMatrix(
    sparse_matrix_type const& matrix, char const* role ) const
{
    auto const* petscMatrix = dynamic_cast<MatrixPetsc<T> const*>( &matrix );
    if ( !petscMatrix )
        throw std::invalid_argument(
            std::string( "SolverOptimizationPetsc " ) + role +
            " must be MatrixPetsc-compatible" );
    if ( !petscMatrix->isInitialized() )
        throw std::invalid_argument(
            std::string( "SolverOptimizationPetsc " ) + role +
            " is an uninitialized PETSc matrix" );

    Mat const native = petscMatrix->mat();
    this->checkCommunicator( reinterpret_cast<PetscObject>( native ), role );
    return native;
}

template<typename T, typename SizeT>
void
SolverOptimizationPetsc<T, SizeT>::checkPetsc(
    PetscErrorCode error, char const* operation )
{
    if ( !error )
        return;

    char const* message = nullptr;
    PetscErrorMessage( error, &message, nullptr );
    std::ostringstream description;
    description << operation << " failed with PETSc error " << error;
    if ( message )
        description << ": " << message;
    throw std::runtime_error( description.str() );
}

template class SolverOptimizationPetsc<double>;

} // namespace Feel

#endif // FEELPP_HAS_PETSC_TAO
