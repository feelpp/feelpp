/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#ifndef FEELPP_FEELOPT_SOLVEROPTIMIZATION_HPP
#define FEELPP_FEELOPT_SOLVEROPTIMIZATION_HPP 1

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <feel/feelcore/commobject.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelopt/solveroptimizationoptions.hpp>

namespace Feel
{

template<typename T, typename SizeT> class Vector;
template<typename T> class MatrixSparse;

/** Normalized backend-independent optimization termination status. */
enum class OptimizationStatus
{
    Converged,                  ///< A solver convergence criterion was satisfied.
    UserStopped,                ///< A user convergence test stopped the solve.
    MaximumIterations,          ///< The iteration limit was reached.
    MaximumFunctionEvaluations, ///< The function-evaluation limit was reached.
    DivergedLineSearch,         ///< The globalization line search failed.
    DivergedTrustRegion,        ///< The trust-region update failed.
    DivergedLinearSolve,        ///< An inner linear solve failed.
    InvalidNumber,              ///< A NaN or another invalid numeric value was found.
    CallbackFailure,            ///< A user callback failed or threw an exception.
    InvalidConfiguration,       ///< Required callbacks or solver data are missing.
    Unsupported,                ///< The requested operation is unsupported by the backend.
    Unknown                     ///< The backend reason has no normalized mapping.
};

/** Structured information delivered to optimization monitor callbacks. */
template<typename T, typename SizeT = uint32_type>
struct OptimizationMonitorRecord
{
    using value_type = T;
    using real_type = typename type_traits<value_type>::real_type;
    using size_type = SizeT;

    size_type iteration = 0; ///< Current optimization iteration.
    real_type objective = std::numeric_limits<real_type>::quiet_NaN(); ///< Current objective value.
    real_type gradientNorm = std::numeric_limits<real_type>::quiet_NaN(); ///< Current gradient norm.
    real_type constraintNorm = std::numeric_limits<real_type>::quiet_NaN(); ///< Current constraint norm.
    real_type stepNorm = std::numeric_limits<real_type>::quiet_NaN(); ///< Current step norm.
    int rawReason = 0; ///< Backend termination code at this iteration.
    OptimizationStatus status = OptimizationStatus::Unknown; ///< Normalized iteration status.
};

/**
 * Algebraic data supplied to a combined objective-and-gradient callback.
 *
 * All references are borrowed for one callback invocation. The callback fills
 * @ref gradient and assigns the computed scalar value to @ref objective.
 *
 * @tparam T optimization scalar type
 * @tparam SizeT algebraic index type
 */
template<typename T, typename SizeT = uint32_type>
struct OptimizationObjectiveGradientData
{
    using value_type = T; ///< Optimization scalar type.
    using real_type = typename type_traits<value_type>::real_type; ///< Objective scalar type.
    using vector_type = Vector<value_type, SizeT>; ///< Distributed algebraic vector type.

    /**
     * Bind one callback invocation to its distributed algebraic storage.
     *
     * @param currentState current optimization variable
     * @param gradientStorage writable gradient vector supplied by the backend
     */
    OptimizationObjectiveGradientData(
        vector_type const& currentState, vector_type& gradientStorage )
        : state( currentState ), gradient( gradientStorage )
    {}

    vector_type const& state; ///< Current distributed optimization variable.
    vector_type& gradient; ///< Writable distributed gradient storage.
    real_type objective = std::numeric_limits<real_type>::quiet_NaN(); ///< Computed objective value.
};

/**
 * Algebraic data supplied to an exact Hessian callback.
 *
 * All references are borrowed for one callback invocation. The Hessian and
 * preconditioner may alias when one assembled matrix serves both roles.
 *
 * @tparam T optimization scalar type
 * @tparam SizeT algebraic index type
 */
template<typename T, typename SizeT = uint32_type>
struct OptimizationHessianData
{
    using value_type = T; ///< Optimization scalar type.
    using vector_type = Vector<value_type, SizeT>; ///< Distributed algebraic vector type.
    using matrix_type = MatrixSparse<value_type>; ///< Distributed sparse-matrix type.

    /**
     * Bind one callback invocation to its distributed algebraic storage.
     *
     * @param currentState current optimization variable
     * @param hessianStorage writable exact Hessian matrix
     * @param preconditionerStorage writable Hessian preconditioning matrix
     */
    OptimizationHessianData(
        vector_type const& currentState, matrix_type& hessianStorage,
        matrix_type& preconditionerStorage )
        : state( currentState ), hessian( hessianStorage ),
          preconditioner( preconditionerStorage )
    {}

    /** @return true when one matrix is used for both Hessian roles */
    [[nodiscard]] bool sharesMatrix() const noexcept
    {
        return std::addressof( hessian ) == std::addressof( preconditioner );
    }

    vector_type const& state; ///< Current distributed optimization variable.
    matrix_type& hessian; ///< Writable distributed exact Hessian.
    matrix_type& preconditioner; ///< Writable distributed Hessian preconditioner.
};

/** Structured result returned by an optimization solve. */
template<typename T, typename SizeT = uint32_type>
struct OptimizationResult
{
    using value_type = T;
    using real_type = typename type_traits<value_type>::real_type;
    using size_type = SizeT;

    bool converged = false; ///< True only when the backend reports a positive convergence reason.
    OptimizationStatus status = OptimizationStatus::Unknown; ///< Normalized termination status.
    int rawReason = 0; ///< Backend termination code represented without backend-specific types.
    std::string reason; ///< Human-readable termination reason.
    size_type iterations = 0; ///< Number of optimization iterations performed.
    size_type functionEvaluations = 0; ///< Number of objective evaluations performed.
    real_type objective = std::numeric_limits<real_type>::quiet_NaN(); ///< Final objective value.
    real_type gradientNorm = std::numeric_limits<real_type>::quiet_NaN(); ///< Final gradient norm.
    real_type constraintNorm = std::numeric_limits<real_type>::quiet_NaN(); ///< Final constraint norm.
    real_type stepNorm = std::numeric_limits<real_type>::quiet_NaN(); ///< Norm of the final step.
    std::string solverType; ///< Effective backend solver type.
    std::string optionsPrefix; ///< Effective runtime-options prefix.
    std::string diagnostic; ///< Additional diagnostic information for unsuccessful solves.
};

/**
 * Backend-independent optimization solver contract.
 *
 * Callback vector references are borrowed and valid only for the duration of
 * the callback. Implementations operate directly on distributed Feel++
 * vectors and must not retain those references.
 */
template<typename T, typename SizeT = uint32_type>
class SolverOptimization : public CommObject
{
public:
    using super = CommObject;
    using self_type = SolverOptimization<T, SizeT>;
    using ptrtype = std::shared_ptr<self_type>; ///< Shared ownership of a generic solver instance.
    using value_type = T;
    using real_type = typename type_traits<value_type>::real_type;
    using size_type = SizeT;
    using vector_type = Vector<value_type, size_type>;
    using vector_ptrtype = std::shared_ptr<vector_type>;
    using sparse_matrix_type = MatrixSparse<value_type>;
    using sparse_matrix_ptrtype = std::shared_ptr<sparse_matrix_type>;
    using result_type = OptimizationResult<value_type, size_type>;
    using monitor_record_type = OptimizationMonitorRecord<value_type, size_type>;
    using objective_gradient_data_type =
        OptimizationObjectiveGradientData<value_type, size_type>;
    using hessian_data_type = OptimizationHessianData<value_type, size_type>;
    using monitor_function_type = std::function<void( monitor_record_type const& )>;
    using monitor_id_type = std::size_t;
    using objective_function_type = std::function<real_type( vector_type const& )>;
    using gradient_function_type = std::function<void( vector_type const&, vector_type& )>;
    using objective_gradient_function_type =
        std::function<void( objective_gradient_data_type& )>;
    using hessian_function_type = std::function<void( hessian_data_type& )>;

    /**
     * Construct an optimization solver contract.
     *
     * @param prefix independent runtime-options prefix
     * @param worldComm communicator used collectively by the solver
     * @param vm Feel++ program-options values used to configure this instance
     */
    explicit SolverOptimization(
        std::string prefix = {},
        worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
        po::variables_map const& vm = Environment::vm() )
        :
        super( worldComm ),
        M_prefix( std::move( prefix ) )
    {
        Environment::addOptions( solveroptimization_options( M_prefix ) );
        this->configureFromFeelppOptions( vm );
    }

    /** Destroy the generic solver state. */
    ~SolverOptimization() override = default;

    /** Copy the generic callback and configuration state. */
    SolverOptimization( SolverOptimization const& ) = default;

    /** Move the generic callback and configuration state. */
    SolverOptimization( SolverOptimization&& ) = default;

    /** Copy-assign the generic callback and configuration state. */
    SolverOptimization& operator=( SolverOptimization const& ) = default;

    /** Move-assign the generic callback and configuration state. */
    SolverOptimization& operator=( SolverOptimization&& ) = default;

    /**
     * Build an optimization solver through the backend-independent contract.
     *
     * This follows the existing algebra solver factory convention while
     * keeping backend-native types out of normal user code.
     *
     * @param backend backend identifier, currently `petsc`
     * @param prefix independent runtime-options prefix
     * @param worldComm communicator used collectively by the solver
     * @param vm Feel++ program-options values
     * @return shared generic optimization solver
     * @throws std::invalid_argument when @p backend is unknown
     * @throws std::runtime_error when the requested backend is unavailable
     */
    [[nodiscard]] static ptrtype build(
        std::string const& backend = "petsc", std::string const& prefix = {},
        worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
        po::variables_map const& vm = Environment::vm() );

    /**
     * Register an objective-only callback.
     *
     * This callback is used with the separate gradient callback only when no
     * combined objective-and-gradient callback is registered.
     *
     * @param callback copy-constructible callable compatible with
     *        `real_type(Vector const&)`
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, vector_type const&> &&
                 std::convertible_to<
                     std::invoke_result_t<std::decay_t<Callable>&, vector_type const&>,
                     real_type>
    void setObjective( Callable&& callback )
    {
        M_objective =
            [fn = std::forward<Callable>( callback )]( vector_type const& x ) mutable
            {
                return static_cast<real_type>( std::invoke( fn, x ) );
            };
    }

    /** Remove the registered objective-only callback. */
    void resetObjective() { M_objective = {}; }

    /** @return true when an objective-only callback is registered */
    [[nodiscard]] bool hasObjective() const noexcept { return bool( M_objective ); }

    /**
     * Register a gradient-only callback.
     *
     * This callback is used with the separate objective callback only when no
     * combined objective-and-gradient callback is registered.
     *
     * @param callback copy-constructible callable compatible with
     *        `void(Vector const&, Vector&)`
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, vector_type const&, vector_type&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&, vector_type const&, vector_type&>,
                     void>
    void setGradient( Callable&& callback )
    {
        M_gradient =
            [fn = std::forward<Callable>( callback )]( vector_type const& x, vector_type& gradient ) mutable
            {
                std::invoke( fn, x, gradient );
            };
    }

    /** Remove the registered gradient-only callback. */
    void resetGradient() { M_gradient = {}; }

    /** @return true when a gradient-only callback is registered */
    [[nodiscard]] bool hasGradient() const noexcept { return bool( M_gradient ); }

    /**
     * Register a structured combined objective-and-gradient callback.
     *
     * @param callback copy-constructible callable compatible with
     *        `void(objective_gradient_data_type&)`
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&,
                                objective_gradient_data_type&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&,
                                          objective_gradient_data_type&>, void>
    void setObjectiveGradient( Callable&& callback )
    {
        M_objectiveGradient =
            [fn = std::forward<Callable>( callback )](
                objective_gradient_data_type& data ) mutable
            {
                std::invoke( fn, data );
            };
    }

    /**
     * Register a combined objective-and-gradient callback.
     *
     * The vector references are borrowed for the duration of the invocation.
     * The callback must fill @p gradient and return the objective value.
     *
     * @param callback copy-constructible callable compatible with
     *        `real_type(Vector const&, Vector&)`
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, vector_type const&, vector_type&> &&
                 std::convertible_to<
                     std::invoke_result_t<std::decay_t<Callable>&, vector_type const&, vector_type&>,
                     real_type>
    void setObjectiveGradient( Callable&& callback )
    {
        M_objectiveGradient =
            [fn = std::forward<Callable>( callback )](
                objective_gradient_data_type& data ) mutable
            {
                data.objective = static_cast<real_type>(
                    std::invoke( fn, data.state, data.gradient ) );
            };
    }

    /** Remove the registered combined objective-and-gradient callback. */
    void resetObjectiveGradient() { M_objectiveGradient = {}; }

    /** @return true when a combined objective-and-gradient callback is registered */
    [[nodiscard]] bool hasObjectiveGradient() const noexcept { return bool( M_objectiveGradient ); }

    /**
     * Register a structured exact Hessian callback with distinct matrices.
     *
     * @param hessian matrix filled with the exact Hessian
     * @param preconditioner matrix filled for Hessian preconditioning
     * @param callback copy-constructible callable compatible with
     *        `void(hessian_data_type&)`
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, hessian_data_type&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&,
                                          hessian_data_type&>, void>
    void setHessian( sparse_matrix_ptrtype hessian,
                     sparse_matrix_ptrtype preconditioner,
                     Callable&& callback )
    {
        if ( !hessian )
            throw std::invalid_argument(
                "SolverOptimization::setHessian received a null Hessian matrix" );
        if ( !preconditioner )
            throw std::invalid_argument(
                "SolverOptimization::setHessian received a null preconditioning matrix" );

        M_hessian = std::move( hessian );
        M_hessianPreconditioner = std::move( preconditioner );
        M_hessianFunction =
            [fn = std::forward<Callable>( callback )]( hessian_data_type& data ) mutable
            {
                std::invoke( fn, data );
            };
    }

    /**
     * Register an exact Hessian callback with distinct operator and
     * preconditioning matrices.
     *
     * The matrices are allocated outside the callback, retained by shared
     * ownership, and reused by the backend between iterations. Callback
     * references are borrowed and must not be retained. This assembled-matrix
     * contract leaves Hessian-vector and matrix-free interfaces available for
     * future, separate extensions.
     *
     * @param hessian matrix filled with the exact Hessian
     * @param preconditioner matrix filled for Hessian preconditioning
     * @param callback copy-constructible callable compatible with
     *        `void(Vector const&, MatrixSparse&, MatrixSparse&)`
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, vector_type const&,
                                sparse_matrix_type&, sparse_matrix_type&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&, vector_type const&,
                                          sparse_matrix_type&, sparse_matrix_type&>,
                     void>
    void setHessian( sparse_matrix_ptrtype hessian,
                     sparse_matrix_ptrtype preconditioner,
                     Callable&& callback )
    {
        if ( !hessian )
            throw std::invalid_argument(
                "SolverOptimization::setHessian received a null Hessian matrix" );
        if ( !preconditioner )
            throw std::invalid_argument(
                "SolverOptimization::setHessian received a null preconditioning matrix" );

        hessian_function_type function =
            [fn = std::forward<Callable>( callback )]( hessian_data_type& data ) mutable
            {
                std::invoke( fn, data.state, data.hessian, data.preconditioner );
            };
        M_hessian = std::move( hessian );
        M_hessianPreconditioner = std::move( preconditioner );
        M_hessianFunction = std::move( function );
    }

    /**
     * Register a structured exact Hessian callback using one shared matrix.
     *
     * @param hessian matrix filled and reused for both Hessian roles
     * @param callback copy-constructible callable compatible with
     *        `void(hessian_data_type&)`
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, hessian_data_type&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&,
                                          hessian_data_type&>, void>
    void setHessian( sparse_matrix_ptrtype hessian, Callable&& callback )
    {
        auto preconditioner = hessian;
        this->setHessian( std::move( hessian ), std::move( preconditioner ),
                          std::forward<Callable>( callback ) );
    }

    /**
     * Register an exact Hessian callback using one matrix for both the Hessian
     * operator and its preconditioner.
     *
     * @param hessian matrix filled and reused for both Hessian roles
     * @param callback copy-constructible callable compatible with
     *        `void(Vector const&, MatrixSparse&, MatrixSparse&)`
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, vector_type const&,
                                sparse_matrix_type&, sparse_matrix_type&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&, vector_type const&,
                                          sparse_matrix_type&, sparse_matrix_type&>,
                     void>
    void setHessian( sparse_matrix_ptrtype hessian, Callable&& callback )
    {
        auto preconditioner = hessian;
        this->setHessian( std::move( hessian ), std::move( preconditioner ),
                          std::forward<Callable>( callback ) );
    }

    /** Remove the exact Hessian callback and release its registered matrices. */
    void resetHessian()
    {
        M_hessianFunction = {};
        M_hessian.reset();
        M_hessianPreconditioner.reset();
    }

    /** @return true when an exact Hessian callback and matrices are registered */
    [[nodiscard]] bool hasHessian() const noexcept { return bool( M_hessianFunction ); }

    /**
     * Add an instance-local C++ iteration monitor.
     *
     * PETSc-native monitors configured through runtime options remain
     * independent. The returned identifier can be passed to removeMonitor().
     * Monitor registration must not be changed from inside a running monitor.
     *
     * @param monitor copy-constructible callable compatible with
     *        `void(OptimizationMonitorRecord const&)`
     * @return stable identifier for this monitor registration
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, monitor_record_type const&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&,
                                          monitor_record_type const&>,
                     void>
    [[nodiscard]] monitor_id_type addMonitor( Callable&& monitor )
    {
        monitor_function_type function =
            [fn = std::forward<Callable>( monitor )]( monitor_record_type const& record ) mutable
            {
                std::invoke( fn, record );
            };
        monitor_id_type const identifier = M_nextMonitorId++;
        M_monitors.emplace_back( identifier, std::move( function ) );
        return identifier;
    }

    /**
     * Remove one C++ monitor registration.
     *
     * @param identifier value returned by addMonitor()
     * @return true when a monitor was removed
     */
    bool removeMonitor( monitor_id_type identifier )
    {
        return std::erase_if(
                   M_monitors,
                   [identifier]( auto const& registration )
                   {
                       return registration.first == identifier;
                   } ) > 0;
    }

    /** Remove all registered C++ iteration monitors. */
    void resetMonitors() noexcept { M_monitors.clear(); }

    /** @return true when at least one C++ monitor is registered */
    [[nodiscard]] bool hasMonitors() const noexcept { return !M_monitors.empty(); }

    /** @return number of registered C++ iteration monitors */
    [[nodiscard]] std::size_t monitorCount() const noexcept { return M_monitors.size(); }

    /** Set the programmatic solver type applied before runtime options. */
    void setType( std::string type ) { M_type = std::move( type ); }

    /** @return the programmatic solver type */
    [[nodiscard]] std::string const& type() const noexcept { return M_type; }

    /** Set the maximum number of optimization iterations. */
    void setMaximumIterations( size_type value ) { M_maxIterations = value; }

    /** Set the maximum number of objective evaluations. */
    void setMaximumFunctionEvaluations( size_type value ) { M_maxFunctionEvaluations = value; }

    /**
     * Set absolute, relative, and initial-gradient stopping tolerances.
     *
     * @param absolute absolute gradient tolerance
     * @param relative objective-relative gradient tolerance
     * @param initial initial-gradient-relative tolerance
     */
    void setGradientTolerances( real_type absolute, real_type relative, real_type initial )
    {
        M_gradientTolerances = GradientTolerances{ absolute, relative, initial };
    }

    /**
     * Set the step or trust-region-radius tolerance.
     *
     * Backend runtime options retain precedence over this programmatic value.
     *
     * @param value nonnegative step tolerance
     * @throws std::invalid_argument when @p value is negative
     */
    void setStepTolerance( real_type value )
    {
        if ( !( value >= real_type( 0 ) ) )
            throw std::invalid_argument(
                "SolverOptimization::setStepTolerance requires a nonnegative value" );
        M_stepTolerance = value;
    }

    /** Remove the optional programmatic step tolerance. */
    void resetStepTolerance() noexcept { M_stepTolerance.reset(); }

    /** Enable or disable the backend's standard text iteration monitor. */
    void setMonitorEnabled( bool enabled = true ) noexcept { M_monitorEnabled = enabled; }

    /** @return true when the standard backend iteration monitor is enabled */
    [[nodiscard]] bool monitorEnabled() const noexcept { return M_monitorEnabled; }

    /** Enable or disable printing the backend convergence reason. */
    void setConvergedReasonEnabled( bool enabled = true ) noexcept
    {
        M_convergedReasonEnabled = enabled;
    }

    /** @return true when backend convergence-reason output is enabled */
    [[nodiscard]] bool convergedReasonEnabled() const noexcept
    {
        return M_convergedReasonEnabled;
    }

    /** Enable or disable printing the configured backend solver. */
    void setViewEnabled( bool enabled = true ) noexcept { M_viewEnabled = enabled; }

    /** @return true when backend solver-view output is enabled */
    [[nodiscard]] bool viewEnabled() const noexcept { return M_viewEnabled; }

    /**
     * Set vector-valued lower and upper variable bounds.
     *
     * Either pointer may be null to request a one-sided bound. Calling this
     * function replaces previously configured constant bounds.
     *
     * @param lower lower-bound vector or null for no explicit lower bound
     * @param upper upper-bound vector or null for no explicit upper bound
     */
    void setVariableBounds( vector_ptrtype lower, vector_ptrtype upper )
    {
        M_lowerBound = std::move( lower );
        M_upperBound = std::move( upper );
        M_constantLowerBound.reset();
        M_constantUpperBound.reset();
    }

    /** Set or clear the vector-valued lower bound. */
    void setLowerBound( vector_ptrtype lower )
    {
        M_lowerBound = std::move( lower );
        M_constantLowerBound.reset();
    }

    /** Set or clear the vector-valued upper bound. */
    void setUpperBound( vector_ptrtype upper )
    {
        M_upperBound = std::move( upper );
        M_constantUpperBound.reset();
    }

    /** Set a constant lower bound materialized with the solution layout at solve time. */
    void setConstantLowerBound( real_type value )
    {
        M_lowerBound.reset();
        M_constantLowerBound = value;
    }

    /** Set a constant upper bound materialized with the solution layout at solve time. */
    void setConstantUpperBound( real_type value )
    {
        M_upperBound.reset();
        M_constantUpperBound = value;
    }

    /** Remove both vector-valued and constant variable bounds. */
    void resetVariableBounds()
    {
        M_lowerBound.reset();
        M_upperBound.reset();
        M_constantLowerBound.reset();
        M_constantUpperBound.reset();
    }

    /** @return true when a lower variable bound is configured */
    [[nodiscard]] bool hasLowerBound() const noexcept
    {
        return bool( M_lowerBound ) || M_constantLowerBound.has_value();
    }

    /** @return true when an upper variable bound is configured */
    [[nodiscard]] bool hasUpperBound() const noexcept
    {
        return bool( M_upperBound ) || M_constantUpperBound.has_value();
    }

    /** @return true when at least one variable bound is configured */
    [[nodiscard]] bool hasVariableBounds() const noexcept
    {
        return this->hasLowerBound() || this->hasUpperBound();
    }

    /**
     * Return the finite sentinel representing an unbounded lower component.
     *
     * @return backend-independent negative-infinity bound value
     */
    [[nodiscard]] static constexpr real_type negativeInfinity() noexcept
    {
        return -std::numeric_limits<real_type>::max() / real_type( 4 );
    }

    /**
     * Return the finite sentinel representing an unbounded upper component.
     *
     * @return backend-independent positive-infinity bound value
     */
    [[nodiscard]] static constexpr real_type positiveInfinity() noexcept
    {
        return std::numeric_limits<real_type>::max() / real_type( 4 );
    }

    /** @return the backend-independent prefix supplied at construction */
    [[nodiscard]] std::string const& prefix() const noexcept { return M_prefix; }

    /** Apply runtime options, which override programmatic defaults. */
    virtual void setFromOptions() = 0;

    /** Release backend resources; repeated calls must be safe. */
    virtual void clear() noexcept = 0;

    /**
     * Solve in place from the initial value stored in @p variable.
     *
     * @return structured convergence and diagnostic information
     */
    [[nodiscard]] virtual result_type solve( vector_type& variable ) = 0;

    /** Shared-pointer convenience overload for solve(). */
    [[nodiscard]] result_type solve( vector_ptrtype const& variable )
    {
        if ( !variable )
            throw std::invalid_argument( "SolverOptimization::solve received a null vector" );
        return this->solve( *variable );
    }

protected:
    /** Internal pairing of a monitor identifier and its type-erased callable. */
    using monitor_registration_type = std::pair<monitor_id_type, monitor_function_type>;

    /** Programmatic gradient stopping tolerances. */
    struct GradientTolerances
    {
        real_type absolute; ///< Absolute gradient tolerance.
        real_type relative; ///< Objective-relative gradient tolerance.
        real_type initial; ///< Initial-gradient-relative tolerance.
    };

    /** @return the registered objective-only callback for backend implementations */
    [[nodiscard]] objective_function_type const& objective() const noexcept
    {
        return M_objective;
    }

    /** @return the registered gradient-only callback for backend implementations */
    [[nodiscard]] gradient_function_type const& gradient() const noexcept
    {
        return M_gradient;
    }

    /** @return the registered combined callback for backend implementations */
    [[nodiscard]] objective_gradient_function_type const& objectiveGradient() const noexcept
    {
        return M_objectiveGradient;
    }

    /** @return the registered exact Hessian callback for backend implementations */
    [[nodiscard]] hessian_function_type const& hessianFunction() const noexcept
    {
        return M_hessianFunction;
    }

    /** @return registered C++ monitors for backend callback bridges */
    [[nodiscard]] std::vector<monitor_registration_type> const& monitors() const noexcept
    {
        return M_monitors;
    }

    /** @return the registered exact Hessian matrix */
    [[nodiscard]] sparse_matrix_ptrtype const& hessianMatrix() const noexcept
    {
        return M_hessian;
    }

    /** @return the registered Hessian preconditioning matrix */
    [[nodiscard]] sparse_matrix_ptrtype const& hessianPreconditionerMatrix() const noexcept
    {
        return M_hessianPreconditioner;
    }

    /** @return the optional programmatic iteration limit */
    [[nodiscard]] std::optional<size_type> maximumIterations() const noexcept
    {
        return M_maxIterations;
    }

    /** @return the optional programmatic function-evaluation limit */
    [[nodiscard]] std::optional<size_type> maximumFunctionEvaluations() const noexcept
    {
        return M_maxFunctionEvaluations;
    }

    /** @return the optional programmatic gradient tolerances */
    [[nodiscard]] std::optional<GradientTolerances> gradientTolerances() const noexcept
    {
        return M_gradientTolerances;
    }

    /** @return the optional programmatic step tolerance */
    [[nodiscard]] std::optional<real_type> stepTolerance() const noexcept
    {
        return M_stepTolerance;
    }

    /** @return the optional vector-valued lower bound */
    [[nodiscard]] vector_ptrtype const& lowerBound() const noexcept { return M_lowerBound; }

    /** @return the optional vector-valued upper bound */
    [[nodiscard]] vector_ptrtype const& upperBound() const noexcept { return M_upperBound; }

    /** @return the optional constant lower bound */
    [[nodiscard]] std::optional<real_type> constantLowerBound() const noexcept
    {
        return M_constantLowerBound;
    }

    /** @return the optional constant upper bound */
    [[nodiscard]] std::optional<real_type> constantUpperBound() const noexcept
    {
        return M_constantUpperBound;
    }

private:
    /** Initialize generic configuration from prefixed Feel++ options. */
    void configureFromFeelppOptions( po::variables_map const& vm )
    {
        auto const key = [this]( std::string const& name )
        {
            return prefixvm( M_prefix, name );
        };

        if ( auto const it = vm.find( key( "tao-type" ) ); it != vm.end() )
            M_type = it->second.template as<std::string>();
        if ( auto const it = vm.find( key( "tao-maxit" ) ); it != vm.end() )
        {
            auto const value = it->second.template as<int>();
            if ( value <= 0 )
                throw std::invalid_argument(
                    "SolverOptimization option tao-maxit must be positive" );
            M_maxIterations = static_cast<size_type>( value );
        }
        if ( auto const it = vm.find( key( "tao-maxfcn" ) ); it != vm.end() )
        {
            auto const value = it->second.template as<int>();
            if ( value <= 0 )
                throw std::invalid_argument(
                    "SolverOptimization option tao-maxfcn must be positive" );
            M_maxFunctionEvaluations = static_cast<size_type>( value );
        }

        auto const absolute = vm.find( key( "tao-gatol" ) );
        auto const relative = vm.find( key( "tao-grtol" ) );
        auto const initial = vm.find( key( "tao-gttol" ) );
        bool const anyGradientTolerance =
            absolute != vm.end() || relative != vm.end() || initial != vm.end();
        bool const allGradientTolerances =
            absolute != vm.end() && relative != vm.end() && initial != vm.end();
        if ( anyGradientTolerance && !allGradientTolerances )
            throw std::invalid_argument(
                "SolverOptimization requires tao-gatol, tao-grtol, and tao-gttol together" );
        if ( allGradientTolerances )
            this->setGradientTolerances(
                absolute->second.template as<real_type>(),
                relative->second.template as<real_type>(),
                initial->second.template as<real_type>() );

        if ( auto const it = vm.find( key( "tao-steptol" ) ); it != vm.end() )
            this->setStepTolerance( it->second.template as<real_type>() );
        if ( auto const it = vm.find( key( "tao-monitor" ) ); it != vm.end() )
            M_monitorEnabled = it->second.template as<bool>();
        if ( auto const it = vm.find( key( "tao-converged-reason" ) ); it != vm.end() )
            M_convergedReasonEnabled = it->second.template as<bool>();
        if ( auto const it = vm.find( key( "tao-view" ) ); it != vm.end() )
            M_viewEnabled = it->second.template as<bool>();
    }

    std::string M_prefix; ///< Backend-independent runtime-options prefix.
    std::string M_type = "lmvm"; ///< Programmatic solver type.
    objective_function_type M_objective; ///< Registered objective-only callback.
    gradient_function_type M_gradient; ///< Registered gradient-only callback.
    objective_gradient_function_type M_objectiveGradient; ///< Registered combined callback.
    hessian_function_type M_hessianFunction; ///< Registered exact Hessian callback.
    sparse_matrix_ptrtype M_hessian; ///< User-allocated exact Hessian matrix.
    sparse_matrix_ptrtype M_hessianPreconditioner; ///< User-allocated Hessian preconditioner.
    std::vector<monitor_registration_type> M_monitors; ///< Instance-local C++ monitor callbacks.
    monitor_id_type M_nextMonitorId = 0; ///< Next unique monitor registration identifier.
    std::optional<size_type> M_maxIterations; ///< Optional iteration limit.
    std::optional<size_type> M_maxFunctionEvaluations; ///< Optional evaluation limit.
    std::optional<GradientTolerances> M_gradientTolerances; ///< Optional stopping tolerances.
    std::optional<real_type> M_stepTolerance; ///< Optional step or trust-radius tolerance.
    bool M_monitorEnabled = false; ///< Standard backend iteration-monitor state.
    bool M_convergedReasonEnabled = false; ///< Backend convergence-reason output state.
    bool M_viewEnabled = false; ///< Backend solver-view output state.
    vector_ptrtype M_lowerBound; ///< Optional vector-valued lower bound.
    vector_ptrtype M_upperBound; ///< Optional vector-valued upper bound.
    std::optional<real_type> M_constantLowerBound; ///< Optional constant lower bound.
    std::optional<real_type> M_constantUpperBound; ///< Optional constant upper bound.
};

/**
 * Build a backend-independent optimization solver with Feel++ named arguments.
 *
 * Supported named arguments are:
 *
 * - `_backend`: backend identifier, defaulting to `petsc`;
 * - `_name`: independent solver and runtime-options prefix;
 * - `_worldcomm`: communicator used collectively by the solver.
 * - `_vm`: Feel++ program-options values.
 *
 * @tparam T optimization scalar type
 * @tparam SizeT algebra index type
 * @tparam Ts named-argument types
 * @param arguments Feel++ named arguments
 * @return shared generic optimization solver
 */
template<typename T = double, typename SizeT = uint32_type, typename... Ts>
[[nodiscard]] typename SolverOptimization<T, SizeT>::ptrtype
optimizationSolver( Ts&&... arguments )
{
    auto args = NA::make_arguments( std::forward<Ts>( arguments )... );
    std::string const backend = args.get_else( _backend, std::string( "petsc" ) );
    std::string const name = args.get_else( _name, std::string{} );
    worldcomm_ptr_t const worldComm =
        args.get_else( _worldcomm, Environment::worldCommPtr() );
    po::variables_map const& vm = args.get_else( _vm, Environment::vm() );
    return SolverOptimization<T, SizeT>::build( backend, name, worldComm, vm );
}

} // namespace Feel

#endif // FEELPP_FEELOPT_SOLVEROPTIMIZATION_HPP
