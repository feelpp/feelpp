/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#ifndef FEELPP_FEELOPT_SOLVEROPTIMIZATIONFIELD_HPP
#define FEELPP_FEELOPT_SOLVEROPTIMIZATIONFIELD_HPP 1

#include <concepts>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <feel/feelalg/backend.hpp>
#include <feel/feelopt/solveroptimization.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/form.hpp>

namespace Feel
{

namespace detail
{

/** True when a space uses the static or dynamic Feel++ product-space API. */
template<typename SpaceType>
inline constexpr bool isOptimizationProductSpace =
    StaticProductSpacesType<SpaceType> || DynamicProductSpaceType<SpaceType>;

/**
 * Return the communicator owned by a scalar or recursively nested product space.
 *
 * @tparam SpacePtrType shared function-space or product-space pointer type
 * @param space space whose original Feel++ communicator is requested
 * @return communicator shared by the product-space entries
 */
template<typename SpacePtrType>
[[nodiscard]] worldcomm_ptr_t
optimizationSpaceWorldComm( SpacePtrType const& space )
{
    using space_type = typename std::remove_cvref_t<SpacePtrType>::element_type;
    if constexpr ( StaticProductSpacesType<space_type> )
        return optimizationSpaceWorldComm( hana::front( space->tupleSpaces() ) );
    else if constexpr ( DynamicProductSpaceType<space_type> )
        return optimizationSpaceWorldComm( ( *space )[0] );
    else
        return space->worldCommPtr();
}

} // namespace detail

/**
 * Writable variational linear form backed by an optimization gradient vector.
 *
 * This non-owning proxy lets a field callback assemble a gradient with the
 * natural Feel++ notation `gradient = integrate(...)`. The borrowed algebraic
 * vector remains valid only during the callback that receives this proxy.
 *
 * @tparam SpaceType Feel++ function-space type
 */
template<typename SpaceType,
         bool IsProductSpace = Feel::detail::isOptimizationProductSpace<SpaceType>>
class OptimizationGradientForm
;

/** Scalar function-space specialization of the optimization gradient proxy. */
template<typename SpaceType>
class OptimizationGradientForm<SpaceType, false>
{
public:
    using space_type = SpaceType; ///< Function-space type used by the form.
    using space_ptrtype = std::shared_ptr<space_type>; ///< Shared function-space ownership type.
    using value_type = typename space_type::value_type; ///< Scalar value type.
    using vector_type = Vector<value_type>; ///< Algebraic gradient vector type.
    using vector_ptrtype = std::shared_ptr<vector_type>; ///< Borrowed vector view type.

    /**
     * Construct a linear-form proxy around a borrowed gradient vector.
     *
     * @param space test function space
     * @param vector vector filled by subsequent variational assignments
     */
    OptimizationGradientForm( space_ptrtype space, vector_type& vector )
        :
        M_space( std::move( space ) ),
        M_vector( std::addressof( vector ), []( vector_type* ) noexcept {} )
    {}

    /**
     * Construct a scalar proxy with the uniform adapter constructor signature.
     *
     * @tparam BackendPtrType shared algebra-backend pointer type
     * @param space test function space
     * @param vector vector filled by subsequent variational assignments
     * @param backend unused scalar-space backend
     */
    template<typename BackendPtrType>
    OptimizationGradientForm( space_ptrtype space, vector_type& vector,
                              BackendPtrType const& backend )
        : OptimizationGradientForm( std::move( space ), vector )
    {
        static_cast<void>( backend );
    }

    /**
     * Assemble an expression into the gradient vector.
     *
     * @tparam Expression Feel++ variational expression type
     * @param expression expression assigned to a `form1`
     * @return this proxy, enabling ordinary assignment semantics
     */
    template<typename Expression>
    OptimizationGradientForm& operator=( Expression&& expression )
    {
        form1( _test = M_space, _vector = M_vector, _init = true ) =
            std::forward<Expression>( expression );
        return *this;
    }

    /** @return borrowed algebraic vector for advanced assembly operations */
    [[nodiscard]] vector_type& vector() const noexcept { return *M_vector; }

    /** Close the assembled gradient vector. */
    void close() { M_vector->close(); }

private:
    space_ptrtype M_space; ///< Test space retained while the proxy is alive.
    vector_ptrtype M_vector; ///< Non-owning shared view of the callback gradient.
};

/**
 * Writable block linear form backed by a product-space optimization gradient.
 *
 * Component access forwards to `blockform1`, so callbacks use expressions such
 * as `gradient( 0_c ) = integrate(...)` while assembly writes directly into the
 * monolithic vector supplied by the algebraic optimization callback.
 *
 * @tparam SpaceType Feel++ static or dynamic product-space type
 */
template<typename SpaceType>
class OptimizationGradientForm<SpaceType, true>
{
public:
    using space_type = SpaceType; ///< Product-space type used by the form.
    using space_ptrtype = std::shared_ptr<space_type>; ///< Shared product-space ownership type.
    using value_type = typename space_type::value_type; ///< Scalar value type.
    using vector_type = Vector<value_type>; ///< Algebraic gradient vector type.
    using vector_ptrtype = std::shared_ptr<vector_type>; ///< Borrowed vector view type.
    using condensed_vector_type = VectorCondensed<value_type>; ///< Monolithic block-vector wrapper type.
    using condensed_vector_ptrtype = std::shared_ptr<condensed_vector_type>; ///< Shared block-vector wrapper type.
    using block_form_type = BlockLinearForm<space_type const&>; ///< Product-space linear form type.

    /**
     * Construct a product-space linear-form proxy.
     *
     * @tparam BackendPtrType shared algebra-backend pointer type
     * @param space product test space
     * @param vector monolithic gradient vector borrowed from the callback
     * @param backend backend used to build the product-space block layout
     */
    template<typename BackendPtrType>
    OptimizationGradientForm( space_ptrtype space, vector_type& vector,
                              BackendPtrType const& backend )
        :
        M_space( std::move( space ) ),
        M_vector( std::addressof( vector ), []( vector_type* ) noexcept {} ),
        M_condensed( std::make_shared<condensed_vector_type>(
            solve::strategy::monolithic, blockVector( *M_space, backend ),
            backend, false ) ),
        M_form( std::as_const( *M_space ), M_condensed )
    {
        M_condensed->vectorPtr() = M_vector;
        M_form.zero();
    }

    /**
     * Access one product-space gradient block.
     *
     * @tparam Args compile-time or runtime block-index argument types
     * @param args block index and optional dynamic subspace index
     * @return writable scalar `form1` view for the requested block
     */
    template<typename... Args>
    decltype(auto) operator()( Args&&... args )
    {
        return M_form( std::forward<Args>( args )... );
    }

    /** Set all gradient blocks and the monolithic callback vector to zero. */
    void zero() { M_form.zero(); }

    /** Close the block form and its monolithic callback vector. */
    void close() { M_form.close(); }

    /** @return borrowed monolithic algebraic gradient vector */
    [[nodiscard]] vector_type& vector() const noexcept { return *M_vector; }

    /** @return writable underlying product-space block form */
    [[nodiscard]] block_form_type& blockForm() noexcept { return M_form; }

private:
    space_ptrtype M_space; ///< Product test space retained by the proxy.
    vector_ptrtype M_vector; ///< Non-owning shared view of the callback gradient.
    condensed_vector_ptrtype M_condensed; ///< Block layout bound to the callback vector.
    block_form_type M_form; ///< Writable block linear form.
};

/**
 * Writable variational bilinear form backed by an optimization Hessian matrix.
 *
 * This non-owning proxy lets a field callback assemble an exact Hessian with
 * `hessian = integrate(...)`. The borrowed matrix remains valid only during
 * the callback that receives this proxy.
 *
 * @tparam SpaceType Feel++ function-space type
 */
template<typename SpaceType,
         bool IsProductSpace = Feel::detail::isOptimizationProductSpace<SpaceType>>
class OptimizationHessianForm
;

/** Scalar function-space specialization of the optimization Hessian proxy. */
template<typename SpaceType>
class OptimizationHessianForm<SpaceType, false>
{
public:
    using space_type = SpaceType; ///< Function-space type used by the form.
    using space_ptrtype = std::shared_ptr<space_type>; ///< Shared function-space ownership type.
    using value_type = typename space_type::value_type; ///< Scalar value type.
    using matrix_type = MatrixSparse<value_type>; ///< Algebraic Hessian matrix type.
    using matrix_ptrtype = std::shared_ptr<matrix_type>; ///< Borrowed matrix view type.

    /**
     * Construct a bilinear-form proxy around a borrowed Hessian matrix.
     *
     * @param space trial and test function space
     * @param matrix matrix filled by subsequent variational assignments
     */
    OptimizationHessianForm( space_ptrtype space, matrix_type& matrix )
        :
        M_space( std::move( space ) ),
        M_matrix( std::addressof( matrix ), []( matrix_type* ) noexcept {} )
    {}

    /**
     * Construct a scalar proxy with the uniform adapter constructor signature.
     *
     * @tparam BackendPtrType shared algebra-backend pointer type
     * @param space trial and test function space
     * @param matrix matrix filled by subsequent variational assignments
     * @param backend unused scalar-space backend
     */
    template<typename BackendPtrType>
    OptimizationHessianForm( space_ptrtype space, matrix_type& matrix,
                             BackendPtrType const& backend )
        : OptimizationHessianForm( std::move( space ), matrix )
    {
        static_cast<void>( backend );
    }

    /**
     * Assemble an expression into the Hessian matrix.
     *
     * @tparam Expression Feel++ variational expression type
     * @param expression expression assigned to a `form2`
     * @return this proxy, enabling ordinary assignment semantics
     */
    template<typename Expression>
    OptimizationHessianForm& operator=( Expression&& expression )
    {
        form2( _trial = M_space, _test = M_space, _matrix = M_matrix,
               _init = true ) = std::forward<Expression>( expression );
        return *this;
    }

    /** @return borrowed algebraic matrix for advanced assembly operations */
    [[nodiscard]] matrix_type& matrix() const noexcept { return *M_matrix; }

    /** Close the assembled scalar Hessian. */
    void close() { M_matrix->close(); }

private:
    space_ptrtype M_space; ///< Trial and test space retained by the proxy.
    matrix_ptrtype M_matrix; ///< Non-owning shared view of the callback Hessian.
};

/**
 * Writable block bilinear form backed by a product-space optimization Hessian.
 *
 * Component access forwards to `blockform2`. After user assembly, `close()`
 * materializes the block form into the matrix supplied by the common algebraic
 * Hessian callback.
 *
 * @tparam SpaceType Feel++ static or dynamic product-space type
 */
template<typename SpaceType>
class OptimizationHessianForm<SpaceType, true>
{
public:
    using space_type = SpaceType; ///< Product-space type used by the form.
    using space_ptrtype = std::shared_ptr<space_type>; ///< Shared product-space ownership type.
    using value_type = typename space_type::value_type; ///< Scalar value type.
    using matrix_type = MatrixSparse<value_type>; ///< Algebraic Hessian matrix type.
    using matrix_ptrtype = std::shared_ptr<matrix_type>; ///< Borrowed matrix view type.
    using condensed_matrix_type = MatrixCondensed<value_type>; ///< Monolithic block-matrix wrapper type.
    using condensed_matrix_ptrtype = std::shared_ptr<condensed_matrix_type>; ///< Shared block-matrix wrapper type.
    using block_form_type = BlockBilinearForm<space_type const&>; ///< Product-space bilinear form type.

    /**
     * Construct a product-space bilinear-form proxy.
     *
     * @tparam BackendPtrType shared algebra-backend pointer type
     * @param space product trial and test space
     * @param matrix monolithic Hessian matrix borrowed from the callback
     * @param backend backend used to allocate the block form
     */
    template<typename BackendPtrType>
    OptimizationHessianForm( space_ptrtype space, matrix_type& matrix,
                             BackendPtrType const& backend )
        :
        M_space( std::move( space ) ),
        M_matrix( std::addressof( matrix ), []( matrix_type* ) noexcept {} ),
        M_condensed( std::make_shared<condensed_matrix_type>(
            solve::strategy::monolithic,
            csrGraphBlocks( std::as_const( *M_space ),
                            std::as_const( *M_space ), Pattern::COUPLED ),
            backend, true ) ),
        M_form( std::as_const( *M_space ), M_condensed )
    {
        M_form.matrix().zero();
    }

    /**
     * Access one product-space Hessian block.
     *
     * @tparam Args compile-time or runtime row/column block-index types
     * @param args row and column block indices
     * @return writable scalar `form2` view for the requested block
     */
    template<typename... Args>
    decltype(auto) operator()( Args&&... args )
    {
        return M_form( std::forward<Args>( args )... );
    }

    /** Set every Hessian block to zero. */
    void zero() { M_form.matrix().zero(); }

    /**
     * Close the block form and copy it into the algebraic callback matrix.
     */
    void close()
    {
        M_form.close();
        M_matrix->zero();
        M_matrix->addMatrix( value_type( 1 ), *M_form.baseMatrixPtr(),
                             SAME_NONZERO_PATTERN );
        M_matrix->close();
    }

    /** @return borrowed monolithic algebraic Hessian matrix */
    [[nodiscard]] matrix_type& matrix() const noexcept { return *M_matrix; }

    /** @return writable underlying product-space block form */
    [[nodiscard]] block_form_type& blockForm() noexcept { return M_form; }

private:
    space_ptrtype M_space; ///< Product trial and test space retained by the proxy.
    matrix_ptrtype M_matrix; ///< Non-owning shared view of the callback Hessian.
    condensed_matrix_ptrtype M_condensed; ///< Block layout used for Hessian assembly.
    block_form_type M_form; ///< Writable block bilinear form.
};

/**
 * Field-native adapter for the backend-independent optimization solver.
 *
 * The adapter translates algebraic TAO callbacks into callbacks receiving
 * Feel++ finite-element fields and writable variational-form proxies. It also
 * owns the algebraic work vectors needed to synchronize ghosts and solves
 * directly from a function-space element.
 *
 * The adapter is a convenience layer around `SolverOptimization`; it is not a
 * PDE problem or sensitivity abstraction. Use algebraicSolver() whenever
 * backend-independent low-level controls are required.
 *
 * @tparam SpaceType Feel++ function-space type
 */
template<typename SpaceType>
class SolverOptimizationField
{
public:
    using self_type = SolverOptimizationField<SpaceType>; ///< Current adapter type.
    using space_type = SpaceType; ///< Optimized function-space type.
    using space_ptrtype = std::shared_ptr<space_type>; ///< Shared function-space ownership type.
    using value_type = typename space_type::value_type; ///< Optimization scalar type.
    using real_type = typename type_traits<value_type>::real_type; ///< Real scalar type.
    using element_type = typename space_type::element_type; ///< Optimized finite-element field type.
    using backend_type = Backend<value_type>; ///< Feel++ algebra backend type.
    using backend_ptrtype = std::shared_ptr<backend_type>; ///< Shared algebra backend type.
    using solver_type = SolverOptimization<value_type>; ///< Backend-independent solver type.
    using solver_ptrtype = typename solver_type::ptrtype; ///< Shared optimization solver type.
    using vector_type = typename solver_type::vector_type; ///< Algebraic optimization vector type.
    using vector_ptrtype = typename solver_type::vector_ptrtype; ///< Shared optimization vector type.
    using matrix_type = typename solver_type::sparse_matrix_type; ///< Exact Hessian matrix type.
    using matrix_ptrtype = typename solver_type::sparse_matrix_ptrtype; ///< Shared Hessian matrix type.
    using result_type = typename solver_type::result_type; ///< Structured solve result type.
    using monitor_record_type = typename solver_type::monitor_record_type; ///< Structured iteration record type.
    using monitor_id_type = typename solver_type::monitor_id_type; ///< Stable monitor registration identifier.
    using objective_gradient_data_type =
        typename solver_type::objective_gradient_data_type; ///< Algebraic objective-gradient data.
    using hessian_data_type = typename solver_type::hessian_data_type; ///< Algebraic Hessian data.
    using gradient_form_type = OptimizationGradientForm<space_type>; ///< Writable gradient form proxy.
    using hessian_form_type = OptimizationHessianForm<space_type>; ///< Writable Hessian form proxy.

    /**
     * Construct a field-native optimization adapter.
     *
     * @param space function space containing the optimization variable
     * @param backendName algebra and optimization backend identifier
     * @param prefix independent runtime-options prefix
     * @param worldComm communicator used collectively by the solver
     * @param vm Feel++ program-options values
     */
    explicit SolverOptimizationField(
        space_ptrtype space, std::string backendName = "petsc",
        std::string prefix = {},
        worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
        po::variables_map const& vm = Environment::vm() )
        :
        M_space( requireSpace( std::move( space ) ) ),
        M_backend( Feel::backend( _kind = backendName,
                                  _worldcomm = worldComm, _rebuild = false ) ),
        M_solver( solver_type::build( backendName, prefix, worldComm, vm ) ),
        M_workspace( std::make_shared<FieldWorkspace>(
            M_space, M_backend, newFieldVector( M_space, M_backend ) ) ),
        M_variable( newFieldVector( M_space, M_backend ) )
    {}

    /** @return optimized function space */
    [[nodiscard]] space_ptrtype const& space() const noexcept { return M_space; }

    /** @return algebra backend used for vectors and matrices */
    [[nodiscard]] backend_ptrtype const& algebraBackend() const noexcept
    {
        return M_backend;
    }

    /** @return wrapped backend-independent optimization solver */
    [[nodiscard]] solver_ptrtype const& algebraicSolver() const noexcept
    {
        return M_solver;
    }

    /**
     * Select the optimization algorithm.
     *
     * @param type backend solver type such as `lmvm`, `nls`, or `blmvm`
     * @return this adapter for fluent configuration
     */
    self_type& algorithm( std::string type )
    {
        M_solver->setType( std::move( type ) );
        return *this;
    }

    /**
     * Set the optimization iteration limit.
     *
     * @param value maximum iteration count
     * @return this adapter for fluent configuration
     */
    self_type& maxIterations( typename solver_type::size_type value )
    {
        M_solver->setMaximumIterations( value );
        return *this;
    }

    /**
     * Set one value for all three gradient convergence tolerances.
     *
     * @param value absolute, relative, and initial-gradient tolerance
     * @return this adapter for fluent configuration
     */
    self_type& gradientTolerance( real_type value )
    {
        M_solver->setGradientTolerances( value, value, value );
        return *this;
    }

    /**
     * Set distinct gradient convergence tolerances.
     *
     * @param absolute absolute gradient tolerance
     * @param relative objective-relative gradient tolerance
     * @param initial initial-gradient-relative tolerance
     * @return this adapter for fluent configuration
     */
    self_type& gradientTolerances( real_type absolute, real_type relative,
                                   real_type initial )
    {
        M_solver->setGradientTolerances( absolute, relative, initial );
        return *this;
    }

    /**
     * Set constant lower and upper bounds for every field degree of freedom.
     *
     * @param lower lower bound
     * @param upper upper bound
     * @return this adapter for fluent configuration
     */
    self_type& bounds( real_type lower, real_type upper )
    {
        M_solver->setConstantLowerBound( lower );
        M_solver->setConstantUpperBound( upper );
        return *this;
    }

    /**
     * Enable or disable the standard backend iteration monitor.
     *
     * @param enabled true to print iteration records
     * @return this adapter for fluent configuration
     */
    self_type& monitor( bool enabled = true )
    {
        M_solver->setMonitorEnabled( enabled );
        return *this;
    }

    /**
     * Register a structured C++ convergence monitor.
     *
     * PETSc invokes the monitor collectively during the algebraic TAO solve.
     * The callback receives objective, gradient residual, step norm, and status
     * information without exposing a native `Tao` handle.
     *
     * @tparam Callable copyable callable accepting `monitor_record_type const&`
     * @param callback callback invoked for every reported TAO iteration
     * @return stable identifier accepted by removeMonitor()
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&,
                                monitor_record_type const&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&,
                                          monitor_record_type const&>,
                     void>
    [[nodiscard]] monitor_id_type addMonitor( Callable&& callback )
    {
        return M_solver->addMonitor( std::forward<Callable>( callback ) );
    }

    /**
     * Remove one structured C++ convergence monitor.
     *
     * @param identifier value returned by addMonitor()
     * @return true when a registered monitor was removed
     */
    bool removeMonitor( monitor_id_type identifier )
    {
        return M_solver->removeMonitor( identifier );
    }

    /** Remove every structured C++ convergence monitor. */
    void resetMonitors() noexcept { M_solver->resetMonitors(); }

    /** @return true when at least one structured C++ monitor is registered */
    [[nodiscard]] bool hasMonitors() const noexcept
    {
        return M_solver->hasMonitors();
    }

    /**
     * Register a field-native objective callback.
     *
     * @tparam Callable copyable callable receiving `element_type const&`
     * @param callback callback returning the objective value
     * @return this adapter for fluent configuration
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, element_type const&> &&
                 std::convertible_to<
                     std::invoke_result_t<std::decay_t<Callable>&,
                                          element_type const&>, real_type>
    self_type& objective( Callable&& callback )
    {
        M_solver->setObjective(
            [workspace = M_workspace,
             fn = std::forward<Callable>( callback )]( vector_type const& state ) mutable
            {
                auto field = workspace->field( state );
                return static_cast<real_type>( std::invoke( fn, field ) );
            } );
        return *this;
    }

    /**
     * Register a field-native gradient callback.
     *
     * The callback receives the current field and a writable proxy supporting
     * `gradient = integrate(...)`.
     *
     * @tparam Callable copyable field-gradient callable
     * @param callback callback assembling the gradient form
     * @return this adapter for fluent configuration
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, element_type const&,
                                gradient_form_type&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&,
                                          element_type const&,
                                          gradient_form_type&>, void>
    self_type& gradient( Callable&& callback )
    {
        M_solver->setGradient(
            [workspace = M_workspace, space = M_space,
             fn = std::forward<Callable>( callback )](
                vector_type const& state, vector_type& gradientVector ) mutable
            {
                auto field = workspace->field( state );
                gradient_form_type gradientForm( space, gradientVector,
                                                 workspace->backend );
                std::invoke( fn, field, gradientForm );
                gradientForm.close();
            } );
        return *this;
    }

    /**
     * Register a combined field-native objective-and-gradient callback.
     *
     * The callback receives the current field and a writable proxy supporting
     * `gradient = integrate(...)`, then returns the scalar objective.
     *
     * @tparam Callable copyable combined callback
     * @param callback callback evaluating the objective and assembling its gradient
     * @return this adapter for fluent configuration
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, element_type const&,
                                gradient_form_type&> &&
                 std::convertible_to<
                     std::invoke_result_t<std::decay_t<Callable>&,
                                          element_type const&,
                                          gradient_form_type&>, real_type>
    self_type& objectiveGradient( Callable&& callback )
    {
        M_solver->setObjectiveGradient(
            [workspace = M_workspace, space = M_space,
             fn = std::forward<Callable>( callback )](
                objective_gradient_data_type& data ) mutable
            {
                auto field = workspace->field( data.state );
                gradient_form_type gradientForm( space, data.gradient,
                                                 workspace->backend );
                data.objective = static_cast<real_type>(
                    std::invoke( fn, field, gradientForm ) );
                gradientForm.close();
            } );
        return *this;
    }

    /**
     * Register an exact field-native Hessian callback with an explicit matrix.
     *
     * The callback receives the current field and a writable proxy supporting
     * `hessian = integrate(...)`.
     *
     * @tparam Callable copyable field-Hessian callable
     * @param matrix matrix retained and reused by the optimization solver
     * @param callback callback assembling the exact Hessian form
     * @return this adapter for fluent configuration
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, element_type const&,
                                hessian_form_type&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&,
                                          element_type const&,
                                          hessian_form_type&>, void>
    self_type& hessian( matrix_ptrtype matrix, Callable&& callback )
    {
        if constexpr ( Feel::detail::isOptimizationProductSpace<space_type> )
        {
            auto hessianForm = std::make_shared<hessian_form_type>(
                M_space, *matrix, M_backend );
            M_solver->setHessian(
                std::move( matrix ),
                [workspace = M_workspace, hessianForm,
                 fn = std::forward<Callable>( callback )](
                    hessian_data_type& data ) mutable
                {
                    auto field = workspace->field( data.state );
                    hessianForm->zero();
                    std::invoke( fn, field, *hessianForm );
                    hessianForm->close();
                } );
        }
        else
        {
            M_solver->setHessian(
                std::move( matrix ),
                [workspace = M_workspace, space = M_space,
                 fn = std::forward<Callable>( callback )](
                    hessian_data_type& data ) mutable
                {
                    auto field = workspace->field( data.state );
                    hessian_form_type hessianForm( space, data.hessian,
                                                   workspace->backend );
                    std::invoke( fn, field, hessianForm );
                    hessianForm.close();
                } );
        }
        return *this;
    }

    /**
     * Register an exact Hessian and allocate its matrix from this field space.
     *
     * This documented convenience overload allocates one matrix and reuses it
     * as both the exact Hessian and its preconditioner.
     *
     * @tparam Callable copyable field-Hessian callable
     * @param callback callback assembling the exact Hessian form
     * @return this adapter for fluent configuration
     */
    template<typename Callable>
        requires std::copy_constructible<std::decay_t<Callable>> &&
                 std::invocable<std::decay_t<Callable>&, element_type const&,
                                hessian_form_type&> &&
                 std::same_as<
                     std::invoke_result_t<std::decay_t<Callable>&,
                                          element_type const&,
                                          hessian_form_type&>, void>
    self_type& hessian( Callable&& callback )
    {
        auto matrix = newFieldMatrix( M_space, M_backend );
        return this->hessian( std::move( matrix ),
                              std::forward<Callable>( callback ) );
    }

    /**
     * Solve in place from an initial finite-element field.
     *
     * The initial field is copied into the backend vector before TAO starts;
     * the final vector is synchronized and copied back into the same field.
     *
     * @param variable initial field on input and optimized field on output
     * @return structured optimization result
     */
    [[nodiscard]] result_type solve( element_type& variable )
    {
        copyElementToVector( variable, *M_variable );
        auto result = M_solver->solve( M_variable );
        M_variable->close();
        copyVectorToElement( M_space, M_variable, variable );
        return result;
    }

private:
    /**
     * Validate the function space before dependent work vectors are created.
     *
     * @param space candidate optimization space
     * @return the validated space
     * @throws std::invalid_argument when @p space is null
     */
    [[nodiscard]] static space_ptrtype requireSpace( space_ptrtype space )
    {
        if ( !space )
            throw std::invalid_argument(
                "SolverOptimizationField requires a non-null function space" );
        return space;
    }

    /**
     * Allocate monolithic algebraic storage for a scalar or product space.
     *
     * @param space optimized scalar or product function space
     * @param backend algebra backend used for allocation
     * @return vector with the exact distributed layout of the field
     */
    [[nodiscard]] static vector_ptrtype
    newFieldVector( space_ptrtype const& space, backend_ptrtype const& backend )
    {
        if constexpr ( Feel::detail::isOptimizationProductSpace<space_type> )
        {
            auto vector = std::make_shared<VectorCondensed<value_type>>(
                solve::strategy::monolithic, blockVector( *space, backend ),
                backend, false );
            return vector->vectorPtr();
        }
        else
            return backend->newVector( space );
    }

    /**
     * Allocate an exact Hessian matrix for a scalar or product space.
     *
     * @param space optimized scalar or product function space
     * @param backend algebra backend used for allocation
     * @return square matrix matching the monolithic optimization vector
     */
    [[nodiscard]] static matrix_ptrtype
    newFieldMatrix( space_ptrtype const& space, backend_ptrtype const& backend )
    {
        if constexpr ( Feel::detail::isOptimizationProductSpace<space_type> )
        {
            auto matrix = std::make_shared<MatrixCondensed<value_type>>(
                solve::strategy::monolithic,
                csrGraphBlocks( std::as_const( *space ),
                                std::as_const( *space ), Pattern::COUPLED ),
                backend, true );
            return matrix->getSparseMatrix();
        }
        else
            return backend->newMatrix( _trial = space, _test = space );
    }

    /**
     * Pack a scalar or block field into its monolithic algebraic vector.
     *
     * @param field finite-element optimization variable
     * @param vector destination algebraic vector
     */
    static void copyElementToVector( element_type const& field,
                                     vector_type& vector )
    {
        if constexpr ( Feel::detail::isOptimizationProductSpace<space_type> )
            field.updateVectorFromSubVectors( vector );
        else
        {
            vector = field;
            vector.close();
        }
    }

    /**
     * Unpack a monolithic vector into a scalar or block field.
     *
     * @param space function space used to construct scalar elements
     * @param vector source monolithic algebraic vector
     * @param field destination finite-element field
     */
    static void copyVectorToElement( space_ptrtype const& space,
                                     vector_ptrtype const& vector,
                                     element_type& field )
    {
        if constexpr ( Feel::detail::isOptimizationProductSpace<space_type> )
            field.localize( vector );
        else
            field = space->element( static_cast<vector_type const&>( *vector ) );
    }

    /** Shared state that safely survives moves of the fluent adapter. */
    struct FieldWorkspace
    {
        /**
         * Construct callback synchronization storage.
         *
         * @param fieldSpace optimized function space
         * @param algebraBackend backend used by callback form proxies
         * @param stateVector ghosted backend work vector
         */
        FieldWorkspace( space_ptrtype fieldSpace, backend_ptrtype algebraBackend,
                        vector_ptrtype stateVector )
            :
            space( std::move( fieldSpace ) ),
            backend( std::move( algebraBackend ) ),
            state( std::move( stateVector ) )
        {}

        /**
         * Synchronize an algebraic state and expose it as a field snapshot.
         *
         * @param source algebraic optimization state
         * @return finite-element field containing current local and ghost values
         */
        [[nodiscard]] element_type field( vector_type const& source )
        {
            *state = source;
            state->close();
            if constexpr ( Feel::detail::isOptimizationProductSpace<space_type> )
            {
                auto result = space->element();
                result.localize( state );
                return result;
            }
            else
                return space->element( static_cast<vector_type const&>( *state ) );
        }

        space_ptrtype space; ///< Function space used to create field snapshots.
        backend_ptrtype backend; ///< Algebra backend used by variational proxies.
        vector_ptrtype state; ///< Ghost-synchronized callback work vector.
    };

    space_ptrtype M_space; ///< Function space containing the optimized field.
    backend_ptrtype M_backend; ///< Algebra backend used for field storage.
    solver_ptrtype M_solver; ///< Wrapped backend-independent optimization solver.
    std::shared_ptr<FieldWorkspace> M_workspace; ///< Shared callback synchronization state.
    vector_ptrtype M_variable; ///< Persistent algebraic vector used by solve().
};

/**
 * Build a field-native optimization solver with Feel++ named arguments.
 *
 * Supported named arguments are `_space` (required), `_backend` (default
 * `petsc`), `_name`, `_worldcomm` (defaulting to the space communicator), and
 * `_vm`.
 *
 * @tparam Ts named-argument types
 * @param arguments Feel++ named arguments
 * @return fluent field-native optimization adapter
 */
template<typename... Ts>
[[nodiscard]] auto optimizationField( Ts&&... arguments )
{
    auto args = NA::make_arguments( std::forward<Ts>( arguments )... );
    auto space = args.get( _space );
    using space_ptrtype = std::remove_cvref_t<decltype( space )>;
    using space_type = typename space_ptrtype::element_type;

    std::string const backendName =
        args.get_else( _backend, std::string( "petsc" ) );
    std::string const name = args.get_else( _name, std::string{} );
    worldcomm_ptr_t const worldComm = args.get_else_invocable(
        _worldcomm,
        [&space]() { return Feel::detail::optimizationSpaceWorldComm( space ); } );
    po::variables_map const& vm = args.get_else( _vm, Environment::vm() );

    return SolverOptimizationField<space_type>(
        std::move( space ), backendName, name, worldComm, vm );
}

} // namespace Feel

#endif // FEELPP_FEELOPT_SOLVEROPTIMIZATIONFIELD_HPP
