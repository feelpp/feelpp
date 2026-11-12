/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_NONLINEAR_HPP
#define FEELPP_FEELSYSTEM_NONLINEAR_HPP 1

/**
 * @file nonlinear.hpp
 * @brief Residual/Jacobian nonlinear systems backed by Feel++ nonlinear solvers.
 */

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

#include <feel/feelalg/backend.hpp>
#include <feel/feelsystem/backend.hpp>
#include <feel/feelsystem/concepts.hpp>
#include <feel/feelsystem/namedarguments.hpp>
#include <feel/feelsystem/result.hpp>

namespace Feel::systems
{

namespace detail
{
/**
 * @brief RAII guard restoring nonlinear solver callbacks on scope exit.
 */
template <typename SolverType>
class NonlinearCallbackGuard
{
  public:
    /** @brief Concrete nonlinear solver type whose callbacks are saved. */
    using solver_type = std::remove_reference_t<SolverType>;
    /** @brief Stored residual callback type. */
    using residual_function_type = std::remove_cvref_t<decltype( std::declval<solver_type&>().residual )>;
    /** @brief Stored Jacobian callback type. */
    using jacobian_function_type = std::remove_cvref_t<decltype( std::declval<solver_type&>().jacobian )>;

    /**
     * @brief Save the current callbacks of `solver`.
     */
    explicit NonlinearCallbackGuard( solver_type& solver )
        : M_solver( &solver ),
          M_residual( solver.residual ),
          M_jacobian( solver.jacobian )
    {
    }

    NonlinearCallbackGuard( NonlinearCallbackGuard const& ) = delete;
    NonlinearCallbackGuard& operator=( NonlinearCallbackGuard const& ) = delete;

    NonlinearCallbackGuard( NonlinearCallbackGuard&& other ) noexcept
        : M_solver( std::exchange( other.M_solver, nullptr ) ),
          M_residual( std::move( other.M_residual ) ),
          M_jacobian( std::move( other.M_jacobian ) )
    {
    }

    NonlinearCallbackGuard& operator=( NonlinearCallbackGuard&& ) = delete;

    ~NonlinearCallbackGuard()
    {
        this->restore();
    }

    /**
     * @brief Restore saved callbacks early.
     */
    void restore()
    {
        if ( !M_solver )
            return;
        M_solver->residual = std::move( M_residual );
        M_solver->jacobian = std::move( M_jacobian );
        M_solver = nullptr;
    }

  private:
    /** @brief Solver whose callbacks will be restored, or null after restore. */
    solver_type* M_solver = nullptr;
    /** @brief Saved residual callback. */
    residual_function_type M_residual;
    /** @brief Saved Jacobian callback. */
    jacobian_function_type M_jacobian;
};

template <typename SolverType>
NonlinearCallbackGuard( SolverType& ) -> NonlinearCallbackGuard<SolverType>;

/**
 * @brief Minimal scope-exit helper for callback lifetime cleanup.
 */
template <typename Function>
class ScopeExit
{
  public:
    /** @brief Store a cleanup function to run at scope exit. */
    explicit ScopeExit( Function function )
        : M_function( std::move( function ) )
    {
    }

    ScopeExit( ScopeExit const& ) = delete;
    ScopeExit& operator=( ScopeExit const& ) = delete;

    ScopeExit( ScopeExit&& other ) noexcept
        : M_function( std::move( other.M_function ) ),
          M_active( std::exchange( other.M_active, false ) )
    {
    }

    ScopeExit& operator=( ScopeExit&& ) = delete;

    ~ScopeExit()
    {
        if ( M_active )
            M_function();
    }

    /** @brief Disable the cleanup function. */
    void release() noexcept { M_active = false; }

  private:
    /** @brief Cleanup function. */
    Function M_function;
    /** @brief True while cleanup should still run. */
    bool M_active = true;
};

/**
 * @brief Build a `ScopeExit` object with decayed function storage.
 */
template <typename Function>
ScopeExit<std::remove_cvref_t<Function>>
scopeExit( Function&& function )
{
    return ScopeExit<std::remove_cvref_t<Function>>{ std::forward<Function>( function ) };
}
} // namespace detail

/**
 * @brief Stateful assembled nonlinear system using residual and Jacobian callbacks.
 *
 * `NonlinearSystem` combines affine matrix/rhs contributions with explicit
 * nonlinear residual and Jacobian terms. During `solve()`, it temporarily
 * installs callbacks into the backend nonlinear solver and restores any
 * previous callbacks before returning.
 *
 * @tparam ProblemType Problem descriptor satisfying `NonlinearProblemLike`.
 */
template <NonlinearProblemLike ProblemType>
class NonlinearSystem
{
  public:
    /** @brief Decayed problem descriptor type. */
    using problem_type = std::remove_cvref_t<ProblemType>;
    /** @brief Function-space type used by the problem. */
    using space_type = typename problem_type::space_type;
    /** @brief Scalar value type used by algebraic objects. */
    using value_type = typename problem_type::value_type;
    /** @brief Finite-element state type. */
    using element_type = typename problem_type::element_type;
    /** @brief Feel++ backend type used for allocation and nonlinear solves. */
    using backend_type = Backend<value_type>;
    /** @brief Shared backend pointer type. */
    using backend_ptrtype = std::shared_ptr<backend_type>;
    /** @brief Sparse matrix pointer type allocated by the backend. */
    using matrix_ptrtype = typename backend_type::sparse_matrix_ptrtype;
    /** @brief Vector pointer type allocated by the backend. */
    using vector_ptrtype = typename backend_type::vector_ptrtype;
    /** @brief Concrete sparse matrix type. */
    using matrix_type = typename matrix_ptrtype::element_type;
    /** @brief Concrete vector type. */
    using vector_type = typename vector_ptrtype::element_type;

    /**
     * @brief Build a nonlinear system from a reusable problem descriptor.
     *
     * If `backend` is empty, `assemble()` creates or reuses a backend whose
     * prefix is derived from `backendName` or from the problem name.
     */
    NonlinearSystem( problem_type const& problem, backend_ptrtype backend = {}, std::string backendName = {} )
        : M_problem( problem ),
          M_backend( std::move( backend ) ),
          M_backendName( detail::defaultBackendName( std::move( backendName ), M_problem.name(), "nonlinear_system" ) )
    {
        if ( M_backend )
            M_backendName = M_backend->prefix();
    }

    /**
     * @brief Allocate algebraic storage and assemble affine contributions.
     */
    NonlinearSystem& assemble()
    {
        if ( !M_backend )
            M_backend = detail::makeBackend<value_type>( M_backendName, M_problem.space()->worldCommPtr() );
        if ( !M_jacobian )
            M_jacobian = M_backend->newMatrix( _test = M_problem.space(), _trial = M_problem.space() );
        if ( !M_residual )
            M_residual = M_backend->newVector( _test = M_problem.space() );
        if ( !M_affineMatrix )
            M_affineMatrix = M_backend->newMatrix( _test = M_problem.space(), _trial = M_problem.space() );
        if ( !M_affineRhs )
            M_affineRhs = M_backend->newVector( _test = M_problem.space() );

        M_affineMatrix->zero();
        M_affineRhs->zero();
        M_problem.assembleAffine( M_affineMatrix, M_affineRhs );
        M_affineMatrix->close();
        M_affineRhs->close();
        M_assembled = true;
        return *this;
    }

    /**
     * @brief Assemble `F(x)` in the owned residual vector.
     */
    NonlinearSystem& assembleResidual( vector_ptrtype const& x )
    {
        if ( !M_assembled )
            this->assemble();
        this->updateResidual( x, M_residual );
        return *this;
    }

    /**
     * @brief Assemble `J(x)` in the owned Jacobian matrix.
     */
    NonlinearSystem& assembleJacobian( vector_ptrtype const& x )
    {
        if ( !M_assembled )
            this->assemble();
        this->updateJacobian( x, M_jacobian );
        return *this;
    }

    /** @brief Return true once affine storage has been allocated and assembled. */
    bool assembled() const noexcept { return M_assembled; }
    /** @brief Problem descriptor copied into this system. */
    problem_type const& problem() const noexcept { return M_problem; }
    /** @brief Backend used for allocation and solving. */
    backend_ptrtype const& backend() const noexcept { return M_backend; }
    /** @brief Backend prefix used for options and default backend lookup. */
    std::string const& backendName() const noexcept { return M_backendName; }

    /** @brief Owned Jacobian matrix storage. */
    matrix_ptrtype const& jacobianPtr() const noexcept { return M_jacobian; }
    /** @brief Owned residual vector storage. */
    vector_ptrtype const& residualPtr() const noexcept { return M_residual; }
    /** @brief Dereferenced owned Jacobian matrix. */
    matrix_type const& jacobian() const noexcept { return *M_jacobian; }
    /** @brief Dereferenced owned residual vector. */
    vector_type const& residual() const noexcept { return *M_residual; }
    /** @brief Affine matrix contribution assembled before nonlinear callbacks. */
    matrix_ptrtype const& affineMatrixPtr() const noexcept { return M_affineMatrix; }
    /** @brief Affine rhs contribution assembled before nonlinear callbacks. */
    vector_ptrtype const& affineRhsPtr() const noexcept { return M_affineRhs; }

    /**
     * @brief Backend nonlinear residual callback.
     */
    void updateResidual( vector_ptrtype const& x, vector_ptrtype& residual )
    {
        this->updateExternalState( x );
        auto state = M_problem.space()->element();
        state = *x;

        residual->zero();
        if ( M_affineMatrix )
            M_affineMatrix->multVector( x, residual );
        if ( M_affineRhs )
            residual->add( -1.0, *M_affineRhs );
        M_problem.assembleResidualTerms( state, residual );
        residual->close();
    }

    /**
     * @brief Backend nonlinear Jacobian callback.
     */
    void updateJacobian( vector_ptrtype const& x, matrix_ptrtype& jacobian )
    {
        this->updateExternalState( x );
        auto state = M_problem.space()->element();
        state = *x;

        jacobian->zero();
        if ( M_affineMatrix )
            jacobian->addMatrix( 1.0, M_affineMatrix );
        M_problem.assembleJacobianTerms( state, jacobian );
        jacobian->close();
    }

    /**
     * @brief Solve the nonlinear system with the backend nonlinear solver.
     *
     * Newton-type methods use explicit residual/Jacobian callbacks. Affine
     * terms contribute automatically as `A u - b` and `A`.
     */
    template <typename... Args>
    solve_result solve( Args&&... args )
    {
        if ( !M_assembled )
            this->assemble();

        static_assert( detail::hasNamedArgument<na::solution, Args...>,
                       "Feel::systems::NonlinearSystem::solve requires _solution" );
        auto&& solution = detail::namedArgumentValue<na::solution>( args... );
        M_updateExternalState = [&solution]( vector_type const& x )
        {
            Feel::detail::ref( solution ) = x;
            Feel::detail::ref( solution ).close();
        };
        auto clearExternalState = detail::scopeExit( [this]()
                                                     { M_updateExternalState = {}; } );

        auto& nonlinearSolver = *M_backend->nlSolver();
        auto callbackGuard = detail::NonlinearCallbackGuard{ nonlinearSolver };
        nonlinearSolver.residual =
            [this]( vector_ptrtype const& x, vector_ptrtype& residual )
        {
            this->updateResidual( x, residual );
        };
        nonlinearSolver.jacobian =
            [this]( vector_ptrtype const& x, matrix_ptrtype& jacobian )
        {
            this->updateJacobian( x, jacobian );
        };

        auto result = M_backend->nlSolve( _jacobian = M_jacobian,
                                          _residual = M_residual,
                                          std::forward<Args>( args )... );
        return solve_result{ result.isConverged(),
                             static_cast<int>( result.nIterations() ),
                             static_cast<double>( result.residual() ),
                             result.isConverged()
                                 ? std::string{ "feelsystem NonlinearSystem solve converged" }
                                 : std::string{ "feelsystem NonlinearSystem solve did not converge" } };
    }

  private:
    /** @brief Synchronize an external finite-element solution with backend vector state. */
    void updateExternalState( vector_ptrtype const& x )
    {
        if ( M_updateExternalState )
            M_updateExternalState( *x );
    }

    /** @brief Problem descriptor copied by value for replayable assembly. */
    problem_type M_problem;
    /** @brief Backend used for allocation and solving. */
    backend_ptrtype M_backend;
    /** @brief Backend prefix selected from user input, backend, or problem name. */
    std::string M_backendName;
    /** @brief Owned Jacobian matrix storage. */
    matrix_ptrtype M_jacobian;
    /** @brief Owned residual vector storage. */
    vector_ptrtype M_residual;
    /** @brief Affine matrix assembled once and reused in residual/Jacobian callbacks. */
    matrix_ptrtype M_affineMatrix;
    /** @brief Affine rhs assembled once and reused in residual callbacks. */
    vector_ptrtype M_affineRhs;
    /** @brief Optional callback updating the user solution object during nonlinear iterations. */
    std::function<void( vector_type const& )> M_updateExternalState;
    /** @brief True once affine storage has been allocated and assembled. */
    bool M_assembled = false;
};

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_NONLINEAR_HPP
