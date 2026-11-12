/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_PICARD_HPP
#define FEELPP_FEELSYSTEM_PICARD_HPP 1

/**
 * @file picard.hpp
 * @brief Picard fixed-point systems with one linear solve per iteration.
 */

#include <memory>
#include <string>
#include <utility>

#include <feel/feelalg/backend.hpp>
#include <feel/feelsystem/backend.hpp>
#include <feel/feelsystem/concepts.hpp>
#include <feel/feelsystem/namedarguments.hpp>
#include <feel/feelsystem/result.hpp>
#include <feel/feelsystem/strategy.hpp>

namespace Feel::systems
{

/**
 * @brief Fixed-point system assembling `A(u_k) u_{k+1} = b(u_k)`.
 *
 * Picard systems are nonlinear algorithms with a linear solve at each
 * iteration. The problem stores callbacks through `.iteration(...)`; the system
 * owns the backend, matrix, and rhs reused across iterations.
 */
template <PicardProblemLike ProblemType>
class PicardSystem
{
public:
    /** @brief Decayed problem descriptor type. */
    using problem_type = std::remove_cvref_t<ProblemType>;
    /** @brief Scalar value type used by algebraic objects. */
    using value_type = typename problem_type::value_type;
    /** @brief Finite-element state type. */
    using element_type = typename problem_type::element_type;
    /** @brief Feel++ backend type used for allocation and linear solves. */
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
     * @brief Build a Picard system from a problem descriptor.
     */
    PicardSystem( problem_type const& problem, backend_ptrtype backend = {}, std::string backendName = {} )
        : M_problem( problem ),
          M_backend( std::move( backend ) ),
          M_backendName( detail::defaultBackendName( std::move( backendName ), M_problem.name(), "picard_system" ) )
    {
        if ( M_backend )
            M_backendName = M_backend->prefix();
    }

    /**
     * @brief Allocate reusable matrix and vector storage.
     */
    PicardSystem& assemble()
    {
        if ( !M_backend )
            M_backend = detail::makeBackend<value_type>( M_backendName, M_problem.space()->worldCommPtr() );
        if ( !M_matrix )
            M_matrix = M_backend->newMatrix( _test = M_problem.space(), _trial = M_problem.space() );
        if ( !M_rhs )
            M_rhs = M_backend->newVector( _test = M_problem.space() );
        M_assembled = true;
        return *this;
    }

    /** @brief Return true once matrix/vector storage has been allocated. */
    bool assembled() const noexcept { return M_assembled; }
    /** @brief Problem descriptor copied into this system. */
    problem_type const& problem() const noexcept { return M_problem; }
    /** @brief Backend used for allocation and linear solves. */
    backend_ptrtype const& backend() const noexcept { return M_backend; }
    /** @brief Backend prefix used for options and default backend lookup. */
    std::string const& backendName() const noexcept { return M_backendName; }

    /** @brief Matrix assembled for the current Picard iteration. */
    matrix_ptrtype const& matrixPtr() const noexcept { return M_matrix; }
    /** @brief Rhs assembled for the current Picard iteration. */
    vector_ptrtype const& rhsPtr() const noexcept { return M_rhs; }
    /** @brief Dereferenced current Picard matrix. */
    matrix_type const& matrix() const noexcept { return *M_matrix; }
    /** @brief Dereferenced current Picard rhs. */
    vector_type const& rhs() const noexcept { return *M_rhs; }

    /**
     * @brief Solve the fixed-point problem with repeated linear solves.
     *
     * Pass `_strategy = systems::picard(...)` to configure `_maxit`,
     * `_rtolerance`, and `_initial_theta` relaxation.
     */
    template <typename... Args>
    solve_result solve( Args&&... args )
    {
        if ( !M_assembled )
            this->assemble();

        static_assert( detail::hasNamedArgument<na::solution, Args...>,
                       "Feel::systems::PicardSystem::solve requires _solution" );
        auto&& solution = detail::namedArgumentValue<na::solution>( args... );
        auto strategy = picard_strategy{};
        if constexpr ( detail::hasNamedArgument<na::strategy, Args...> )
            strategy = detail::namedArgumentValue<na::strategy>( args... );

        auto current = Feel::detail::ref( solution );
        auto next = current;
        auto delta = current;
        double residual = 0.0;
        bool converged = false;
        int iteration = 0;

        for ( ; iteration < strategy.maxit; ++iteration )
        {
            M_matrix->zero();
            M_rhs->zero();
            M_problem.assembleAffine( M_matrix, M_rhs );
            M_problem.assemblePicardIteration( current, M_matrix, M_rhs );
            M_matrix->close();
            M_rhs->close();

            next = current;
            auto linearResult = M_backend->solve( _matrix = M_matrix,
                                                  _rhs = M_rhs,
                                                  _solution = next );
            if ( !linearResult.isConverged() )
            {
                return solve_result{ false,
                                     iteration,
                                     static_cast<double>( linearResult.residual() ),
                                     "feelsystem PicardSystem inner linear solve did not converge" };
            }

            delta = next;
            delta -= current;
            residual = delta.l2Norm();

            if ( strategy.relaxation != 1.0 )
            {
                delta.scale( strategy.relaxation );
                next = current;
                next += delta;
            }

            Feel::detail::ref( solution ) = next;
            Feel::detail::ref( solution ).close();
            current = next;

            if ( residual <= strategy.rtolerance )
            {
                converged = true;
                ++iteration;
                break;
            }
        }

        return solve_result{ converged,
                             iteration,
                             residual,
                             converged
                                 ? std::string{ "feelsystem PicardSystem solve converged" }
                                 : std::string{ "feelsystem PicardSystem solve did not converge" } };
    }

private:
    /** @brief Problem descriptor copied by value for replayable iteration assembly. */
    problem_type M_problem;
    /** @brief Backend used for allocation and linear solves. */
    backend_ptrtype M_backend;
    /** @brief Backend prefix selected from user input, backend, or problem name. */
    std::string M_backendName;
    /** @brief Matrix reused and reassembled at every Picard iteration. */
    matrix_ptrtype M_matrix;
    /** @brief Rhs reused and reassembled at every Picard iteration. */
    vector_ptrtype M_rhs;
    /** @brief True once reusable matrix and rhs storage have been allocated. */
    bool M_assembled = false;
};

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_PICARD_HPP
