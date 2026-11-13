/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_EVOLUTION_HPP
#define FEELPP_FEELSYSTEM_EVOLUTION_HPP 1

/**
 * @file evolution.hpp
 * @brief Stateful time-step systems for first- and second-order evolution problems.
 */

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include <feel/feelalg/backend.hpp>
#include <feel/feelsystem/backend.hpp>
#include <feel/feelsystem/concepts.hpp>
#include <feel/feelsystem/namedarguments.hpp>
#include <feel/feelsystem/problem.hpp>
#include <feel/feelsystem/result.hpp>
#include <feel/feelsystem/terms.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelts/newmark.hpp>

namespace Feel::systems
{

/**
 * @brief Stateful algebraic system for one evolution problem time loop.
 *
 * `TimeStepSystem` keeps the `feelts` scheme attached to the problem and
 * assembles the current step as a standard linear algebraic system. BDF
 * first-order problems produce `alpha0*M + K`; Newmark second-order problems
 * produce `a0*M + a1*C + K`. Time history contributions are assembled through
 * matrix-vector products with the scheme polynomials.
 *
 * Restart and file-save behavior remain owned by the underlying time scheme.
 */
template <EvolutionProblemLike ProblemType>
class TimeStepSystem
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
    /** @brief Feel++ backend type used for allocation and solving. */
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
    /** @brief Matrix/vector pair actually passed to the solver. */
    using active_system_type = std::pair<matrix_ptrtype, vector_ptrtype>;

    /**
     * @brief Build a time-step system from an evolution problem descriptor.
     *
     * If `backend` is empty, assembly creates or reuses a backend whose prefix
     * is derived from `backendName` or from the problem name.
     */
    TimeStepSystem( problem_type const& problem, backend_ptrtype backend = {}, std::string backendName = {} )
        : M_problem( problem ),
          M_backend( std::move( backend ) ),
          M_backendName( detail::defaultBackendName( std::move( backendName ), M_problem.name(), "time_step_system" ) )
    {
        if ( M_backend )
            M_backendName = M_backend->prefix();
    }

    /**
     * @brief Store the initial condition used by the next `start()` call.
     */
    TimeStepSystem& initialize( element_type const& initialCondition )
    {
        M_initialCondition.emplace( initialCondition );
        return *this;
    }

    /**
     * @brief Start the attached time scheme.
     */
    double start()
    {
        M_assembled = false;
        M_started = true;
        if ( auto scheme = M_problem.bdfTimeScheme() )
        {
            if ( M_initialCondition )
                return scheme->start( *M_initialCondition );
            return scheme->start();
        }
        if ( auto scheme = M_problem.newmarkTimeScheme() )
        {
            if ( M_initialCondition )
                return scheme->start( *M_initialCondition );
            return scheme->start();
        }
        throw std::logic_error( "Feel::systems::TimeStepSystem requires a BDF or Newmark time scheme" );
    }

    /**
     * @brief Store an initial condition and start the attached time scheme.
     */
    double start( element_type const& initialCondition )
    {
        this->initialize( initialCondition );
        return this->start();
    }

    /**
     * @brief Advance the time scheme using the accepted solution.
     */
    double advance( element_type const& solution )
    {
        M_assembled = false;
        if ( auto scheme = M_problem.bdfTimeScheme() )
            return scheme->next( solution );
        if ( auto scheme = M_problem.newmarkTimeScheme() )
            return scheme->next( solution );
        throw std::logic_error( "Feel::systems::TimeStepSystem requires a BDF or Newmark time scheme" );
    }

    /**
     * @brief Alias returning the current step system for fluent loops.
     */
    TimeStepSystem& step() noexcept { return *this; }
    TimeStepSystem const& step() const noexcept { return *this; }

    /**
     * @brief Assemble the current time step.
     */
    TimeStepSystem& assemble()
    {
        this->ensureStorage();
        M_matrix->zero();
        M_rhs->zero();

        auto context = this->context();
        M_problem.updateForTimeStep( context );
        if constexpr ( FirstOrderEvolutionProblemLike<problem_type> )
        {
            if ( auto scheme = M_problem.bdfTimeScheme() )
                this->assembleBdfStep( context, *scheme );
            else
                throw std::logic_error( "Feel::systems first-order evolution requires a BDF time scheme" );
        }
        else if constexpr ( SecondOrderEvolutionProblemLike<problem_type> )
        {
            if ( auto scheme = M_problem.newmarkTimeScheme() )
                this->assembleNewmarkStep( context, *scheme );
            else
                throw std::logic_error( "Feel::systems second-order evolution requires a Newmark time scheme" );
        }

        auto activeSystem = M_problem.applyConstraints( context, M_matrix, M_rhs );
        M_matrix->close();
        M_rhs->close();
        activeSystem = detail::completeActiveSystem( std::move( activeSystem ), M_matrix, M_rhs );
        M_activeMatrix = activeSystem.first;
        M_activeRhs = activeSystem.second;
        if ( M_activeMatrix )
            M_activeMatrix->close();
        if ( M_activeRhs )
            M_activeRhs->close();
        M_assembled = true;
        return *this;
    }

    bool assembled() const noexcept { return M_assembled; }
    bool started() const noexcept { return M_started; }
    problem_type const& problem() const noexcept { return M_problem; }
    backend_ptrtype const& backend() const noexcept { return M_backend; }
    std::string const& backendName() const noexcept { return M_backendName; }

    matrix_ptrtype const& matrixPtr() const noexcept { return M_matrix; }
    vector_ptrtype const& rhsPtr() const noexcept { return M_rhs; }
    matrix_type const& matrix() const noexcept { return *M_matrix; }
    vector_type const& rhs() const noexcept { return *M_rhs; }

    active_system_type activeSystem() const noexcept { return { M_activeMatrix, M_activeRhs }; }
    matrix_ptrtype activeMatrixPtr() const noexcept { return M_activeMatrix; }
    vector_ptrtype activeVectorPtr() const noexcept { return M_activeRhs; }
    bool constrained() const noexcept { return M_activeMatrix != M_matrix || M_activeRhs != M_rhs; }

    /**
     * @brief Return true when the attached time scheme has reached its final time.
     */
    bool isFinished() const
    {
        if ( auto scheme = M_problem.bdfTimeScheme() )
            return scheme->isFinished();
        if ( auto scheme = M_problem.newmarkTimeScheme() )
            return scheme->isFinished();
        throw std::logic_error( "Feel::systems::TimeStepSystem requires a BDF or Newmark time scheme" );
    }

    /**
     * @brief Current target time for the assembled step.
     */
    double time() const { return this->context().time; }

    /**
     * @brief Current time step size.
     */
    double timeStep() const { return this->context().timeStep; }

    /**
     * @brief Current time iteration index.
     */
    int iteration() const { return this->context().iteration; }

    /**
     * @brief Metadata used by internal time-step assembly.
     */
    time_step_context context() const
    {
        if ( auto scheme = M_problem.bdfTimeScheme() )
        {
            return time_step_context{ time_scheme_kind::bdf,
                                      scheme->time(),
                                      scheme->timeStep(),
                                      scheme->iteration(),
                                      scheme->polyDerivCoefficient( 0 ),
                                      0.0,
                                      0.0 };
        }
        if ( auto scheme = M_problem.newmarkTimeScheme() )
        {
            return time_step_context{ time_scheme_kind::newmark,
                                      scheme->time(),
                                      scheme->timeStep(),
                                      scheme->iteration(),
                                      scheme->polyDerivCoefficient(),
                                      scheme->polyFirstDerivCoefficient(),
                                      scheme->polySecondDerivCoefficient() };
        }
        throw std::logic_error( "Feel::systems::TimeStepSystem requires a BDF or Newmark time scheme" );
    }

    /**
     * @brief Solve the assembled step with the configured backend.
     */
    template <typename... Args>
    solve_result solve( Args&&... args )
    {
        if ( !M_assembled )
            this->assemble();

        auto result = M_backend->solve( _matrix = M_activeMatrix,
                                        _rhs = M_activeRhs,
                                        std::forward<Args>( args )... );
        return solve_result{ result.isConverged(),
                             static_cast<int>( result.nIterations() ),
                             static_cast<double>( result.residual() ),
                             result.isConverged()
                                 ? std::string{ "feelsystem TimeStepSystem solve converged" }
                                 : std::string{ "feelsystem TimeStepSystem solve did not converge" } };
    }

  private:
    void ensureStorage()
    {
        if ( !M_backend )
            M_backend = detail::makeBackend<value_type>( M_backendName, M_problem.space()->worldCommPtr() );
        if ( !M_matrix )
            M_matrix = M_backend->newMatrix( _test = M_problem.space(), _trial = M_problem.space() );
        if ( !M_rhs )
            M_rhs = M_backend->newVector( _test = M_problem.space() );
    }

    void addMatrixRole( time_step_context const& context, term_role role, double scale )
    {
        auto roleMatrix = M_backend->newMatrix( _test = M_problem.space(), _trial = M_problem.space() );
        roleMatrix->zero();
        M_problem.assembleMatrixTerms( context, role, roleMatrix );
        roleMatrix->close();
        M_matrix->addMatrix( scale, roleMatrix );
    }

    void addVectorRole( time_step_context const& context, term_role role, double scale )
    {
        auto roleVector = M_backend->newVector( _test = M_problem.space() );
        roleVector->zero();
        M_problem.assembleVectorTerms( context, role, roleVector );
        roleVector->close();
        M_rhs->add( scale, *roleVector );
    }

    void addHistoryRole( time_step_context const& context, term_role role, element_type const& history, double scale )
    {
        auto roleMatrix = M_backend->newMatrix( _test = M_problem.space(), _trial = M_problem.space() );
        auto roleVector = M_backend->newVector( _test = M_problem.space() );
        roleMatrix->zero();
        M_problem.assembleMatrixTerms( context, role, roleMatrix );
        roleMatrix->close();
        roleVector->zero();
        roleMatrix->multVector( history, *roleVector );
        roleVector->close();
        M_rhs->add( scale, *roleVector );
    }

    /** @brief Assemble one BDF step from matrix roles and BDF history. */
    template <typename BdfType>
    void assembleBdfStep( time_step_context const& context, BdfType const& scheme )
    {
        this->addMatrixRole( context, term_role::lhs, 1.0 );
        this->addMatrixRole( context, term_role::stiffness, 1.0 );
        this->addMatrixRole( context, term_role::mass, context.derivativeCoefficient );
        this->addVectorRole( context, term_role::rhs, 1.0 );
        this->addHistoryRole( context, term_role::mass, scheme.polyDeriv(), 1.0 );
    }

    /** @brief Assemble one Newmark step from matrix roles and Newmark history. */
    template <typename NewmarkType>
    void assembleNewmarkStep( time_step_context const& context, NewmarkType const& scheme )
    {
        this->addMatrixRole( context, term_role::lhs, 1.0 );
        this->addMatrixRole( context, term_role::stiffness, 1.0 );
        this->addMatrixRole( context, term_role::damping, context.firstDerivativeCoefficient );
        this->addMatrixRole( context, term_role::mass, context.secondDerivativeCoefficient );
        this->addVectorRole( context, term_role::rhs, 1.0 );
        this->addHistoryRole( context, term_role::damping, scheme.polyFirstDeriv(), 1.0 );
        this->addHistoryRole( context, term_role::mass, scheme.polySecondDeriv(), 1.0 );
    }

    /** @brief Problem descriptor copied by value for replayable time-step assembly. */
    problem_type M_problem;
    /** @brief Backend used for allocation and solving. */
    backend_ptrtype M_backend;
    /** @brief Backend prefix selected from user input, backend, or problem name. */
    std::string M_backendName;
    /** @brief Owned step matrix. */
    matrix_ptrtype M_matrix;
    /** @brief Owned step rhs vector. */
    vector_ptrtype M_rhs;
    /** @brief Matrix passed to the solver after constraints are applied. */
    matrix_ptrtype M_activeMatrix;
    /** @brief Rhs passed to the solver after constraints are applied. */
    vector_ptrtype M_activeRhs;
    /** @brief Optional initial condition used by the next `start()` call. */
    std::optional<element_type> M_initialCondition;
    /** @brief True after the attached time scheme has been started. */
    bool M_started = false;
    /** @brief True when the current step has been assembled. */
    bool M_assembled = false;
};

/**
 * @brief Function object creating `TimeStepSystem` instances.
 */
struct evolve_t
{
    /**
     * @brief Create a stateful time-step system from an evolution problem.
     */
    template <EvolutionProblemLike ProblemType, typename... Args>
    auto operator()( ProblemType const& problem, Args&&... args ) const
    {
        using system_type = TimeStepSystem<ProblemType>;
        auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
        auto backend = namedArgs.get_else( _backend, typename system_type::backend_ptrtype{} );
        std::string backendName = namedArgs.get_else( _name, std::string{} );
        return system_type( problem, backend, std::move( backendName ) );
    }
};

/**
 * @brief Global FeelSystem evolution customization-point object.
 */
inline constexpr evolve_t evolve{};

/**
 * @brief Member implementation of `Problem::evolve()`.
 */
template <typename SpaceType, typename... Tags>
template <typename... Args>
auto Problem<SpaceType, Tags...>::evolve( Args&&... args ) const
{
    return systems::evolve( *this, std::forward<Args>( args )... );
}

} // namespace Feel::systems

namespace Feel
{
inline constexpr auto const& evolve = systems::evolve;
}

#endif // FEELPP_FEELSYSTEM_EVOLUTION_HPP
