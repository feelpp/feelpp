/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_LINEAR_HPP
#define FEELPP_FEELSYSTEM_LINEAR_HPP 1

/**
 * @file linear.hpp
 * @brief Assembled linear algebraic systems built from `Feel::systems::Problem`.
 */

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <feel/feelalg/backend.hpp>
#include <feel/feelsystem/auxiliary.hpp>
#include <feel/feelsystem/backend.hpp>
#include <feel/feelsystem/concepts.hpp>
#include <feel/feelsystem/namedarguments.hpp>
#include <feel/feelsystem/result.hpp>
#include <feel/feelsystem/terms.hpp>

namespace Feel::systems
{

/**
 * @brief Stateful assembled form of a stationary linear problem.
 *
 * `LinearSystem` owns the backend, sparse matrix, rhs vector, and active
 * constrained algebraic pair produced by a `LinearProblemLike` descriptor. The
 * descriptor remains immutable and can be assembled into multiple systems with
 * different backend prefixes or explicit backend instances.
 *
 * @tparam ProblemType Problem descriptor satisfying `LinearProblemLike`.
 */
template <LinearProblemLike ProblemType>
class LinearSystem
{
  public:
    /** @brief Decayed problem descriptor type. */
    using problem_type = std::remove_cvref_t<ProblemType>;
    /** @brief Function-space type used by the problem. */
    using space_type = typename problem_type::space_type;
    /** @brief Scalar value type used by algebraic objects. */
    using value_type = typename problem_type::value_type;
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
    /** @brief Per-rank auxiliary operator result. */
    using auxiliary_operator_type = AssembledAuxiliaryOperator<value_type>;

    /**
     * @brief Build a linear system from a reusable problem descriptor.
     *
     * @param problem Problem descriptor whose terms are replayed by `assemble()`.
     * @param backend Optional backend to reuse. When empty, `assemble()` creates
     * a backend from `backendName`, the problem name, or `"linear_system"`.
     * @param backendName Optional backend prefix used when `backend` is empty.
     */
    LinearSystem( problem_type const& problem, backend_ptrtype backend = {}, std::string backendName = {} )
        : M_problem( problem ),
          M_backend( std::move( backend ) ),
          M_backendName( detail::defaultBackendName( std::move( backendName ), M_problem.name(), "linear_system" ) )
    {
        if ( M_backend )
            M_backendName = M_backend->prefix();
    }

    /**
     * @brief Allocate storage and assemble matrix, rhs, and constraints.
     *
     * The base matrix and vector are available through `matrixPtr()` and
     * `rhsPtr()`. The active pair returned by the deferred constraint machinery
     * is available through `activeSystem()`.
     *
     * @return Reference to this system for fluent use.
     */
    LinearSystem& assemble()
    {
        if ( !M_backend )
            M_backend = detail::makeBackend<value_type>( M_backendName, M_problem.space()->worldCommPtr() );
        if ( !M_matrix )
            M_matrix = M_backend->newMatrix( _test = M_problem.space(), _trial = M_problem.space() );
        M_matrix->zero();
        if ( !M_rhs )
            M_rhs = M_backend->newVector( _test = M_problem.space() );
        M_rhs->zero();

        auto activeSystem = M_problem.assembleLinear( M_matrix, M_rhs );
        M_matrix->close();
        M_rhs->close();
        activeSystem = detail::completeActiveSystem( std::move( activeSystem ), M_matrix, M_rhs );
        M_activeMatrix = activeSystem.first;
        M_activeRhs = activeSystem.second;
        if ( M_activeMatrix )
            M_activeMatrix->close();
        if ( M_activeRhs )
            M_activeRhs->close();
        this->assembleAuxiliaryOperators();
        M_assembled = true;
        return *this;
    }

    /** @brief Return true after `assemble()` has completed. */
    bool assembled() const noexcept { return M_assembled; }
    /** @brief Problem descriptor copied into this system. */
    problem_type const& problem() const noexcept { return M_problem; }
    /** @brief Backend used for allocation and solving. */
    backend_ptrtype const& backend() const noexcept { return M_backend; }
    /** @brief Backend prefix used for options and default backend lookup. */
    std::string const& backendName() const noexcept { return M_backendName; }

    /** @brief Assembled base matrix before deferred constraint replacement. */
    matrix_ptrtype const& matrixPtr() const noexcept { return M_matrix; }
    /** @brief Assembled base rhs before deferred constraint replacement. */
    vector_ptrtype const& rhsPtr() const noexcept { return M_rhs; }
    /** @brief Dereferenced assembled base matrix. */
    matrix_type const& matrix() const noexcept { return *M_matrix; }
    /** @brief Dereferenced assembled base rhs. */
    vector_type const& rhs() const noexcept { return *M_rhs; }

    /** @brief Matrix/vector pair actually used by the backend solve. */
    active_system_type activeSystem() const noexcept { return { M_activeMatrix, M_activeRhs }; }
    /** @brief Active matrix pointer used by the backend solve. */
    matrix_ptrtype activeMatrixPtr() const noexcept { return M_activeMatrix; }
    /** @brief Active rhs vector pointer used by the backend solve. */
    vector_ptrtype activeVectorPtr() const noexcept { return M_activeRhs; }
    /** @brief Return true when constraints replaced the base matrix or rhs. */
    bool constrained() const noexcept { return M_activeMatrix != M_matrix || M_activeRhs != M_rhs; }

    /** @brief Assembled auxiliary operators in declaration order. */
    std::vector<auxiliary_operator_type> const& auxiliaryOperators() const noexcept
    {
        return M_auxiliaryOperators;
    }

    /** @brief Return true when a named auxiliary operator has been assembled. */
    bool hasAuxiliaryOperator( std::string const& name ) const noexcept
    {
        return std::any_of( M_auxiliaryOperators.begin(), M_auxiliaryOperators.end(),
                            [&name]( auto const& op ) { return op.name() == name; } );
    }

    /** @brief Look up a named assembled auxiliary operator. */
    auxiliary_operator_type const& auxiliaryOperator( std::string const& name ) const
    {
        auto const it = std::find_if( M_auxiliaryOperators.begin(), M_auxiliaryOperators.end(),
                                      [&name]( auto const& op ) { return op.name() == name; } );
        if ( it == M_auxiliaryOperators.end() )
            throw std::out_of_range( "Feel::systems::LinearSystem auxiliary operator not found: " + name );
        return *it;
    }

    /**
     * @brief Solve the assembled linear system.
     *
     * Calls `assemble()` lazily when needed, then forwards all named arguments
     * to `Backend::solve()`.
     *
     * @param args Backend solve named arguments such as `_solution`.
     * @return Convergence, iteration, residual, and message summary.
     */
    template <typename... Args>
    solve_result solve( Args&&... args )
    {
        if ( !M_assembled )
            this->assemble();

        auto result = this->solveWithAuxiliaryOperators( std::forward<Args>( args )... );
        return solve_result{ result.isConverged(),
                             static_cast<int>( result.nIterations() ),
                             static_cast<double>( result.residual() ),
                             result.isConverged()
                                 ? std::string{ "feelsystem LinearSystem solve converged" }
                                 : std::string{ "feelsystem LinearSystem solve did not converge" } };
    }

  private:
    template <typename... Args>
    auto solveWithAuxiliaryOperators( Args&&... args )
    {
        std::string pc = M_backend->pcType();
        if constexpr ( detail::hasNamedArgument<na::pc, Args...> )
            pc = detail::namedArgumentValue<na::pc>( args... );

        if ( !this->hasAuxiliaryOperator( "geneo" ) &&
             this->needsInheritedLocalOperator( pc ) )
            this->assembleInheritedLocalOperator();

        if ( M_auxiliaryOperators.empty() )
            return M_backend->solve( _matrix = M_activeMatrix,
                                     _rhs = M_activeRhs,
                                     std::forward<Args>( args )... );

        if constexpr ( detail::hasNamedArgument<na::prec, Args...> )
        {
            auto preconditioner = detail::namedArgumentValue<na::prec>( args... );
            detail::attachAuxiliaryOperators( preconditioner, M_auxiliaryOperators );
            return M_backend->solve( _matrix = M_activeMatrix,
                                     _rhs = M_activeRhs,
                                     std::forward<Args>( args )... );
        }
        else
        {
            std::string factorPackage = M_backend->pcFactorMatSolverPackageType();
            if constexpr ( detail::hasNamedArgument<na::pcfactormatsolverpackage, Args...> )
                factorPackage = detail::namedArgumentValue<na::pcfactormatsolverpackage>( args... );

            auto preconditioner = detail::makeAuxiliaryPreconditioner(
                M_backend, M_activeMatrix, pc, factorPackage );
            detail::attachAuxiliaryOperators( preconditioner, M_auxiliaryOperators );
            return M_backend->solve( _matrix = M_activeMatrix,
                                     _rhs = M_activeRhs,
                                     _prec = preconditioner,
                                     std::forward<Args>( args )... );
        }
    }

    bool needsInheritedLocalOperator( std::string const& pc ) const
    {
        if ( pc != "hpddm" || !M_problem.hasInferredLocalOperator() )
            return false;

        // HPDDM-specific options are parsed lazily by ConfigurePCHPDDM and are
        // deliberately not part of the generic backend variables map. Mirror
        // that parsing here so the inherited operator is assembled only when
        // the HPDDM Neumann-matrix path has actually been requested.
        auto const optionName = prefixvm( M_backendName, "pc-hpddm-has-neumann" );
        po::options_description options( "FeelSystem inherited HPDDM operator option" );
        options.add_options()
            ( optionName.c_str(), po::value<bool>()->default_value( false ),
              "assemble and attach the inherited local problem operator" );

        po::variables_map vm;
        auto parser = Environment::commandLineParser();
        auto parsed = parser.options( options )
                          .style( po::command_line_style::allow_long |
                                  po::command_line_style::long_allow_adjacent |
                                  po::command_line_style::long_allow_next )
                          .allow_unregistered()
                          .run();
        po::store( parsed, vm );
        for ( auto& configFile : Environment::configFiles() )
            po::store( po::parse_config_file( std::get<1>( configFile ), options, true ), vm );
        if ( !vm.count( optionName ) )
        {
            auto const optionIt = Environment::vm().find( optionName );
            if ( optionIt != Environment::vm().end() )
                vm.insert( *optionIt );
        }
        po::notify( vm );
        return vm[optionName].template as<bool>();
    }

    void assembleInheritedLocalOperator()
    {
        auto layout = detail::makeAuxiliaryDataMap( M_matrix->mapRowPtr() );
        auto localBackend = detail::makeBackendLike<value_type>(
            M_backend, M_backendName + ".auxiliary.geneo",
            Environment::worldCommSeqPtr() );
        auto localMatrix = localBackend->newMatrix( layout.localMap, layout.localMap,
                                                    NON_HERMITIAN, false );
        auto const n = layout.localMap->nDof();
        localMatrix->init( n, n, n, n, 64, 0 );
        localMatrix->zero();
        auto localVector = localBackend->newVector( layout.localMap );
        localVector->zero();
        auto active = M_problem.assembleInferredLocalOperator( localMatrix, localVector );
        localMatrix = active.first;
        localMatrix->close();

        AuxiliaryOperatorDescriptor descriptor{
            "geneo", AuxiliaryOperatorKind::local_neumann,
            AuxiliaryOperatorSource::inherited_problem };
        M_auxiliaryOperators.emplace_back( std::move( descriptor ),
                                           std::move( localBackend ),
                                           std::move( localMatrix ),
                                           std::move( layout.globalDofIds ) );
    }

    void assembleAuxiliaryOperators()
    {
        M_auxiliaryOperators.clear();
        for ( auto const& descriptor : M_problem.auxiliaryOperatorDescriptors() )
        {
            auto layout = detail::makeAuxiliaryDataMap( M_matrix->mapRowPtr() );
            auto localBackend = detail::makeBackendLike<value_type>(
                M_backend,
                M_backendName + ".auxiliary." + descriptor.name,
                Environment::worldCommSeqPtr() );
            auto localMatrix = localBackend->newMatrix( layout.localMap, layout.localMap,
                                                        NON_HERMITIAN, false );
            auto const n = layout.localMap->nDof();
            localMatrix->init( n, n, n, n, 64, 0 );
            localMatrix->zero();
            M_problem.assembleAuxiliaryOperator( descriptor.name, localMatrix );
            localMatrix->close();
            M_auxiliaryOperators.emplace_back( descriptor,
                                               std::move( localBackend ),
                                               std::move( localMatrix ),
                                               std::move( layout.globalDofIds ) );
        }
    }

    /** @brief Problem descriptor copied by value for replayable assembly. */
    problem_type M_problem;
    /** @brief Backend used for allocation and solving. */
    backend_ptrtype M_backend;
    /** @brief Backend prefix selected from user input, backend, or problem name. */
    std::string M_backendName;
    /** @brief Owned assembled base matrix. */
    matrix_ptrtype M_matrix;
    /** @brief Owned assembled base rhs vector. */
    vector_ptrtype M_rhs;
    /** @brief Matrix passed to the solver after constraints are applied. */
    matrix_ptrtype M_activeMatrix;
    /** @brief Rhs passed to the solver after constraints are applied. */
    vector_ptrtype M_activeRhs;
    /** @brief True once storage has been assembled at least once. */
    bool M_assembled = false;
    /** @brief Sequential per-rank auxiliary matrices and overlap mappings. */
    std::vector<auxiliary_operator_type> M_auxiliaryOperators;
};

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_LINEAR_HPP
