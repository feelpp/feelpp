/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_EIGEN_HPP
#define FEELPP_FEELSYSTEM_EIGEN_HPP 1

/**
 * @file eigen.hpp
 * @brief Assembled standard and generalized eigenvalue systems.
 */

#include <algorithm>
#include <cstddef>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelsystem/backend.hpp>
#include <feel/feelsystem/concepts.hpp>
#include <feel/feelsystem/namedarguments.hpp>
#include <feel/feelsystem/result.hpp>
#include <feel/feelsystem/strategy.hpp>

namespace Feel::systems
{

/**
 * @brief Stateful assembled form of a stationary eigenvalue problem.
 *
 * The system assembles `operatorA()` and, when present, `operatorB()` from a
 * `GeneralizedEigenvalueProblemLike` descriptor and delegates the eigensolve to
 * the existing Feel++ SLEPc backend.
 *
 * @tparam ProblemType Problem descriptor carrying `operatorA` and optional
 * `operatorB` terms.
 */
template <GeneralizedEigenvalueProblemLike ProblemType>
class EigenSystem
{
public:
    /** @brief Decayed problem descriptor type. */
    using problem_type = std::remove_cvref_t<ProblemType>;
    /** @brief Scalar value type used by the eigensolver. */
    using value_type = typename problem_type::value_type;
    /** @brief Finite-element eigenmode type. */
    using element_type = typename problem_type::element_type;
    /** @brief Feel++ backend type used for matrix allocation. */
    using backend_type = Backend<value_type>;
    /** @brief Shared backend pointer type. */
    using backend_ptrtype = std::shared_ptr<backend_type>;
    /** @brief Sparse matrix pointer type allocated by the backend. */
    using matrix_ptrtype = typename backend_type::sparse_matrix_ptrtype;
    /** @brief Vector pointer type allocated by the backend. */
    using vector_ptrtype = typename backend_type::vector_ptrtype;
    /** @brief Stored eigenpair with real eigenvalue part and eigenmode field. */
    using eigenmode_type = std::pair<value_type, element_type>;
    /** @brief Collection of stored eigenmodes. */
    using eigenmodes_type = std::vector<eigenmode_type>;

    /**
     * @brief Build an eigen system from a reusable problem descriptor.
     */
    EigenSystem( problem_type const& problem, backend_ptrtype backend = {}, std::string backendName = {} )
        : M_problem( problem ),
          M_backend( std::move( backend ) ),
          M_backendName( detail::defaultBackendName( std::move( backendName ), M_problem.name(), "eigen_system" ) )
    {
        if ( M_backend )
            M_backendName = M_backend->prefix();
    }

    /**
     * @brief Allocate and assemble operator matrices.
     *
     * `operatorB()` is allocated only when the problem has `operator_b` terms.
     *
     * @return Reference to this system for fluent use.
     */
    EigenSystem& assemble()
    {
        if ( !M_backend )
            M_backend = detail::makeBackend<value_type>( M_backendName, M_problem.space()->worldCommPtr() );
        if ( !M_A )
            M_A = M_backend->newMatrix( _test = M_problem.space(), _trial = M_problem.space() );
        M_A->zero();
        M_hasOperatorB = !M_problem.operatorBTerms().empty();
        if ( M_hasOperatorB )
        {
            if ( !M_B )
                M_B = M_backend->newMatrix( _test = M_problem.space(), _trial = M_problem.space() );
            M_B->zero();
        }
        else
        {
            M_B.reset();
        }

        M_problem.assembleEigen( M_A, M_B );
        M_A->close();
        if ( M_B )
            M_B->close();
        M_assembled = true;
        return *this;
    }

    /** @brief Return true once operator matrices have been assembled. */
    bool assembled() const noexcept { return M_assembled; }
    /** @brief Problem descriptor copied into this system. */
    problem_type const& problem() const noexcept { return M_problem; }
    /** @brief Backend used for matrix allocation. */
    backend_ptrtype const& backend() const noexcept { return M_backend; }
    /** @brief Backend/eigensolver option prefix. */
    std::string const& backendName() const noexcept { return M_backendName; }

    /** @brief Assembled left/operator A matrix. */
    matrix_ptrtype const& operatorA() const noexcept { return M_A; }
    /** @brief Assembled right/operator B matrix, or empty for standard problems. */
    matrix_ptrtype const& operatorB() const noexcept { return M_B; }
    /** @brief Return true when `operatorB()` participates in the solve. */
    bool generalized() const noexcept { return M_hasOperatorB; }
    /** @brief Real parts of the eigenvalues from the last solve. */
    std::vector<value_type> const& eigenvalues() const noexcept { return M_eigenvalues; }
    /** @brief Imaginary parts of the eigenvalues from the last solve. */
    std::vector<value_type> const& eigenvaluesImaginary() const noexcept { return M_eigenvaluesImaginary; }
    /** @brief Eigenmodes copied into finite-element fields by the last solve. */
    eigenmodes_type const& eigenModes() const noexcept { return M_modes; }

    /**
     * @brief Solve the assembled eigenvalue problem.
     *
     * Named arguments may include `_strategy = systems::eigen(...)`; unset
     * fields fall back to problem-level settings and backend options.
     *
     * @return Eigensolver convergence information and stored eigenvalues.
     */
    template <typename... Args>
    eigen_result solve( Args&&... args )
    {
        if ( !M_assembled )
            this->assemble();

        auto strategy = this->makeStrategy( std::forward<Args>( args )... );

        try
        {
            auto solver = SolverEigen<value_type>::build( SOLVERS_SLEPC );
            solver->setEigenSolverType( strategy.solver );
            solver->setEigenProblemType( this->problemType( strategy ) );
            solver->setPositionOfSpectrum( strategy.spectrum );
            solver->setSpectralTransform( strategy.transform );
            solver->setNumberOfEigenValues( strategy.nev );
            solver->setNumberOfEigenValuesConverged( strategy.ncv );
            solver->setMaximumProjectedDimension( strategy.mpd > 0 ? strategy.mpd : invalid_v<typename SolverEigen<value_type>::size_type> );
            solver->setInterval( strategy.intervalA, strategy.intervalB );
            solver->setMaxIterations( strategy.maxit );
            solver->setTolerance( strategy.tolerance );
            solver->setMapRow( M_A->mapRowPtr() );
            solver->setMapCol( M_A->mapColPtr() );

            typename SolverEigen<value_type>::solve_return_type rawResult;
            if ( M_hasOperatorB )
                rawResult = solver->solve( M_A, M_B, strategy.nev, strategy.ncv,
                                           strategy.tolerance, strategy.maxit );
            else
                rawResult = solver->solve( M_A, strategy.nev, strategy.ncv,
                                           strategy.tolerance, strategy.maxit );

            int converged = 0;
            int iterations = 0;
            std::vector<double> errors;
            boost::tie( converged, iterations, errors ) = rawResult;
            this->storeModes( solver->eigenModes() );

            auto result = eigen_result{};
            result.converged = converged;
            result.iterations = iterations;
            result.errors = std::move( errors );
            result.eigenvalues = M_eigenvalues;
            result.eigenvaluesImaginary = M_eigenvaluesImaginary;
            std::ostringstream os;
            os << "feelsystem EigenSystem solve converged " << converged
               << " eigenpair(s) in " << iterations << " iteration(s)";
            result.message = os.str();
            return result;
        }
        catch ( std::exception const& e )
        {
            M_modes.clear();
            M_eigenvalues.clear();
            M_eigenvaluesImaginary.clear();
            return eigen_result{ 0, 0, {}, {}, {},
                                 std::string{ "feelsystem EigenSystem solve unavailable or failed: " } + e.what() };
        }
    }

private:
    /** @brief Compose a backend-prefixed option key. */
    static std::string optionName( std::string prefix, std::string const& name )
    {
        if ( !prefix.empty() && prefix.back() != '-' )
            prefix += "-";
        return prefix + name;
    }

    /** @brief Read an integer eigensolver option with backend-prefix fallback. */
    int optionInt( std::string const& name, int fallback ) const
    {
        auto key = optionName( M_backendName, name );
        auto const& vm = Environment::vm();
        if ( vm.count( key ) )
            return vm[key].template as<int>();
        if ( vm.count( name ) )
            return vm[name].template as<int>();
        return fallback;
    }

    /** @brief Read a floating-point eigensolver option with backend-prefix fallback. */
    double optionDouble( std::string const& name, double fallback ) const
    {
        auto key = optionName( M_backendName, name );
        auto const& vm = Environment::vm();
        if ( vm.count( key ) )
            return vm[key].template as<double>();
        if ( vm.count( name ) )
            return vm[name].template as<double>();
        return fallback;
    }

    /** @brief Read a boolean eigensolver option with backend-prefix fallback. */
    bool optionBool( std::string const& name, bool fallback ) const
    {
        auto key = optionName( M_backendName, name );
        auto const& vm = Environment::vm();
        if ( vm.count( key ) )
            return vm[key].template as<bool>();
        if ( vm.count( name ) )
            return vm[name].template as<bool>();
        return fallback;
    }

    /** @brief Read and map an enum-valued eigensolver option. */
    template <typename EnumType>
    EnumType optionEnum( std::string const& name, EnumType fallback ) const
    {
        auto key = optionName( M_backendName, name );
        auto const& vm = Environment::vm();
        std::string value;
        if ( vm.count( key ) )
            value = vm[key].template as<std::string>();
        else if ( vm.count( name ) )
            value = vm[name].template as<std::string>();
        else
            return fallback;

        auto it = EigenMap.find( value );
        if ( it == EigenMap.end() )
            return fallback;
        return static_cast<EnumType>( it->second );
    }

    /** @brief Merge explicit strategy, problem defaults, and backend options. */
    template <typename... Args>
    eigen_strategy makeStrategy( Args&&... args ) const
    {
        auto strategy = eigen_strategy{};
        if constexpr ( detail::hasNamedArgument<na::strategy, Args...> )
            strategy = detail::namedArgumentValue<na::strategy>( std::forward<Args>( args )... );

        strategy.solver = this->optionEnum( "solvereigen.solver", strategy.solver );
        strategy.nev = strategy.nev > 0
                           ? strategy.nev
                           : ( M_problem.hasNev()
                                   ? M_problem.nev()
                                   : this->optionInt( "solvereigen.nev", M_problem.nev() ) );
        strategy.ncv = strategy.ncv > 0 ? strategy.ncv : M_problem.ncv();
        strategy.mpd = strategy.mpd > 0 ? strategy.mpd : M_problem.mpd();
        strategy.maxit = strategy.maxit > 0 ? strategy.maxit : M_problem.maxit();
        strategy.tolerance = strategy.tolerance > 0.0 ? strategy.tolerance : M_problem.tolerance();
        if ( strategy.ncv <= 0 )
            strategy.ncv = this->optionInt( "solvereigen.ncv", std::max( 2 * strategy.nev, strategy.nev + 2 ) );
        strategy.ncv = std::max( strategy.ncv, strategy.nev + 1 );
        if ( strategy.mpd <= 0 )
            strategy.mpd = this->optionInt( "solvereigen.mpd", -2 );
        if ( strategy.maxit <= 0 )
            strategy.maxit = this->optionInt( "solvereigen.maxiter", 10000 );
        if ( strategy.tolerance <= 0.0 )
            strategy.tolerance = this->optionDouble( "solvereigen.tolerance", 1e-10 );
        if ( strategy.spectrum == INVALID_Postion_of_Spectrum )
            strategy.spectrum = M_problem.hasSpectrum()
                                    ? M_problem.spectrum()
                                    : this->optionEnum( "solvereigen.spectrum", LARGEST_MAGNITUDE );
        if ( !strategy.hasTransform && M_problem.hasTransform() )
            strategy.transform = M_problem.transform();
        else if ( !strategy.hasTransform )
            strategy.transform = this->optionEnum( "solvereigen.transform", SHIFT );
        strategy.intervalA = this->optionDouble( "solvereigen.interval-a", strategy.intervalA );
        strategy.intervalB = this->optionDouble( "solvereigen.interval-b", strategy.intervalB );
        strategy.verbose = this->optionBool( "solvereigen.verbose", strategy.verbose );
        return strategy;
    }

    /** @brief Select the backend eigenproblem type from strategy/problem shape. */
    EigenProblemType problemType( eigen_strategy const& strategy ) const
    {
        if ( strategy.problem != INVALID_EIGENPROBLEMTYPE )
            return strategy.problem;
        return M_hasOperatorB ? GHEP : HEP;
    }

    /** @brief Copy raw eigensolver modes into finite-element fields and value arrays. */
    void storeModes( typename SolverEigen<value_type>::eigenmodes_type const& rawModes )
    {
        M_modes.clear();
        M_eigenvalues.clear();
        M_eigenvaluesImaginary.clear();
        M_modes.reserve( rawModes.size() );
        M_eigenvalues.reserve( rawModes.size() );
        M_eigenvaluesImaginary.reserve( rawModes.size() );

        for ( auto const& mode : rawModes )
        {
            auto eigenmode = M_problem.space()->element();
            eigenmode = *mode.second.template get<2>();
            eigenmode.close();
            M_eigenvalues.push_back( mode.second.template get<0>() );
            M_eigenvaluesImaginary.push_back( mode.second.template get<1>() );
            M_modes.emplace_back( mode.second.template get<0>(), std::move( eigenmode ) );
        }
    }

    /** @brief Problem descriptor copied by value for replayable assembly. */
    problem_type M_problem;
    /** @brief Backend used to allocate operator matrices. */
    backend_ptrtype M_backend;
    /** @brief Backend/eigensolver option prefix. */
    std::string M_backendName;
    /** @brief Assembled operator A matrix. */
    matrix_ptrtype M_A;
    /** @brief Assembled operator B matrix for generalized solves. */
    matrix_ptrtype M_B;
    /** @brief Finite-element eigenmodes from the last solve. */
    eigenmodes_type M_modes;
    /** @brief Real eigenvalue parts from the last solve. */
    std::vector<value_type> M_eigenvalues;
    /** @brief Imaginary eigenvalue parts from the last solve. */
    std::vector<value_type> M_eigenvaluesImaginary;
    /** @brief True when the current problem has an operator B. */
    bool M_hasOperatorB = false;
    /** @brief True after `assemble()` has completed. */
    bool M_assembled = false;
};

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_EIGEN_HPP
