/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_STRATEGY_HPP
#define FEELPP_FEELSYSTEM_STRATEGY_HPP 1

/**
 * @file strategy.hpp
 * @brief Solver strategy records and factories for FeelSystem solve facades.
 */

#include <utility>

#include <feel/feelalg/enums.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelsystem/namedarguments.hpp>

namespace Feel::systems
{

/**
 * @brief Marker strategy for backend Newton/SNES nonlinear solves.
 */
struct newton_strategy
{};

/**
 * @brief Build a Newton strategy marker.
 */
inline newton_strategy
newton()
{
    return {};
}

/**
 * @brief Parameters for Picard fixed-point iterations.
 */
struct picard_strategy
{
    /** @brief Maximum number of fixed-point iterations. */
    int maxit = 25;
    /** @brief Stopping tolerance on the Picard update norm. */
    double rtolerance = 1e-8;
    /** @brief Relaxation factor applied to the update, where `1` is unrelaxed. */
    double relaxation = 1.0;
};

/**
 * @brief Select eigenvalues with smallest magnitude.
 */
struct smallestMagnitude {};

/**
 * @brief Select eigenvalues with largest magnitude.
 */
struct largestMagnitude {};

/**
 * @brief Select eigenvalues with smallest real part.
 */
struct smallestReal {};

/**
 * @brief Select eigenvalues with largest real part.
 */
struct largestReal {};

/**
 * @brief Select a shift-invert spectral transform.
 */
struct shiftInvert {};

/**
 * @brief Select a plain shift spectral transform.
 */
struct shiftTransform {};

namespace detail
{
/** @brief Convert a smallest-magnitude tag to the backend enum. */
inline PositionOfSpectrum eigenSpectrum( smallestMagnitude ) noexcept { return SMALLEST_MAGNITUDE; }
/** @brief Convert a largest-magnitude tag to the backend enum. */
inline PositionOfSpectrum eigenSpectrum( largestMagnitude ) noexcept { return LARGEST_MAGNITUDE; }
/** @brief Convert a smallest-real tag to the backend enum. */
inline PositionOfSpectrum eigenSpectrum( smallestReal ) noexcept { return SMALLEST_REAL; }
/** @brief Convert a largest-real tag to the backend enum. */
inline PositionOfSpectrum eigenSpectrum( largestReal ) noexcept { return LARGEST_REAL; }
/** @brief Forward an already concrete spectrum enum. */
inline PositionOfSpectrum eigenSpectrum( PositionOfSpectrum spectrum ) noexcept { return spectrum; }

/** @brief Convert a shift-invert tag to the backend enum. */
inline SpectralTransformType eigenTransform( shiftInvert ) noexcept { return SINVERT; }
/** @brief Convert a shift-transform tag to the backend enum. */
inline SpectralTransformType eigenTransform( shiftTransform ) noexcept { return SHIFT; }
/** @brief Forward an already concrete transform enum. */
inline SpectralTransformType eigenTransform( SpectralTransformType transform ) noexcept { return transform; }
} // namespace detail

/**
 * @brief Solver intent for standard and generalized eigenvalue problems.
 */
struct eigen_strategy
{
    /** @brief Number of requested eigenvalues, or negative to use options/defaults. */
    int nev = -1;
    /** @brief Number of basis vectors, or negative to use options/defaults. */
    int ncv = -1;
    /** @brief Maximum projected dimension, or negative to use options/defaults. */
    int mpd = -2;
    /** @brief Maximum eigensolver iterations, or negative to use options/defaults. */
    int maxit = -1;
    /** @brief Eigensolver tolerance, or negative to use options/defaults. */
    double tolerance = -1.0;
    /** @brief Lower interval bound for interval-based eigensolver modes. */
    double intervalA = 0.0;
    /** @brief Upper interval bound for interval-based eigensolver modes. */
    double intervalB = 0.0;
    /** @brief Backend eigensolver algorithm. */
    EigenSolverType solver = KRYLOVSCHUR;
    /** @brief Standard or generalized eigenproblem type override. */
    EigenProblemType problem = INVALID_EIGENPROBLEMTYPE;
    /** @brief Requested part of the spectrum. */
    PositionOfSpectrum spectrum = INVALID_Postion_of_Spectrum;
    /** @brief Spectral transform used by the eigensolver. */
    SpectralTransformType transform = SHIFT;
    /** @brief True when `transform` was explicitly set by named arguments. */
    bool hasTransform = false;
    /** @brief True when backend eigensolver verbosity should be enabled. */
    bool verbose = false;
};

/**
 * @brief Build Picard iteration parameters from Feel++ named arguments.
 *
 * Supported arguments are `_maxit`, `_rtolerance`, and `_initial_theta`
 * for the relaxation factor.
 */
template <typename... Args>
picard_strategy
picard( Args&&... args )
{
    auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
    return picard_strategy{
        static_cast<int>( namedArgs.get_else( _maxit, 25 ) ),
        static_cast<double>( namedArgs.get_else( _rtolerance, 1e-8 ) ),
        static_cast<double>( namedArgs.get_else( _initial_theta, 1.0 ) )
    };
}

/**
 * @brief Build eigen solver intent from Feel++ named arguments.
 *
 * Supported arguments are `_nev`, `_ncv`, `_mpd`, `_maxit`, `_tolerance`,
 * `_solver`, `_problem`, `_spectrum`, `_transform`, `_interval_a`,
 * `_interval_b`, and `_verbose`.
 */
template <typename... Args>
eigen_strategy
eigen( Args&&... args )
{
    auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
    auto strategy = eigen_strategy{};
    strategy.nev = static_cast<int>( namedArgs.get_else( _nev, -1 ) );
    strategy.ncv = static_cast<int>( namedArgs.get_else( _ncv, -1 ) );
    strategy.mpd = static_cast<int>( namedArgs.get_else( _mpd, -2 ) );
    strategy.maxit = static_cast<int>( namedArgs.get_else( _maxit, -1 ) );
    strategy.tolerance = static_cast<double>( namedArgs.get_else( _tolerance, -1.0 ) );
    strategy.intervalA = static_cast<double>( namedArgs.get_else( _interval_a, 0.0 ) );
    strategy.intervalB = static_cast<double>( namedArgs.get_else( _interval_b, 0.0 ) );
    strategy.solver = namedArgs.get_else( _solver, KRYLOVSCHUR );
    strategy.problem = namedArgs.get_else( _problem, INVALID_EIGENPROBLEMTYPE );
    strategy.spectrum = namedArgs.get_else( _spectrum, INVALID_Postion_of_Spectrum );
    strategy.transform = namedArgs.get_else( _transform, SHIFT );
    strategy.hasTransform = detail::hasNamedArgument<na::transform, Args...>;
    strategy.verbose = static_cast<bool>( namedArgs.get_else( _verbose, false ) );
    return strategy;
}

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_STRATEGY_HPP
