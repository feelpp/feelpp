/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_RESULT_HPP
#define FEELPP_FEELSYSTEM_RESULT_HPP 1

/**
 * @file result.hpp
 * @brief Lightweight result records returned by FeelSystem solve facades.
 */

#include <cstddef>
#include <string>
#include <vector>

namespace Feel::systems
{

/**
 * @brief Summary returned by linear, nonlinear, Picard, block, and HDG solves.
 */
struct solve_result
{
    /** @brief True when the backend or algorithm reported convergence. */
    bool converged = false;
    /** @brief Iteration count reported by the solver, or zero when unavailable. */
    int iterations = 0;
    /** @brief Final residual norm reported by the solver, or zero when unavailable. */
    double residual = 0.0;
    /** @brief Human-readable convergence or failure message. */
    std::string message;
};

/**
 * @brief Summary returned by eigenvalue solves.
 */
struct eigen_result
{
    /** @brief Number of converged eigenpairs reported by the eigensolver. */
    int converged = 0;
    /** @brief Iteration count reported by the eigensolver. */
    int iterations = 0;
    /** @brief Per-eigenpair error estimates. */
    std::vector<double> errors;
    /** @brief Real parts of computed eigenvalues. */
    std::vector<double> eigenvalues;
    /** @brief Imaginary parts of computed eigenvalues. */
    std::vector<double> eigenvaluesImaginary;
    /** @brief Human-readable convergence or failure message. */
    std::string message;
};

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_RESULT_HPP
