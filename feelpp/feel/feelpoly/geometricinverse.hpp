/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file geometricinverse.hpp
 * @brief Structured Newton and Gauss--Newton geometric-map inversion.
 */
#if !defined( FEELPP_FEELPOLY_GEOMETRICINVERSE_HPP )
#define FEELPP_FEELPOLY_GEOMETRICINVERSE_HPP 1

#include <Eigen/Core>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include <string_view>
#include <utility>

namespace Feel::geometry
{
/** @enum InverseStatus
 * @brief Numerical termination status of a geometric inverse operation.
 */
enum class InverseStatus
{
    /** The physical residual satisfies the requested tolerance. */
    converged_exact,
    /** A stationary closest point or constrained boundary point was found. */
    converged_projected,
    /** The configured iteration limit was reached. */
    maximum_iterations,
    /** The pointwise geometric Jacobian lost reference-space rank. */
    rank_deficient,
    /** An input or an evaluated mapping quantity was not finite. */
    non_finite,
    /** The globalization step could not reduce the distance functional. */
    line_search_failed,
    /** The iteration stopped without satisfying an inverse contract. */
    stagnated
};

/** @enum InverseLocation
 * @brief Location of the recovered coordinate relative to the reference cell.
 */
enum class InverseLocation
{
    /** Strictly inside the reference cell. */
    inside,
    /** On the reference-cell boundary within the requested tolerance. */
    boundary,
    /** Outside the admissible reference cell. */
    outside,
    /** Not classified by the numerical kernel. */
    indeterminate
};

/** @enum InverseDomainPolicy
 * @brief Treatment of iterates outside the admissible reference cell.
 */
enum class InverseDomainPolicy
{
    /** Permit polynomial continuation outside the reference cell. */
    unrestricted,
    /** Project every trial point onto the closed reference cell. */
    constrained
};

/** @brief Return a stable textual name for an inverse status. */
[[nodiscard]] constexpr std::string_view
inverseStatusName( InverseStatus status ) noexcept
{
    switch ( status )
    {
    case InverseStatus::converged_exact:
        return "converged_exact";
    case InverseStatus::converged_projected:
        return "converged_projected";
    case InverseStatus::maximum_iterations:
        return "maximum_iterations";
    case InverseStatus::rank_deficient:
        return "rank_deficient";
    case InverseStatus::non_finite:
        return "non_finite";
    case InverseStatus::line_search_failed:
        return "line_search_failed";
    case InverseStatus::stagnated:
        return "stagnated";
    }
    return "unknown";
}

/** @brief Return a stable textual name for a recovered-point location. */
[[nodiscard]] constexpr std::string_view
inverseLocationName( InverseLocation location ) noexcept
{
    switch ( location )
    {
    case InverseLocation::inside:
        return "inside";
    case InverseLocation::boundary:
        return "boundary";
    case InverseLocation::outside:
        return "outside";
    case InverseLocation::indeterminate:
        return "indeterminate";
    }
    return "unknown";
}

/**
 * @struct InverseOptions
 * @brief Numerical and admissible-domain policy for geometric inversion.
 * @tparam Scalar floating-point scalar type
 */
template <typename Scalar>
struct InverseOptions
{
    /** Maximum Newton or Gauss--Newton iterations. */
    int maximum_iterations = 25;
    /** Maximum backtracking reductions per nonlinear iteration. */
    int maximum_line_search_iterations = 12;
    /** Absolute tolerance for the physical residual norm. */
    Scalar absolute_residual_tolerance = Scalar( 1e-12 );
    /** Target-scaled tolerance for the physical residual norm. */
    Scalar relative_residual_tolerance = Scalar( 1e-12 );
    /** Tolerance for the least-squares gradient norm. */
    Scalar gradient_tolerance = Scalar( 1e-12 );
    /** Relative tolerance for a vanishing reference-space update. */
    Scalar step_tolerance = Scalar( 1e-13 );
    /** Tolerance used to classify reference-cell boundaries. */
    Scalar reference_tolerance = Scalar( 1e-10 );
    /** Rank threshold passed to the column-pivoted QR factorization. */
    Scalar rank_tolerance = Scalar( 1e-12 );
    /** Armijo coefficient used by the backtracking globalization. */
    Scalar armijo_coefficient = Scalar( 1e-4 );
    /** Treatment of reference points outside the reference cell. */
    InverseDomainPolicy domain_policy = InverseDomainPolicy::unrestricted;
    /** Reuse the previous converged coordinate in batched operations. */
    bool warm_start = true;
};

/**
 * @struct InverseMapEvaluation
 * @brief Physical value and Jacobian returned by an inverse-map evaluator.
 * @tparam Scalar scalar type
 * @tparam ReferenceDim reference-space dimension
 * @tparam PhysicalDim physical-space dimension
 */
template <typename Scalar, int ReferenceDim, int PhysicalDim>
struct InverseMapEvaluation
{
    /** Physical image of the evaluated reference point. */
    Eigen::Matrix<Scalar, PhysicalDim, 1> physical_point;
    /** Pointwise physical Jacobian, with reference derivatives by columns. */
    Eigen::Matrix<Scalar, PhysicalDim, ReferenceDim> jacobian;

    /** @brief Return whether both value and Jacobian are finite. */
    [[nodiscard]] bool allFinite() const noexcept
    {
        return physical_point.allFinite() && jacobian.allFinite();
    }
};

/**
 * @struct InverseResult
 * @brief Complete result of one geometric inverse or closest-point operation.
 * @tparam Scalar scalar type
 * @tparam ReferenceDim reference-space dimension
 * @tparam PhysicalDim physical-space dimension
 */
template <typename Scalar, int ReferenceDim, int PhysicalDim>
struct InverseResult
{
    /** Recovered reference coordinate. */
    Eigen::Matrix<Scalar, ReferenceDim, 1> reference_point =
        Eigen::Matrix<Scalar, ReferenceDim, 1>::Zero();
    /** Physical image of `reference_point`. */
    Eigen::Matrix<Scalar, PhysicalDim, 1> physical_point =
        Eigen::Matrix<Scalar, PhysicalDim, 1>::Zero();
    /** Difference `physical_point-target_point`. */
    Eigen::Matrix<Scalar, PhysicalDim, 1> residual =
        Eigen::Matrix<Scalar, PhysicalDim, 1>::Zero();
    /** Euclidean physical distance to the target. */
    Scalar residual_norm = std::numeric_limits<Scalar>::infinity();
    /**
     * Least-squares stationarity norm. This is `|K^T residual|` for an
     * unrestricted solve and the projected-gradient norm for a constrained
     * solve.
     */
    Scalar gradient_norm = std::numeric_limits<Scalar>::infinity();
    /** Number of accepted nonlinear iterations. */
    int iterations = 0;
    /** Numerical termination status. */
    InverseStatus status = InverseStatus::maximum_iterations;
    /** Reference-cell classification filled by the geometric-map facade. */
    InverseLocation location = InverseLocation::indeterminate;
    /** Numerical rank of the last pointwise Jacobian. */
    int jacobian_rank = 0;
    /** Whether reference-domain projection was active. */
    bool constrained = false;
    /** Whether the initial coordinate came from a previous batch result. */
    bool used_warm_start = false;

    /** @brief Return whether a valid inverse or closest point was obtained. */
    [[nodiscard]] bool converged() const noexcept
    {
        return status == InverseStatus::converged_exact ||
               status == InverseStatus::converged_projected;
    }

    /** @brief Return whether the target lies on the mapped cell to tolerance. */
    [[nodiscard]] bool exact() const noexcept
    {
        return status == InverseStatus::converged_exact;
    }

    /** @brief Return whether the coordinate is admissible for the cell. */
    [[nodiscard]] bool admissible() const noexcept
    {
        return location == InverseLocation::inside ||
               location == InverseLocation::boundary;
    }

    /** @brief Return the Euclidean distance from the target to the result. */
    [[nodiscard]] Scalar distance() const noexcept
    {
        return residual_norm;
    }

    /** @brief Return whether the closest point lies on the cell boundary. */
    [[nodiscard]] bool onBoundary() const noexcept
    {
        return location == InverseLocation::boundary;
    }
};

namespace detail
{
/** @brief Identity reference-domain projector used by unrestricted solves. */
struct IdentityReferenceProjection
{
    template <typename VectorType>
    [[nodiscard]] VectorType operator()( VectorType const& point ) const
    {
        return point;
    }
};
} // namespace detail

/**
 * @brief Solve one geometric inverse or closest-point problem.
 *
 * Square maps use a Newton step. Embedded maps use the same rank-revealing
 * least-squares kernel, which is the Gauss--Newton step for
 * `0.5*|F(xHat)-x|^2`. A backtracking line search globalizes both cases.
 * Reference-domain projection is supplied by the geometric representation,
 * keeping topology policy out of this numerical kernel.
 *
 * @tparam ReferenceDim reference-space dimension
 * @tparam PhysicalDim physical-space dimension
 * @tparam Scalar scalar type
 * @tparam Evaluator callable returning `InverseMapEvaluation`
 * @tparam ReferenceProjection callable projecting a reference coordinate
 * @param target physical target point
 * @param initial_guess initial reference coordinate
 * @param evaluator physical map and Jacobian evaluator
 * @param options numerical and domain policy
 * @param project_reference reference-cell projection callable
 * @return structured inverse result
 */
template <int ReferenceDim, int PhysicalDim, typename Scalar,
          typename Evaluator,
          typename ReferenceProjection = detail::IdentityReferenceProjection>
    requires std::invocable<Evaluator,
                            Eigen::Matrix<Scalar, ReferenceDim, 1> const&> &&
             std::invocable<ReferenceProjection,
                            Eigen::Matrix<Scalar, ReferenceDim, 1> const&>
[[nodiscard]] InverseResult<Scalar, ReferenceDim, PhysicalDim>
solveInverseMap(
    Eigen::Matrix<Scalar, PhysicalDim, 1> const& target,
    Eigen::Matrix<Scalar, ReferenceDim, 1> const& initial_guess,
    Evaluator&& evaluator,
    InverseOptions<Scalar> const& options = {},
    ReferenceProjection&& project_reference = {} )
{
    using reference_vector_type = Eigen::Matrix<Scalar, ReferenceDim, 1>;
    using physical_vector_type = Eigen::Matrix<Scalar, PhysicalDim, 1>;
    using jacobian_type =
        Eigen::Matrix<Scalar, PhysicalDim, ReferenceDim>;
    using result_type = InverseResult<Scalar, ReferenceDim, PhysicalDim>;

    result_type result;
    result.constrained =
        options.domain_policy == InverseDomainPolicy::constrained;
    if ( !target.allFinite() || !initial_guess.allFinite() )
    {
        result.status = InverseStatus::non_finite;
        return result;
    }

    reference_vector_type reference_point = initial_guess;
    if ( result.constrained )
        reference_point = std::invoke( project_reference, reference_point );

    const Scalar residualTolerance =
        options.absolute_residual_tolerance +
        options.relative_residual_tolerance *
            std::max( Scalar( 1 ), target.norm() );

    auto evaluation = std::invoke( evaluator, reference_point );
    auto updateResult = [&]( auto const& currentEvaluation )
    {
        result.reference_point = reference_point;
        result.physical_point = currentEvaluation.physical_point;
        result.residual = result.physical_point - target;
        result.residual_norm = result.residual.norm();
        result.gradient_norm =
            ( currentEvaluation.jacobian.transpose() * result.residual ).norm();
    };

    if ( !evaluation.allFinite() )
    {
        result.status = InverseStatus::non_finite;
        return result;
    }
    updateResult( evaluation );
    Eigen::ColPivHouseholderQR<jacobian_type> initialFactorization(
        evaluation.jacobian );
    initialFactorization.setThreshold( options.rank_tolerance );
    result.jacobian_rank = initialFactorization.rank();
    if ( result.jacobian_rank < ReferenceDim )
    {
        result.status = InverseStatus::rank_deficient;
        return result;
    }
    if ( result.residual_norm <= residualTolerance )
    {
        result.status = InverseStatus::converged_exact;
        return result;
    }

    for ( int iteration = 0; iteration < options.maximum_iterations;
          ++iteration )
    {
        jacobian_type const& jacobian = evaluation.jacobian;
        Eigen::ColPivHouseholderQR<jacobian_type> factorization( jacobian );
        factorization.setThreshold( options.rank_tolerance );
        result.jacobian_rank = factorization.rank();
        if ( result.jacobian_rank < ReferenceDim )
        {
            result.status = InverseStatus::rank_deficient;
            return result;
        }

        physical_vector_type const residual =
            evaluation.physical_point - target;
        reference_vector_type const gradient =
            jacobian.transpose() * residual;
        if ( result.constrained )
        {
            reference_vector_type const projectedGradientStep =
                std::invoke( project_reference,
                             reference_point - gradient ) -
                reference_point;
            result.gradient_norm = projectedGradientStep.norm();
            Scalar const referenceScale =
                std::max( Scalar( 1 ), reference_point.norm() );
            Scalar const stationarityBound = std::max(
                options.gradient_tolerance * referenceScale,
                std::sqrt( std::numeric_limits<Scalar>::epsilon() ) *
                    referenceScale );
            if ( projectedGradientStep.norm() <= stationarityBound )
            {
                result.status = InverseStatus::converged_projected;
                return result;
            }
        }
        if constexpr ( PhysicalDim > ReferenceDim )
        {
            if ( gradient.norm() <= options.gradient_tolerance )
            {
                result.status = InverseStatus::converged_projected;
                return result;
            }
        }

        reference_vector_type const newtonStep =
            factorization.solve( -residual );
        if ( !newtonStep.allFinite() )
        {
            result.status = InverseStatus::non_finite;
            return result;
        }

        Scalar const objective = Scalar( 0.5 ) * residual.squaredNorm();
        bool accepted = false;
        bool vanishingStep = false;
        reference_vector_type acceptedPoint = reference_point;
        decltype( evaluation ) acceptedEvaluation = evaluation;
        auto tryDirection = [&]( reference_vector_type const& direction )
        {
            Scalar scale = Scalar( 1 );
            for ( int lineSearch = 0;
                  lineSearch < options.maximum_line_search_iterations;
                  ++lineSearch, scale *= Scalar( 0.5 ) )
            {
                reference_vector_type candidate =
                    reference_point + scale * direction;
                if ( result.constrained )
                    candidate = std::invoke( project_reference, candidate );
                reference_vector_type const actualStep =
                    candidate - reference_point;
                Scalar const stepBound = options.step_tolerance *
                                         std::max( Scalar( 1 ), reference_point.norm() );
                if ( actualStep.norm() <= stepBound )
                {
                    vanishingStep = true;
                    break;
                }

                auto candidateEvaluation = std::invoke( evaluator, candidate );
                if ( !candidateEvaluation.allFinite() )
                    continue;
                physical_vector_type const candidateResidual =
                    candidateEvaluation.physical_point - target;
                Scalar const candidateObjective =
                    Scalar( 0.5 ) * candidateResidual.squaredNorm();
                Scalar const armijoBound = objective +
                                           options.armijo_coefficient * gradient.dot( actualStep );
                if ( candidateObjective <= armijoBound ||
                     candidateObjective < objective )
                {
                    accepted = true;
                    acceptedPoint = candidate;
                    acceptedEvaluation = std::move( candidateEvaluation );
                    return;
                }
            }
        };

        tryDirection( newtonStep );
        if ( !accepted && result.constrained )
            tryDirection( -gradient );

        if ( !accepted )
        {
            result.status = vanishingStep
                                ? InverseStatus::stagnated
                                : InverseStatus::line_search_failed;
            return result;
        }

        reference_point = acceptedPoint;
        evaluation = std::move( acceptedEvaluation );
        result.iterations = iteration + 1;
        updateResult( evaluation );
        if ( result.residual_norm <= residualTolerance )
        {
            result.status = InverseStatus::converged_exact;
            return result;
        }
        if ( result.gradient_norm <= options.gradient_tolerance )
        {
            result.status = InverseStatus::converged_projected;
            return result;
        }
    }

    result.status = InverseStatus::maximum_iterations;
    return result;
}

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_GEOMETRICINVERSE_HPP
