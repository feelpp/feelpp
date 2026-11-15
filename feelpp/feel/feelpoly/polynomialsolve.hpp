/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file polynomialsolve.hpp
 * @brief Shared Eigen LU solves with Release shape, singularity and residual diagnostics.
 */
#ifndef FEELPP_FEELPOLY_POLYNOMIALSOLVE_HPP
#define FEELPP_FEELPOLY_POLYNOMIALSOLVE_HPP 1
#include <Eigen/Core>
#include <Eigen/LU>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <string_view>

namespace Feel::detail
{
/** @brief Reuse a partial-pivot LU factorization for polynomial interpolation systems.
 * @tparam T Scalar type.
 *
 * This preserves the existing LU algorithm family. It owns the scaled system
 * for residual checks and the Eigen factorization; no inverse is formed to
 * apply a solve. QR-versus-LU selection remains a separate numerical decision.
 */
template<typename T>
class PolynomialLuSolver
{
public:
    /** Dynamic owning matrix used by Eigen's factorization workspace. */
    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    /** Real scalar used for norms and residual bounds. */
    using Real = typename Eigen::NumTraits<T>::Real;

    /** @brief Validate and factor a finite square polynomial system. */
    template<typename Derived>
    explicit PolynomialLuSolver( Eigen::MatrixBase<Derived> const& system, std::string_view context )
        : M_system( system ), M_context( context )
    {
        if ( M_system.rows() != M_system.cols() )
            throw std::invalid_argument( M_context + ": interpolation matrix must be square" );
        if ( !M_system.allFinite() )
            throw std::invalid_argument( M_context + ": interpolation matrix must be finite" );
        if ( M_system.rows() == 0 ) return;
        M_scale = M_system.cwiseAbs().maxCoeff();
        if ( M_scale == Real( 0 ) )
            throw std::runtime_error( M_context + ": interpolation matrix is singular" );
        M_system /= M_scale;
        M_decomposition.compute( M_system );
        if ( !M_decomposition.matrixLU().allFinite() ||
             ( M_decomposition.matrixLU().diagonal().array().abs() == Real( 0 ) ).any() )
            throw std::runtime_error( M_context + ": interpolation matrix is singular" );
    }

    /** @brief Solve one or more right-hand sides and verify scaled backward error.
     * @throws std::invalid_argument For incompatible or non-finite right-hand sides.
     * @throws std::runtime_error If the computed solution fails its residual check.
     */
    template<typename Derived>
    [[nodiscard]] Matrix solve( Eigen::MatrixBase<Derived> const& rhs ) const
    { return solveSystem<false>( rhs ); }

    /** @brief Solve X A = B using the same factorization, without forming an inverse. */
    template<typename Derived>
    [[nodiscard]] Matrix solveRight( Eigen::MatrixBase<Derived> const& rhs ) const
    { return solveSystem<true>( rhs.transpose() ).transpose().eval(); }

private:
    /** @brief Share scaling and backward-error checks for ordinary and transposed solves. */
    template<bool Transposed, typename Derived>
    [[nodiscard]] Matrix solveSystem( Eigen::MatrixBase<Derived> const& rhs ) const
    {
        if ( rhs.rows() != M_system.rows() )
            throw std::invalid_argument( M_context + ": interpolation right-hand side has wrong row count" );
        if ( !rhs.allFinite() )
            throw std::invalid_argument( M_context + ": interpolation right-hand side must be finite" );
        if ( M_system.rows() == 0 ) return Matrix( 0, rhs.cols() );
        if ( rhs.cols() == 0 ) return Matrix( M_system.rows(), 0 );
        Matrix scaledRhs = rhs / M_scale;
        if ( !scaledRhs.allFinite() )
            throw std::runtime_error( M_context + ": right-hand-side scaling exceeds the scalar range" );
        Matrix result;
        if constexpr ( Transposed ) result = M_decomposition.transpose().solve( scaledRhs );
        else result = M_decomposition.solve( scaledRhs );
        if ( !result.allFinite() )
            throw std::runtime_error( M_context + ": interpolation solution is not finite" );
        Real const magnitude = std::max( result.cwiseAbs().maxCoeff(), scaledRhs.cwiseAbs().maxCoeff() );
        if ( magnitude == Real( 0 ) ) return result;
        // Evaluate backward error after normalization, avoiding overflowing norms.
        Matrix normalizedResult = result / magnitude;
        Matrix normalizedRhs = scaledRhs / magnitude;
        Matrix residual;
        if constexpr ( Transposed ) residual = M_system.transpose() * normalizedResult - normalizedRhs;
        else residual = M_system * normalizedResult - normalizedRhs;
        Real const scale = M_system.stableNorm() * normalizedResult.stableNorm() + normalizedRhs.stableNorm();
        Real const tolerance = Real( 128 ) * Eigen::NumTraits<Real>::epsilon() *
            Real( std::max<Eigen::Index>( 1, M_system.rows() ) ) * scale;
        if ( !residual.allFinite() ||
             !( residual.stableNorm() <= tolerance ) )
            throw std::runtime_error( M_context + ": interpolation solve failed its residual check" );
        return result;
    }

    /** Scaled original system retained for independent residual evaluation. */
    Matrix M_system;
    /** Reusable Eigen partial-pivot LU decomposition. */
    Eigen::PartialPivLU<Matrix> M_decomposition;
    /** Global scaling factor, preserving the LU pivoting algorithm. */
    Real M_scale = Real( 1 );
    /** Owned diagnostic context, independent of caller string lifetimes. */
    std::string M_context;
};
} // namespace Feel::detail
#endif // FEELPP_FEELPOLY_POLYNOMIALSOLVE_HPP
