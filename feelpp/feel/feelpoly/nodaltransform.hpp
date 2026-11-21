/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file nodaltransform.hpp
 * @brief Checked scalar modal-to-nodal transforms and bounded immutable reference sharing.
 */
#ifndef FEELPP_FEELPOLY_NODALTRANSFORM_HPP
#define FEELPP_FEELPOLY_NODALTRANSFORM_HPP 1

#include <Eigen/QR>
#include <array>
#include <cmath>
#include <feel/feelpoly/eigenextent.hpp>
#include <feel/feelpoly/polynomialsolve.hpp>
#include <memory>
#include <mutex>

namespace Feel::detail
{
/** @brief Solver strategy scoped to nodal interpolation, independent of moment elements. */
enum class NodalSolveMethod
{
    Automatic,    ///< Try LU; retry with QR on conditioning or cardinal-residual concerns.
    PartialPivLU, ///< Checked partial-pivot LU, also available for matched comparisons.
    ColumnPivQR   ///< Rank-revealing column-pivoted QR.
};

/** @brief Numerical evidence retained with a solved nodal representation. */
template <typename T>
struct NodalSolveDiagnostics
{
    /** Selected factorization; Automatic is never a completed result. */
    NodalSolveMethod method = NodalSolveMethod::Automatic;
    /** Numerical rank; accepted LU has a checked nonsingular cardinal representation. */
    Eigen::Index rank = 0;
    /** Reciprocal one-norm condition number computed from the required solved transform. */
    T reciprocalCondition = 0;
    /** Relative Frobenius backward error of C V^T = I. */
    T backwardError = 0;
    /** Largest absolute entry of C V^T - I. */
    T kroneckerDefect = 0;
    /** Conditioning diagnostic; alone this does not reject a full-rank solve. */
    bool illConditioned = false;
};

/** @brief Own C = V^(-T), so every derivative table uses the same ell = C psi product.
 * @tparam T Real floating-point scalar type.
 * @note Solving against identity constructs the required basis representation once.
 * Evaluation never refactors V or forms another inverse. LU remains the fast path;
 * QR is retried below sqrt(epsilon) reciprocal conditioning or on a failed check.
 */
template <std::floating_point T>
class NodalPolynomialTransform
{
  public:
    /** Owned row-major coefficients, also suitable for a static-extent borrowed facade. */
    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    /** Numerical result accompanying the immutable coefficients. */
    using Diagnostics = NodalSolveDiagnostics<T>;

    /** @brief Factor a finite nonempty square V(i,j) = psi_j(x_i).
     * @throws std::invalid_argument For incompatible shape or non-finite data.
     * @throws std::runtime_error For numerical rank loss or an unacceptable residual.
     */
    template <typename Derived>
    explicit NodalPolynomialTransform( Eigen::MatrixBase<Derived> const& vandermonde,
                                       NodalSolveMethod method = NodalSolveMethod::Automatic )
        : M_vandermonde( vandermonde )
    {
        // Normalize before factorization and diagnostics to avoid scalar-range artifacts.
        if ( vandermonde.rows() == 0 || vandermonde.rows() != vandermonde.cols() || !vandermonde.allFinite() )
            throw std::invalid_argument( "nodal Vandermonde must be finite, nonempty and square" );
        T const scale = vandermonde.cwiseAbs().maxCoeff();
        if ( scale == 0 )
            throw std::runtime_error( "nodal Vandermonde is rank deficient" );
        Matrix const scaled = vandermonde / scale;

        // The established checked LU solve supplies the inexpensive candidate.
        bool accepted = false;
        if ( method != NodalSolveMethod::ColumnPivQR )
        {
            try
            {
                PolynomialLuSolver<T> solver( scaled, "nodal Vandermonde" );
                Matrix coefficients = solver.solve( Matrix::Identity( scaled.rows(), scaled.cols() ) ).transpose();
                auto diagnostics = diagnose( scaled, coefficients, NodalSolveMethod::PartialPivLU, scaled.rows() );
                accepted = acceptable( diagnostics, scaled.rows() ) &&
                           ( method == NodalSolveMethod::PartialPivLU || !diagnostics.illConditioned );
                if ( accepted )
                {
                    M_coefficients = std::move( coefficients );
                    M_diagnostics = diagnostics;
                }
            }
            catch ( std::runtime_error const& )
            {
                if ( method == NodalSolveMethod::PartialPivLU ) throw;
            }
            if ( !accepted && method == NodalSolveMethod::PartialPivLU )
                throw std::runtime_error( "nodal LU solve failed its cardinal residual check" );
        }

        // Rank-revealing QR retries difficult systems without a pseudoinverse or regularization.
        if ( !accepted )
        {
            Eigen::ColPivHouseholderQR<Matrix> solver( scaled );
            solver.setThreshold( Eigen::NumTraits<T>::epsilon() * scaled.rows() );
            if ( solver.rank() != scaled.rows() )
                throw std::runtime_error( "nodal Vandermonde is numerically rank deficient" );
            M_coefficients = solver.solve( Matrix::Identity( scaled.rows(), scaled.cols() ) ).transpose();
            M_diagnostics = diagnose( scaled, M_coefficients, NodalSolveMethod::ColumnPivQR, solver.rank() );
            if ( !acceptable( M_diagnostics, scaled.rows() ) )
                throw std::runtime_error( "nodal QR solve failed its cardinal residual check" );
        }

        // Retain the physical modal normalization in the sole published coefficient owner.
        M_coefficients /= scale;
        if ( !M_coefficients.allFinite() )
            throw std::runtime_error( "nodal coefficients exceed the scalar range" );
    }

    /** @brief Borrow solved coefficients until this owner is destroyed. */
    [[nodiscard]] Matrix const& coefficients() const& noexcept { return M_coefficients; }
    /** @brief Reject borrowing from a temporary transform. */
    Matrix const& coefficients() const&& = delete;
    /** @brief Borrow the original modal evaluation matrix for dual application and diagnostics. */
    [[nodiscard]] Matrix const& vandermonde() const& noexcept { return M_vandermonde; }
    /** @brief Reject borrowing from a temporary transform. */
    Matrix const& vandermonde() const&& = delete;
    /** @brief Return numerical diagnostics by value, independent of this owner's lifetime. */
    [[nodiscard]] Diagnostics diagnostics() const noexcept { return M_diagnostics; }

  private:
    /** @brief Compute forward cardinal and backward solve errors from independent products. */
    static Diagnostics diagnose( Matrix const& system, Matrix const& coefficients,
                                 NodalSolveMethod method, Eigen::Index rank )
    {
        Matrix const residual = coefficients * system.transpose() - Matrix::Identity( system.rows(), system.cols() );
        T const condition = system.cwiseAbs().colwise().sum().maxCoeff() *
                            coefficients.cwiseAbs().rowwise().sum().maxCoeff();
        T const reciprocal = T( 1 ) / condition;
        return { method, rank, reciprocal,
                 residual.stableNorm() / ( system.stableNorm() * coefficients.stableNorm() + std::sqrt( T( system.rows() ) ) ),
                 residual.cwiseAbs().maxCoeff(), reciprocal < std::sqrt( Eigen::NumTraits<T>::epsilon() ) };
    }

    /** @brief Require finite diagnostics and a backward-stable solve with a usable cardinal defect. */
    static bool acceptable( Diagnostics const& diagnostics, Eigen::Index size )
    {
        T const tolerance = T( 128 ) * Eigen::NumTraits<T>::epsilon() * size;
        return std::isfinite( diagnostics.reciprocalCondition ) && diagnostics.reciprocalCondition > 0 &&
               std::isfinite( diagnostics.backwardError ) && diagnostics.backwardError <= tolerance &&
               std::isfinite( diagnostics.kroneckerDefect ) &&
               diagnostics.kroneckerDefect <= std::max( tolerance, std::min( std::sqrt( Eigen::NumTraits<T>::epsilon() ), T( 32 ) * Eigen::NumTraits<T>::epsilon() / diagnostics.reciprocalCondition ) );
    }

    /** Original finite scalar Vandermonde in interpolation-node order. */
    Matrix M_vandermonde;
    /** Solved nodal coefficients in the orthonormal modal expansion. */
    Matrix M_coefficients;
    /** Evidence from the factorization that produced the published coefficients. */
    Diagnostics M_diagnostics;
};

/** @brief Reference data shared only by a canonical runtime scalar modal family.
 * @tparam CanonicalPrimal Runtime scalar primal type encoding dimension, cell and normalization.
 * @note Actual coordinates and effective degree complete the key. Field rank, continuity,
 * static/runtime facade and node-family label do not alter identical mathematical data.
 */
template <typename CanonicalPrimal>
struct LagrangeReferenceData
{
    /** Scalar type inherited from the mathematical modal family. */
    using T = typename CanonicalPrimal::value_type;
    /** Matrix owner used for exact coordinate identity checks. */
    using Matrix = typename NodalPolynomialTransform<T>::Matrix;
    /** Effective polynomial degree, independent of facade seed values. */
    uint16_type degree;
    /** Exact interpolation coordinates in local node order. */
    Matrix points;
    /** Immutable-after-publication scalar transform and numerical evidence. */
    NodalPolynomialTransform<T> transform;

    /** @brief Build only after a cache miss, consuming modal data before the caller expires. */
    template <typename Points, typename Build>
    LagrangeReferenceData( RuntimeOrder order, Points const& coordinates, Build&& build )
        : degree( order.value ), points( coordinates ), transform( build().transpose() ) {}
};

/** @brief Share a transform by exact configuration, with bounded weak retention and serialized publication.
 * @note At most sixteen weak handles are indexed per canonical family. Evicted live owners
 * remain valid; expired data are reconstructed. Exceptions publish nothing and release the lock.
 */
template <typename CanonicalPrimal>
class LagrangeReferenceCache
{
  public:
    /** Mathematical reference configuration retained by each published owner. */
    using Data = LagrangeReferenceData<CanonicalPrimal>;

    /** @brief Find or construct one exact effective-degree/coordinate configuration. */
    template <typename Points, typename Build>
    [[nodiscard]] static std::shared_ptr<Data const> get( RuntimeOrder order, Points const& points, Build&& build )
    {
        if ( !points.allFinite() )
            throw std::invalid_argument( "Lagrange interpolation coordinates must be finite" );

        // One lock prevents duplicate first factorization and protects weak-handle replacement.
        std::lock_guard lock( M_mutex );
        std::size_t available = M_next;
        for ( std::size_t i = 0; i < M_entries.size(); ++i )
        {
            if ( auto entry = M_entries[i].lock() )
            {
                if ( entry->degree == order.value && entry->points.rows() == points.rows() &&
                     entry->points.cols() == points.cols() && ( entry->points.array() == points.array() ).all() )
                    return entry;
            }
            else
                available = i;
        }

        // Publish only a fully checked transform; failure leaves all existing owners unchanged.
        auto entry = std::make_shared<Data const>( order, points, std::forward<Build>( build ) );
        M_entries[available] = entry;
        M_next = ( available + 1 ) % M_entries.size();
        return entry;
    }

  private:
    /** Fixed-size weak index; coefficient memory is owned exclusively by live users. */
    inline static std::array<std::weak_ptr<Data const>, 16> M_entries;
    /** Next replacement slot when every indexed configuration is live. */
    inline static std::size_t M_next = 0;
    /** Protect exact lookup, first construction and weak-index updates. */
    inline static std::mutex M_mutex;
};
} // namespace Feel::detail
#endif // FEELPP_FEELPOLY_NODALTRANSFORM_HPP
