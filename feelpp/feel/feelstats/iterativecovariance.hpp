/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileCopyrightText: 2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/
#ifndef FEELPP_STATS_ITERATIVECOVARIANCE_HPP
#define FEELPP_STATS_ITERATIVECOVARIANCE_HPP 1

#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace Feel::Stats
{

/**
 * @brief Online covariance and correlation of paired scalar samples.
 * @tparam ValueType Floating-point type used for both coordinates.
 */
template <std::floating_point ValueType = double>
class IterativeCovariance
{
public:
    /** @brief Type of one coordinate. */
    using value_type = ValueType;

    /** @brief Add one finite pair using a one-pass centered update. */
    void add( ValueType x, ValueType y )
    {
        if ( !std::isfinite( x ) || !std::isfinite( y ) )
            throw std::invalid_argument( "covariance requires finite values" );
        if ( M_Count == std::numeric_limits<std::uint64_t>::max() )
            throw std::overflow_error( "covariance sample count overflow" );
        ++M_Count;
        ValueType count = static_cast<ValueType>( M_Count );
        ValueType deltaX = x - M_MeanX;
        ValueType deltaY = y - M_MeanY;
        M_MeanX += deltaX / count;
        M_MeanY += deltaY / count;
        M_CrossDeviation += deltaX * ( y - M_MeanY );
        M_M2X += deltaX * ( x - M_MeanX );
        M_M2Y += deltaY * ( y - M_MeanY );
    }

    /** @brief Combine two independent paired streams with Chan's formulas. */
    void merge( IterativeCovariance const& other )
    {
        if ( other.empty() )
            return;
        if ( M_Count > std::numeric_limits<std::uint64_t>::max() - other.M_Count )
            throw std::overflow_error( "covariance sample count overflow" );
        if ( empty() )
        {
            *this = other;
            return;
        }
        ValueType leftCount = static_cast<ValueType>( M_Count );
        ValueType rightCount = static_cast<ValueType>( other.M_Count );
        ValueType combinedCount = leftCount + rightCount;
        ValueType factor = leftCount * rightCount / combinedCount;
        ValueType deltaX = other.M_MeanX - M_MeanX;
        ValueType deltaY = other.M_MeanY - M_MeanY;
        M_CrossDeviation += other.M_CrossDeviation + deltaX * deltaY * factor;
        M_M2X += other.M_M2X + deltaX * deltaX * factor;
        M_M2Y += other.M_M2Y + deltaY * deltaY * factor;
        M_MeanX += deltaX * rightCount / combinedCount;
        M_MeanY += deltaY * rightCount / combinedCount;
        M_Count += other.M_Count;
    }

    /** @brief Whether no pairs have been added. */
    bool empty() const { return M_Count == 0; }
    /** @brief Number of pairs. */
    std::uint64_t count() const { return M_Count; }
    /** @brief Mean of the first coordinate. */
    ValueType meanX() const { return M_MeanX; }
    /** @brief Mean of the second coordinate. */
    ValueType meanY() const { return M_MeanY; }
    /** @brief Population covariance, or zero when empty. */
    ValueType covariance() const { return empty() ? 0 : M_CrossDeviation / static_cast<ValueType>( M_Count ); }
    /** @brief Sample covariance, or zero with fewer than two pairs. */
    ValueType sampleCovariance() const { return M_Count < 2 ? 0 : M_CrossDeviation / static_cast<ValueType>( M_Count - 1 ); }
    /** @brief Pearson correlation, or zero when either coordinate is constant. */
    ValueType correlation() const
    {
        return M_M2X <= 0 || M_M2Y <= 0 ? 0 : M_CrossDeviation / std::sqrt( M_M2X * M_M2Y );
    }

    /** @brief Serialize count, means, and centered sums. */
    template <typename Archive>
    void serialize( Archive& archive, unsigned int )
    {
        archive & M_Count & M_MeanX & M_MeanY & M_CrossDeviation & M_M2X & M_M2Y;
    }

private:
    /// Number of paired observations.
    std::uint64_t M_Count = 0;
    /// Mean of the first coordinate.
    ValueType M_MeanX = 0;
    /// Mean of the second coordinate.
    ValueType M_MeanY = 0;
    /// Sum of cross deviations from the means.
    ValueType M_CrossDeviation = 0;
    /// Sum of squared deviations of the first coordinate.
    ValueType M_M2X = 0;
    /// Sum of squared deviations of the second coordinate.
    ValueType M_M2Y = 0;
};

} // namespace Feel::Stats

#endif
