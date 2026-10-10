/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileCopyrightText: 2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/
#ifndef FEELPP_STATS_ITERATIVE_STATISTICS_HPP
#define FEELPP_STATS_ITERATIVE_STATISTICS_HPP 1

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace Feel::Stats
{

/**
 * @brief Streaming moments and extrema, mergeable with Chan's formulas.
 * @tparam ValueType Floating-point type used for samples and statistics.
 */
template <std::floating_point ValueType = double>
class IterativeStatistics
{
public:
    /** @brief Type of one observation. */
    using value_type = ValueType;

    /** @brief Add one finite sample using Welford's update. */
    void add( ValueType value )
    {
        if ( !std::isfinite( value ) )
            throw std::invalid_argument( "statistics require a finite value" );
        if ( M_Count == std::numeric_limits<std::uint64_t>::max() )
            throw std::overflow_error( "statistics sample count overflow" );
        if ( empty() )
            M_Min = M_Max = value;
        else
        {
            M_Min = std::min( M_Min, value );
            M_Max = std::max( M_Max, value );
        }
        ++M_Count;
        addToTotal( value );
        ValueType delta = value - M_Mean;
        M_Mean += delta / static_cast<ValueType>( M_Count );
        M_M2 += delta * ( value - M_Mean );
    }

    /** @brief Combine independent sample sets using Chan's update. */
    void merge( IterativeStatistics const& other )
    {
        if ( other.empty() )
            return;
        if ( M_Count > std::numeric_limits<std::uint64_t>::max() - other.M_Count )
            throw std::overflow_error( "statistics sample count overflow" );
        if ( empty() )
        {
            *this = other;
            return;
        }
        ValueType leftCount = static_cast<ValueType>( M_Count );
        ValueType rightCount = static_cast<ValueType>( other.M_Count );
        ValueType combinedCount = leftCount + rightCount;
        ValueType delta = other.M_Mean - M_Mean;
        M_M2 += other.M_M2 + delta * delta * ( leftCount * rightCount / combinedCount );
        M_Mean += delta * ( rightCount / combinedCount );
        M_Min = std::min( M_Min, other.M_Min );
        M_Max = std::max( M_Max, other.M_Max );
        M_Count += other.M_Count;
        addToTotal( other.total() );
    }

    /** @brief Whether the sample set is empty. */
    bool empty() const { return M_Count == 0; }
    /** @brief Number of samples. */
    std::uint64_t count() const { return M_Count; }
    /** @brief Compensated sum of samples. */
    ValueType total() const { return M_Total + M_TotalCompensation; }
    /** @brief Exact minimum, valid for a nonempty set. */
    ValueType min() const { return M_Min; }
    /** @brief Exact maximum, valid for a nonempty set. */
    ValueType max() const { return M_Max; }
    /** @brief Arithmetic mean, or zero for an empty set. */
    ValueType mean() const { return M_Mean; }
    /** @brief Sum of squared deviations from the mean. */
    ValueType m2() const { return M_M2; }
    /** @brief Population variance, or zero for an empty set. */
    ValueType variance() const { return empty() ? 0 : std::max( M_M2, ValueType{0} ) / static_cast<ValueType>( M_Count ); }
    /** @brief Population standard deviation. */
    ValueType stdDev() const { return std::sqrt( variance() ); }

    /** @brief Serialize state for MPI or persistent storage. */
    template <typename Archive>
    void serialize( Archive& archive, unsigned int )
    {
        archive & M_Count & M_Total & M_TotalCompensation & M_Min & M_Max & M_Mean & M_M2;
    }

private:
    /** @brief Accumulate a value with Neumaier compensation. */
    void addToTotal( ValueType value )
    {
        ValueType updated = M_Total + value;
        M_TotalCompensation += std::abs( M_Total ) >= std::abs( value )
                                   ? ( M_Total - updated ) + value
                                   : ( value - updated ) + M_Total;
        M_Total = updated;
    }

    /// Number of samples.
    std::uint64_t M_Count = 0;
    /// Leading part of the sum.
    ValueType M_Total = 0;
    /// Low-order correction to the sum.
    ValueType M_TotalCompensation = 0;
    /// Smallest sample.
    ValueType M_Min = 0;
    /// Largest sample.
    ValueType M_Max = 0;
    /// Running arithmetic mean.
    ValueType M_Mean = 0;
    /// Sum of squared deviations.
    ValueType M_M2 = 0;
};

} // namespace Feel::Stats

#endif
