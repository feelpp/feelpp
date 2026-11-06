/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileCopyrightText: 2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/
#ifndef FEELPP_STATS_ITERATIVEMOMENTS_HPP
#define FEELPP_STATS_ITERATIVEMOMENTS_HPP 1

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace Feel::Stats
{

/**
 * @brief Streaming central moments through fourth order, with parallel merge.
 *
 * Uses centered one-pass updates and the parallel formulas described in
 * Pébay, SAND2008-6212. Kurtosis is the population fourth standardized moment.
 * @tparam ValueType Floating-point type used for moments.
 */
template <std::floating_point ValueType = double>
class IterativeMoments
{
public:
    /** @brief Type of one observation. */
    using value_type = ValueType;

    /** @brief Add one finite sample using centered online updates. */
    void add( ValueType value )
    {
        if ( !std::isfinite( value ) )
            throw std::invalid_argument( "moments require a finite value" );
        if ( M_Count == std::numeric_limits<std::uint64_t>::max() )
            throw std::overflow_error( "moments sample count overflow" );
        ValueType previousCount = static_cast<ValueType>( M_Count );
        ++M_Count;
        ValueType count = static_cast<ValueType>( M_Count );
        ValueType delta = value - M_Mean;
        ValueType deltaOverCount = delta / count;
        ValueType deltaOverCount2 = deltaOverCount * deltaOverCount;
        ValueType term = delta * deltaOverCount * previousCount;
        M_M4 += term * deltaOverCount2 * ( count * count - 3 * count + 3 )
                + 6 * deltaOverCount2 * M_M2 - 4 * deltaOverCount * M_M3;
        M_M3 += term * deltaOverCount * ( count - 2 ) - 3 * deltaOverCount * M_M2;
        M_M2 += term;
        M_Mean += deltaOverCount;
    }

    /** @brief Merge another sample set using parallel central-moment formulas. */
    void merge( IterativeMoments const& other )
    {
        if ( other.empty() )
            return;
        if ( M_Count > std::numeric_limits<std::uint64_t>::max() - other.M_Count )
            throw std::overflow_error( "moments sample count overflow" );
        if ( empty() )
        {
            *this = other;
            return;
        }
        ValueType a = static_cast<ValueType>( M_Count );
        ValueType b = static_cast<ValueType>( other.M_Count );
        ValueType n = a + b;
        ValueType delta = other.M_Mean - M_Mean;
        ValueType delta2 = delta * delta;
        ValueType delta3 = delta2 * delta;
        ValueType delta4 = delta2 * delta2;
        ValueType mergedM4 = M_M4 + other.M_M4
            + delta4 * a * b * ( a * a - a * b + b * b ) / ( n * n * n )
            + 6 * delta2 * ( a * a * other.M_M2 + b * b * M_M2 ) / ( n * n )
            + 4 * delta * ( a * other.M_M3 - b * M_M3 ) / n;
        ValueType mergedM3 = M_M3 + other.M_M3
            + delta3 * a * b * ( a - b ) / ( n * n )
            + 3 * delta * ( a * other.M_M2 - b * M_M2 ) / n;
        M_M2 += other.M_M2 + delta2 * a * b / n;
        M_M3 = mergedM3;
        M_M4 = mergedM4;
        M_Mean += delta * b / n;
        M_Count += other.M_Count;
    }

    /** @brief Whether no samples have been added. */
    bool empty() const { return M_Count == 0; }
    /** @brief Number of samples. */
    std::uint64_t count() const { return M_Count; }
    /** @brief Arithmetic mean. */
    ValueType mean() const { return M_Mean; }
    /** @brief Sum of second central powers. */
    ValueType m2() const { return M_M2; }
    /** @brief Sum of third central powers. */
    ValueType m3() const { return M_M3; }
    /** @brief Sum of fourth central powers. */
    ValueType m4() const { return M_M4; }
    /** @brief Population variance, or zero when empty. */
    ValueType variance() const { return empty() ? 0 : std::max( M_M2, ValueType{0} ) / static_cast<ValueType>( M_Count ); }
    /** @brief Standardized population skewness, or zero for zero variance. */
    ValueType skewness() const
    {
        return M_M2 <= 0 ? 0 : std::sqrt( static_cast<ValueType>( M_Count ) ) * M_M3 / std::pow( M_M2, ValueType{1.5} );
    }
    /** @brief Population kurtosis, or zero for zero variance. */
    ValueType kurtosis() const
    {
        return M_M2 <= 0 ? 0 : static_cast<ValueType>( M_Count ) * M_M4 / ( M_M2 * M_M2 );
    }
    /** @brief Kurtosis minus the normal-distribution value of three. */
    ValueType excessKurtosis() const { return M_M2 <= 0 ? 0 : kurtosis() - 3; }

    /** @brief Serialize the count and four central-moment accumulators. */
    template <typename Archive>
    void serialize( Archive& archive, unsigned int )
    {
        archive & M_Count & M_Mean & M_M2 & M_M3 & M_M4;
    }

private:
    /// Number of samples.
    std::uint64_t M_Count = 0;
    /// Arithmetic mean.
    ValueType M_Mean = 0;
    /// Sum of second central powers.
    ValueType M_M2 = 0;
    /// Sum of third central powers.
    ValueType M_M3 = 0;
    /// Sum of fourth central powers.
    ValueType M_M4 = 0;
};

} // namespace Feel::Stats

#endif
