/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileCopyrightText: 2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/
#ifndef FEELPP_STATS_ITERATIVEEXTREMA_HPP
#define FEELPP_STATS_ITERATIVEEXTREMA_HPP 1

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace Feel::Stats
{

/**
 * @brief Online minimum and maximum of an ordered sample stream.
 * @tparam ValueType Copyable, totally ordered sample type.
 */
template <typename ValueType>
    requires std::totally_ordered<ValueType> && std::copyable<ValueType> && std::default_initializable<ValueType>
class IterativeExtrema
{
public:
    /** @brief Type of one observation. */
    using value_type = ValueType;

    /** @brief Add one sample. */
    void add( ValueType const& value )
    {
        if constexpr ( std::floating_point<ValueType> )
        {
            if ( !std::isfinite( value ) )
                throw std::invalid_argument( "extrema require finite values" );
        }
        if ( M_Count == std::numeric_limits<std::uint64_t>::max() )
            throw std::overflow_error( "extrema sample count overflow" );
        if ( empty() )
            M_Min = M_Max = value;
        else
        {
            M_Min = std::min( M_Min, value );
            M_Max = std::max( M_Max, value );
        }
        ++M_Count;
    }

    /** @brief Merge another set of extrema. */
    void merge( IterativeExtrema const& other )
    {
        if ( other.empty() )
            return;
        if ( M_Count > std::numeric_limits<std::uint64_t>::max() - other.M_Count )
            throw std::overflow_error( "extrema sample count overflow" );
        if ( empty() )
        {
            *this = other;
            return;
        }
        M_Min = std::min( M_Min, other.M_Min );
        M_Max = std::max( M_Max, other.M_Max );
        M_Count += other.M_Count;
    }

    /** @brief Whether no samples have been added. */
    bool empty() const { return M_Count == 0; }
    /** @brief Number of samples. */
    std::uint64_t count() const { return M_Count; }
    /** @brief Exact minimum, valid for a nonempty set. */
    ValueType const& min() const { return M_Min; }
    /** @brief Exact maximum, valid for a nonempty set. */
    ValueType const& max() const { return M_Max; }

    /** @brief Serialize the extrema state. */
    template <typename Archive>
    void serialize( Archive& archive, unsigned int )
    {
        archive & M_Count & M_Min & M_Max;
    }

private:
    /// Number of input values.
    std::uint64_t M_Count = 0;
    /// Smallest input value.
    ValueType M_Min{};
    /// Largest input value.
    ValueType M_Max{};
};

} // namespace Feel::Stats

#endif
