/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-FileCopyrightText: 2026 University of Strasbourg
   SPDX-License-Identifier: LGPL-2.1-or-later
*/
#ifndef FEELPP_STATS_KLL_HPP
#define FEELPP_STATS_KLL_HPP 1

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <concepts>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/serialization/vector.hpp>

namespace Feel::Stats
{

/**
 * @brief Bounded, mergeable quantile sketch using randomized KLL compaction.
 *
 * Each level retains at most @c capacity values after compaction; a value at
 * level i represents 2^i observations. Quantiles use weighted order statistics
 * and return a retained observation. Exact extrema are preserved separately.
 * This equal-capacity variant follows the KLL compaction idea described at
 * https://datasketches.apache.org/docs/KLL/KLLSketch.html. It uses more space
 * than Apache DataSketches' optimized implementation and does not claim its
 * published normalized-rank error bound.
 * @tparam ValueType Copyable, totally ordered sample type.
 */
template <typename ValueType>
    requires std::totally_ordered<ValueType> && std::copyable<ValueType> && std::default_initializable<ValueType>
class KllSketch
{
public:
    /** @brief Type of one observation. */
    using value_type = ValueType;

    /** @brief Create a sketch with capacity 200 per level. */
    KllSketch() : KllSketch( 200 ) {}

    /** @brief Create a sketch with the given capacity per level. */
    explicit KllSketch( std::size_t capacity )
        : M_Capacity( capacity ), M_Levels( 1 )
    {
        if ( capacity < 8 )
            throw std::invalid_argument( "KLL capacity must be at least 8" );
    }

    /** @brief Insert one observation. */
    void add( ValueType const& value )
    {
        if constexpr ( std::floating_point<ValueType> )
        {
            if ( !std::isfinite( value ) )
                throw std::invalid_argument( "KLL requires finite values" );
        }
        if ( M_Count == std::numeric_limits<std::uint64_t>::max() )
            throw std::overflow_error( "KLL sample count overflow" );
        if ( empty() )
            M_Min = M_Max = value;
        else
        {
            if ( std::less<ValueType>{}( value, M_Min ) )
                M_Min = value;
            if ( std::less<ValueType>{}( M_Max, value ) )
                M_Max = value;
        }
        ++M_Count;
        M_Levels[0].push_back( value );
        compact();
    }

    /** @brief Merge another sketch into this one. */
    void merge( KllSketch const& other )
    {
        if ( this == &other )
        {
            KllSketch copy( other );
            merge( copy );
            return;
        }
        if ( other.empty() )
            return;
        if ( M_Count > std::numeric_limits<std::uint64_t>::max() - other.M_Count )
            throw std::overflow_error( "KLL sample count overflow" );
        if ( empty() )
        {
            M_Min = other.M_Min;
            M_Max = other.M_Max;
        }
        else
        {
            if ( std::less<ValueType>{}( other.M_Min, M_Min ) )
                M_Min = other.M_Min;
            if ( std::less<ValueType>{}( M_Max, other.M_Max ) )
                M_Max = other.M_Max;
        }
        M_Count += other.M_Count;
        if ( M_Levels.size() < other.M_Levels.size() )
            M_Levels.resize( other.M_Levels.size() );
        for ( std::size_t level = 0; level < other.M_Levels.size(); ++level )
            M_Levels[level].insert( M_Levels[level].end(), other.M_Levels[level].begin(), other.M_Levels[level].end() );
        compact();
    }

    /** @brief Whether no observations have been inserted. */
    bool empty() const { return M_Count == 0; }
    /** @brief Number of represented observations. */
    std::uint64_t count() const { return M_Count; }
    /** @brief Maximum number of retained values per level. */
    std::size_t capacity() const { return M_Capacity; }
    /** @brief Exact minimum, valid for a nonempty sketch. */
    ValueType min() const { return M_Min; }
    /** @brief Exact maximum, valid for a nonempty sketch. */
    ValueType max() const { return M_Max; }

    /** @brief Number of currently retained values across all levels. */
    std::size_t retainedCount() const
    {
        std::size_t result = 0;
        for ( auto const& level : M_Levels )
            result += level.size();
        return result;
    }

    /**
     * @brief Estimate a quantile by Apache KLL's inclusive weighted rank.
     * @param probability Rank fraction from zero through one.
     * @return A retained observation, or an exact extremum at zero or one.
     * @throws std::invalid_argument If empty or probability is outside [0, 1].
     */
    ValueType quantile( double probability ) const
    {
        if ( empty() || !( probability >= 0 && probability <= 1 ) )
            throw std::invalid_argument( "KLL quantile requires samples and a probability in [0, 1]" );
        if ( probability == 0 )
            return M_Min;
        if ( probability == 1 )
            return M_Max;

        std::vector<std::pair<ValueType, std::uint64_t>> weighted;
        weighted.reserve( retainedCount() );
        for ( std::size_t level = 0; level < M_Levels.size(); ++level )
        {
            std::uint64_t weight = std::uint64_t{1} << level;
            for ( auto const& value : M_Levels[level] )
                weighted.emplace_back( value, weight );
        }
        std::sort( weighted.begin(), weighted.end(), []( auto const& left, auto const& right )
                   { return std::less<ValueType>{}( left.first, right.first ); } );
        std::uint64_t target = static_cast<std::uint64_t>( std::ceil( probability * static_cast<long double>( M_Count ) ) );
        std::uint64_t accumulated = 0;
        for ( auto const& [value, weight] : weighted )
        {
            accumulated += weight;
            if ( accumulated >= target )
                return value;
        }
        return weighted.back().first;
    }

    /** @brief Serialize the sketch, including levels and random state. */
    template <typename Archive>
    void serialize( Archive& archive, unsigned int )
    {
        archive & M_Capacity & M_Count & M_Min & M_Max & M_Levels & M_RandomState;
    }

private:
    /** @brief Draw one deterministic pseudorandom bit for compaction. */
    bool randomBit()
    {
        M_RandomState ^= M_RandomState << 13;
        M_RandomState ^= M_RandomState >> 7;
        M_RandomState ^= M_RandomState << 17;
        return ( M_RandomState & 1 ) != 0;
    }

    /** @brief Promote randomly selected values until all levels fit. */
    void compact()
    {
        for ( std::size_t level = 0; level < M_Levels.size(); ++level )
        {
            while ( M_Levels[level].size() > M_Capacity )
            {
                if ( level >= 63 )
                    throw std::overflow_error( "KLL level overflow" );
                if ( level + 1 == M_Levels.size() )
                    M_Levels.emplace_back();

                auto& values = M_Levels[level];
                std::sort( values.begin(), values.end(), std::less<ValueType>{} );
                std::size_t first = 0;
                std::size_t last = values.size();
                std::vector<ValueType> leftover;
                if ( values.size() % 2 != 0 )
                {
                    if ( randomBit() )
                        leftover.push_back( values[first++] );
                    else
                        leftover.push_back( values[--last] );
                }
                std::size_t offset = randomBit() ? 1 : 0;
                auto& upper = M_Levels[level + 1];
                for ( std::size_t i = first + offset; i < last; i += 2 )
                    upper.push_back( values[i] );
                values.swap( leftover );
            }
        }
    }

    /// Retained-value capacity of each level.
    std::size_t M_Capacity;
    /// Number of represented input observations.
    std::uint64_t M_Count = 0;
    /// Exact smallest input observation.
    ValueType M_Min{};
    /// Exact largest input observation.
    ValueType M_Max{};
    /// Values retained at each power-of-two weight.
    std::vector<std::vector<ValueType>> M_Levels;
    /// State of the compaction coin generator.
    std::uint64_t M_RandomState = 0x9e3779b97f4a7c15ULL;
};

} // namespace Feel::Stats

#endif
