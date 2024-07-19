/**
 * @file enumerate.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-07-21
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include <tuple>
#include <iterator>

namespace Feel {


/**
 * @brief create an index associated to a iterable object
 * 
 * @tparam T type of the iterable
 * @tparam TIter type of the iterator
 * @param start index start
 * @param iterable iterable object to create the index from
 * @return constexpr auto pair index and iterator
 */
template <typename T,
          typename TIter = decltype(std::begin(std::declval<T>())),
          typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T && iterable, size_t start=0)
{
    struct iterator
    {
        size_t i;
        TIter iter;
        bool operator != (const iterator & other) const { return iter != other.iter; }
        void operator ++ () { ++i; ++iter; }
        auto operator * () const { return std::tie(i, *iter); }
    };
    struct iterable_wrapper
    {
        size_t start;
        T iterable;
        auto begin() { return iterator{ start, std::begin(iterable) }; }
        auto end() { return iterator{ start, std::end(iterable) }; }
    };
    return iterable_wrapper{ start, std::forward<T>(iterable) };
}

} // Feel