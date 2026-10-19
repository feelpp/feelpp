/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

 SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
 SPDX-License-Identifier: LGPL-2.1-or-later
*/

#pragma once

#include <fmt/core.h>
#include <fmt/std.h>     // For std::filesystem::path, std::optional, etc.
#include <fmt/ranges.h>  // For STL containers (vector, map, set, etc.)
#include <fmt/ostream.h> // For types with ostream operators (fallback)

// Forward declare boost::ublas and GiNaC types
#if __has_include(<boost/numeric/ublas/vector.hpp>)
namespace boost { namespace numeric { namespace ublas {
    // Forward declarations matching actual ublas declarations
    template<class T, class A> class unbounded_array;
    template<class L, class T> class layout_base;
     template <class Z, class D>
    struct basic_row_major;
    template <class Z, class D>
    struct basic_column_major;
    
    template<class T, class A> class vector;
    template<class T, class L, class A> class matrix;
}}}

namespace fmt {
// Exclude all boost::ublas vector types from range detection
template <typename T, typename A>
struct is_range<boost::numeric::ublas::vector<T, A>, char> : std::false_type {};

// Exclude all boost::ublas matrix types from range detection  
template <typename T, typename L, typename A>
struct is_range<boost::numeric::ublas::matrix<T, L, A>, char> : std::false_type {};

// Prevent fmt/ranges.h from treating ublas types as container adaptors
#if FMT_VERSION < 120200
namespace detail {
template <typename T, typename A>
struct is_container_adaptor_like<boost::numeric::ublas::vector<T, A>> : std::false_type {};

template <typename T, typename L, typename A>
struct is_container_adaptor_like<boost::numeric::ublas::matrix<T, L, A>> : std::false_type {};
}
#endif
}

#if FMT_VERSION >= 120200
FMT_BEGIN_NAMESPACE
/** @brief Keep Boost.uBLAS vectors out of fmt container-adaptor formatting. */
template <typename T, typename A>
struct is_container_adaptor<boost::numeric::ublas::vector<T, A>> : std::false_type {};

/** @brief Keep Boost.uBLAS matrices out of fmt container-adaptor formatting. */
template <typename T, typename L, typename A>
struct is_container_adaptor<boost::numeric::ublas::matrix<T, L, A>> : std::false_type {};
FMT_END_NAMESPACE
#endif
#endif

#if __has_include(<ginac/ginac.h>)
namespace GiNaC { class ex; class symbol; class numeric; }

namespace fmt {
template <>
struct is_range<GiNaC::ex, char> : std::false_type {};
}
#endif


#include <filesystem>
#include <string>

// Custom formatters for feelcore types
namespace fmt
{

#if defined(FEELPP_HAS_HANA)
#include <boost/hana/integral_constant.hpp>

// Formatter for boost::hana::integral_constant
template <int N>
struct formatter<boost::hana::integral_constant<int, N>> : formatter<int> 
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx) { return formatter<int>::parse(ctx); }

    template <typename FormatContext>
    auto format(const boost::hana::integral_constant<int, N>& value, FormatContext& ctx) const
    {
        return formatter<int>::format(value(), ctx);
    }
};
#endif

#if defined(FEELPP_HAS_GOOGLE_GLOG) || defined(FEELPP_HAS_SPDLOG)
// Formatter for google::Counter_t (from glog or our compatibility shim)
namespace google { class Counter_t; }
template <> 
struct formatter<google::Counter_t> : ostream_formatter {};
#endif

// Provide ostream formatters for boost::ublas types (since we marked them as NOT ranges above)
#if __has_include(<boost/numeric/ublas/vector.hpp>)
template <typename T, typename A>
struct formatter<boost::numeric::ublas::vector<T, A>, char, void> : ostream_formatter {};

template <typename T, typename L, typename A>
struct formatter<boost::numeric::ublas::matrix<T, L, A>, char, void> : ostream_formatter {};
#endif

// Provide ostream formatters for GiNaC types (since we marked them as NOT ranges above)
#if __has_include(<ginac/ginac.h>)
template <>
struct formatter<GiNaC::ex> : ostream_formatter {};

template <>
struct formatter<GiNaC::symbol> : ostream_formatter {};

template <>
struct formatter<GiNaC::numeric> : ostream_formatter {};
#endif

} // namespace fmt
