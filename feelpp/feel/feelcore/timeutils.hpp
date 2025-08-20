/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2025-08-19

  Copyright (C) 2025 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef FEELPP_FEELCORE_TIME_UTILS_HPP
#define FEELPP_FEELCORE_TIME_UTILS_HPP 1

#include <chrono>
#include <ctime>
#include <fmt/chrono.h>

namespace Feel
{
/**
 * @brief Return current time as std::tm in UTC for use with {fmt}.
 *
 * Example:
 * @code{.cpp}
 *   #include <feel/feelcore/time_utils.hpp>
 *   #include <fmt/core.h>
 *   #include <fmt/chrono.h>
 *
 *   int main()
 *   {
 *       fmt::print("[{:%Y-%m-%d %H:%M:%S}] hello from Feel++!\\n", Feel::gmtimeNow());
 *       return 0;
 *   }
 * @endcode
 */
inline std::tm gmtimeNow()
{
    using clock_t = std::chrono::system_clock;
    auto now = clock_t::now();
    auto t   = clock_t::to_time_t(now);
    return fmt::gmtime(t);
}

/**
 * @brief Return a std::tm in UTC from a given system_clock::time_point.
 *
 * Example:
 * @code{.cpp}
 *   auto tp = std::chrono::system_clock::now();
 *   fmt::print("[{:%Y-%m-%d %H:%M:%S}] checkpoint\\n", Feel::gmtimeFrom(tp));
 * @endcode
 */
inline std::tm gmtimeFrom(std::chrono::system_clock::time_point tp)
{
    auto t = std::chrono::system_clock::to_time_t(tp);
    return fmt::gmtime(t);
}

} // namespace Feel

#endif // FEELPP_FEELCORE_TIME_UTILS_HPP