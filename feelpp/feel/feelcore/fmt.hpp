/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 27 Oct 2024

 Copyright (C) 2024 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#pragma once

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <filesystem>
#include <string>

namespace fmt
{
template <>
struct fmt::formatter<std::filesystem::path> : fmt::formatter<std::string> 
{
    template <typename FormatContext>
    auto format(const std::filesystem::path& path, FormatContext& ctx) const
    {
        return fmt::formatter<std::string>::format(path.string(), ctx);
    }
};

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

// add google::Counter to the formatter
template <> struct formatter<google::Counter_t> : ostream_formatter {};



} // namespace fmt