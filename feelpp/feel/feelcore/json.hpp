//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 05 Oct 2020
//! @copyright 2020 Feel++ Consortium
//!
#pragma once
#include <optional>
#include <feel/feelcore/feel.hpp>
#define JSON_DIAGNOSTICS FEELPP_ENABLE_JSON_DIAGNOSTICS
#include <feel/feelcore/_json.hpp>
namespace Feel
{
    namespace nl = nlohmann;
    using json = nl::json;
}

// partial specialization (full specialization works too)
namespace nlohmann
{
template<>
struct adl_serializer<Feel::fs::path>
{
    static void to_json( json& j, const Feel::fs::path& p )
    {
        j=p.string();
    }

    static void from_json( const json& j, Feel::fs::path& p )
    {
        p = j.get<std::string>();
    }
};


//
// std::optional
//

template <class T>
void to_json( nlohmann::json& j, const std::optional<T>& v )
{
    if ( v.has_value() )
        j = *v;
    else
        j = nullptr;
}

template <class T>
void from_json( const nlohmann::json& j, std::optional<T>& v )
{
    if ( j.is_null() )
        v = std::nullopt;
    else
        v = j.get<T>();
}

} // namespace nlohmann
