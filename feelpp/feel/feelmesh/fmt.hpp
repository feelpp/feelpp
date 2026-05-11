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
//! @file feel/feelmesh/fmt.hpp
//! @brief Formatter specialization for Mesh objects
//! @date 1 Nov 2024
//! @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
//! @copyright (C) 2024 Feel++ Consortium
//!
#pragma once

#include <fmt/format.h>

#include <feel/feelmesh/meshbase.hpp>
#include <feel/feelmesh/ranges.hpp>

namespace fmt
{
/**
 * \brief Formatter specialization for Feel::Range.
 *
 * This struct provides a custom formatter for Feel::Range objects, enabling
 * them to be formatted using the {fmt} library.
 *
 * \tparam MeshType The type of the mesh.
 * \tparam MESH_ENTITIES The number of mesh entities.
 *
 * This formatter is enabled only if MeshType is derived from Feel::MeshBase.
 */
template <typename MeshType, int MESH_ENTITIES>
struct formatter<Feel::Range<MeshType, MESH_ENTITIES>, std::enable_if_t<std::is_base_of_v<Feel::MeshBase<>, MeshType>, char>>
{
    /**
     * \brief Parses format specifications.
     *
     * This function is called to parse any format specifications that might
     * be provided. In this case, it simply returns the beginning of the
     * context, as no specific format specifications are needed.
     *
     * \param ctx The format parse context.
     * \return An iterator to the beginning of the context.
     */
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return ctx.begin();
    }

    /**
     * \brief Formats a Feel::Range object.
     *
     * This function is called to format a Feel::Range object. It accesses
     * properties from MeshType and Feel::Range to create a formatted string.
     *
     * \param range The Feel::Range object to format.
     * \param ctx The format context.
     * \return An iterator to the end of the formatted output.
     */
    template <typename FormatContext>
    auto format(const Feel::Range<MeshType, MESH_ENTITIES>& range, FormatContext& ctx) const -> decltype(ctx.out())
    {
        return fmt::format_to(ctx.out(), "Range<{}>(entities={}, size={})",
                                typeid(MeshType).name(), MESH_ENTITIES, range.size());
    }
};

}