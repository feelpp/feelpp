/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file submeshcontext.hpp
 * @brief Shared extraction options for mesh and submesh APIs.
 */
#ifndef FEELPP_DISCR_SUBMESHCONTEXT_HPP
#define FEELPP_DISCR_SUBMESHCONTEXT_HPP 1

#include <feel/feelcore/feeltypes.hpp>

namespace Feel
{
inline constexpr size_type EXTRACTION_KEEP_POINTS_IDS = ( 1 << 0 );
inline constexpr size_type EXTRACTION_KEEP_EDGES_IDS = ( 1 << 1 );
inline constexpr size_type EXTRACTION_KEEP_FACES_IDS = ( 1 << 2 );
inline constexpr size_type EXTRACTION_KEEP_VOLUMES_IDS = ( 1 << 3 );
inline constexpr size_type EXTRACTION_KEEP_ALL_IDS =
    EXTRACTION_KEEP_POINTS_IDS |
    EXTRACTION_KEEP_EDGES_IDS |
    EXTRACTION_KEEP_FACES_IDS |
    EXTRACTION_KEEP_VOLUMES_IDS;
inline constexpr size_type EXTRACTION_KEEP_MESH_RELATION = ( 1 << 4 );
inline constexpr size_type EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT =
    ( 1 << 5 );
} // namespace Feel

#endif // FEELPP_DISCR_SUBMESHCONTEXT_HPP
