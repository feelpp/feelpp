/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#ifndef FEELPP_FEELDISCR_BROKENH_HPP
#define FEELPP_FEELDISCR_BROKENH_HPP 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/broken.hpp>

namespace Feel
{
/** @brief Function-space type for a broken finite-element family. */
template<class Family, class MeshType, typename T = double>
using brokenh_type = FunctionSpace<MeshType, bases<Broken<Family>>, T>;

/** @brief Shared pointer to a broken finite-element function space. */
template<class Family, class MeshType, typename T = double>
using brokenh_ptrtype = std::shared_ptr<brokenh_type<Family, MeshType, T>>;

/**
 * @brief Build a static-order broken function space.
 * @tparam Family wrapped FE family descriptor
 */
template<class Family, class MeshType, typename T = double>
    requires ( !Family::is_order_dynamic )
[[nodiscard]] brokenh_ptrtype<Family, MeshType, T>
brokenh( std::shared_ptr<MeshType> const& mesh,
             DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return brokenh_type<Family, MeshType, T>::New(
        _mesh = mesh,
        _worldscomm = makeWorldsComm( 1, mesh->worldComm() ),
        _extended_doftable = dte );
}

/**
 * @brief Build a runtime-order broken function space.
 * @tparam Family wrapped FE family descriptor, normally using `Dynamic`
 */
template<class Family, class MeshType, typename T = double>
[[nodiscard]] brokenh_ptrtype<Family, MeshType, T>
brokenh( std::shared_ptr<MeshType> const& mesh,
             RuntimeOrder order,
             DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return brokenh_type<Family, MeshType, T>::New(
        _mesh = mesh,
        _worldscomm = makeWorldsComm( 1, mesh->worldComm() ),
        _extended_doftable = dte,
        _runtime_order = order );
}
} // namespace Feel
#endif
