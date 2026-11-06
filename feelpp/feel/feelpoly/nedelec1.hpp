/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file nedelec1.hpp
 * @brief Public descriptor for first-kind Nedelec finite elements.
 */
#ifndef FEELPP_FEELPOLY_NEDELEC1_HPP
#define FEELPP_FEELPOLY_NEDELEC1_HPP 1

#include <feel/feelpoly/nedelec.hpp>

namespace Feel
{
/**
 * @brief First-kind/trimmed Nedelec family descriptor.
 * @tparam Order public static order or `Dynamic`
 * @tparam Tag family tag retained for compatibility
 */
template<int Order, uint16_type Tag = 0>
using NedelecFirstKind = Nedelec<Order, NedelecKind::NED1, Tag>;
} // namespace Feel

#endif // FEELPP_FEELPOLY_NEDELEC1_HPP
