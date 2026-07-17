/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file geometryrepresentation.hpp
 * @brief Public umbrella for polynomial geometric-mapping representations.
 *
 * Include focused headers in implementation code when only one representation
 * or operation is required. This header preserves the complete public geometry
 * representation API for applications.
 */
#if !defined( FEELPP_FEELPOLY_GEOMETRYREPRESENTATION_HPP )
#define FEELPP_FEELPOLY_GEOMETRYREPRESENTATION_HPP 1

#include <feel/feelpoly/bernsteinbezier.hpp>
#include <feel/feelpoly/geometrycoefficients.hpp>
#include <feel/feelpoly/geometrymapping.hpp>
#include <feel/feelpoly/geometrynodallagrange.hpp>
#include <feel/feelpoly/geometryoperations.hpp>

#endif // FEELPP_FEELPOLY_GEOMETRYREPRESENTATION_HPP
