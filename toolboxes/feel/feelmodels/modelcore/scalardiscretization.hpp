/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#pragma once

#include <feel/feelpoly/order.hpp>
#include <charconv>
#include <optional>
#include <regex>

namespace Feel::FeelModels
{
struct ScalarDiscretization
{
    int polynomialOrder = 1;
    int geometryOrder = 1;
};

// Scalar toolbox CLI syntax: P<k> or P<k>G<g>.
inline std::optional<ScalarDiscretization>
parseScalarDiscretization( std::string const& value )
{
    static const std::regex re( R"(^P([0-9]+)(?:G([0-9]+))?$)" );
    std::smatch match;
    if ( !std::regex_match( value, match, re ) )
        return std::nullopt;
    auto parseOrder = []( std::string const& token, int& order ) {
        auto result = std::from_chars( token.data(), token.data()+token.size(), order );
        return result.ec == std::errc{} && order >= 1 && order <= std::numeric_limits<uint16_type>::max();
    };
    ScalarDiscretization spec;
    if ( !parseOrder( match[1].str(), spec.polynomialOrder ) ||
         ( match[2].matched && !parseOrder( match[2].str(), spec.geometryOrder ) ) )
        return std::nullopt;
    return spec;
}
} // namespace Feel::FeelModels
