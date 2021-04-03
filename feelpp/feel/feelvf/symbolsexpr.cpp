/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2018-06-03

  Copyright (C) 2018 Feel++ Consortium

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

#include <feel/feelvf/symbolsexpr.hpp>

namespace Feel
{
namespace vf
{
nl::json symbolExprInformations( std::string const& symbol, std::string const& expr, SymbolExprComponentSuffix const& symbolSuffix, std::string const& name )
{
    nl::json jo( { { "symbol", symbol } } );
    if ( !name.empty() )
         jo.emplace( "name",name );
    if ( !expr.empty() )
         jo.emplace( "expr",expr );

    if ( symbolSuffix.empty() )
    {
        jo.emplace( "shape",  "scalar" );
    }
    else
    {
        std::string shape;
        if ( symbolSuffix.nComp2() == 1 )
            shape = "vectorial [" + std::to_string( symbolSuffix.nComp1() ) + "]";
        else if ( symbolSuffix.nComp1() == 1 )
            shape = "vectorial [" + std::to_string( symbolSuffix.nComp2() ) + "]";
        else
            shape = "tensor2 [" + std::to_string( symbolSuffix.nComp2() ) + "x" + std::to_string( symbolSuffix.nComp2() ) + "]";

        nl::json::array_t jcomp;
        for ( auto const& [_suffix,compArray] : symbolSuffix )
            jcomp.push_back( nl::json( { { "symbol", symbol+ _suffix }, { "indices", json::array({compArray[0],compArray[1]}) } } ) );

        jo.emplace( "shape", shape );
        jo.emplace( "components", jcomp );
    }
    return jo;
}

} // namespace vf
} // namespace Feel
