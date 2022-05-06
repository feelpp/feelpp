/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2021-09-15

  Copyright (C) 2021 Feel++ Consortium

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
#include <fmt/ranges.h>
#include <feel/feelvf/ginacbase.hpp>

namespace Feel::vf
{

GiNaCBase::GiNaCBase( std::vector<GiNaC::symbol> const& syms )
: M_syms( syms ),
  M_params( vec_type::Zero( M_syms.size() ) ),
  M_indexSymbolXYZ(),
  M_indexSymbolN(),
  M_indexSymbolGeom(),
  M_context( 0 ),
  M_isNumericExpression( false )
{
        // detect if symbol x,y,z are present and get index access in M_params
        std::map<int, std::string> lstxyz{{0, "x"}, {1, "y"}, {2, "z"}};
        for ( auto const& str : lstxyz )
        {
            auto itSym = std::find_if( M_syms.begin(), M_syms.end(),
                                       [&str]( GiNaC::symbol const& s ) { return s.get_name() == str.second; } );
            if ( itSym != M_syms.end() )
                M_indexSymbolXYZ.insert( std::make_pair( str.first, std::distance( M_syms.begin(), itSym ) ) );
        }
        std::map<int, std::string> lstN{{3, "nx"}, {4, "ny"}, {5, "nz"}};
        for ( auto const& str : lstN )
        {
            auto itSym = std::find_if( M_syms.begin(), M_syms.end(),
                                       [&str]( GiNaC::symbol const& s ) { return s.get_name() == str.second; } );
            if ( itSym != M_syms.end() )
                M_indexSymbolN.insert( std::make_pair( str.first, std::distance( M_syms.begin(), itSym ) ) );
        }
        std::map<int, std::string> lstGeom{{6, "h"}, {7, "meas"}, {8, "measPEN"}, {9, "nPEN"}, {10, "emarker"}};
        for ( auto const& str : lstGeom )
        {
            auto itSym = std::find_if( M_syms.begin(), M_syms.end(),
                                       [&str]( GiNaC::symbol const& s ) { return s.get_name() == str.second; } );
            if ( itSym != M_syms.end() )
            {
                int d = std::distance( M_syms.begin(), itSym );
                M_indexSymbolGeom.insert( std::make_pair( str.first, d ) );
                VLOG( 1 ) << fmt::format( "[expr relation] {}:{}", str.first, d ) << std::endl;
            }
        }
        for ( auto const& is : M_indexSymbolXYZ )
        {
            VLOG( 1 ) << "index symbol relation:  " << is.first << " -> " << is.second << "\n";
        }
        for ( auto const& is : M_indexSymbolN )
            VLOG( 1 ) << "index symbol relation:  " << is.first << " -> " << is.second << "\n";

        if ( hasSymbol( "x" ) || hasSymbol( "y" ) || hasSymbol( "z" ) )
            M_context = M_context | vm::POINT;
        if ( hasAnySymbolN() )
            M_context = M_context | vm::KB | vm::NORMAL;
        if ( hasSymbol( "h" ) | hasSymbol( "measPEN" ) || hasSymbol( "meas" ) )
            M_context = M_context | vm::MEASURE;
}
} // Feel::vf