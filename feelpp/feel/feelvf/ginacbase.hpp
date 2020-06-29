/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-03-31

  Copyright (C) 2014-2016 Feel++ Consortium

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
/**
   \file ginacbase.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-03-31
 */
#ifndef FEELPP_GINACBASE_HPP
#define FEELPP_GINACBASE_HPP

namespace Feel::vf
{

//class ExprDynamicBase;

class GiNaCBase : public Feel::vf::ExprDynamicBase
{
  public:
    typedef double value_type;
    typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> vec_type;

    GiNaCBase() = default;
    GiNaCBase( std::vector<GiNaC::symbol> const& syms )
        : M_syms( syms ),
          M_params( vec_type::Zero( M_syms.size() ) ),
          M_indexSymbolXYZ(),
          M_indexSymbolN(),
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
    }

    GiNaCBase( GiNaCBase const& g ) = default;

    virtual ~GiNaCBase() {}

    const std::vector<GiNaC::symbol>& symbols() const
    {
        return M_syms;
    }

    bool hasSymbol( std::string const& symb ) const
    {
        for ( auto const& s : M_syms )
            if ( s.get_name() == symb )
                return true;
        return false;
    }

    vec_type const& parameterValue() const { return M_params; }
    value_type parameterValue( int p ) const { return M_params[p]; }

    std::map<std::string, value_type> const& symbolNameToValue() const { return M_symbolNameToValue; }

    std::set<std::pair<uint16_type, uint16_type>> const& indexSymbolXYZ() const { return M_indexSymbolXYZ; }
    std::set<std::pair<uint16_type, uint16_type>> const& indexSymbolN() const { return M_indexSymbolN; }
    bool hasAnySymbolXYZ() const { return !M_indexSymbolXYZ.empty(); }
    bool hasAnySymbolN() const { return !M_indexSymbolN.empty(); }
    size_type dynamicContext() const
    {
        return M_context;
    }

    //! return true if the expression is numeric (i.e. not a symbolic expression)
    bool isNumericExpression() const { return M_isNumericExpression; }

    uint16_type index( std::string const& sname ) const
        {
            auto it = std::find_if( M_syms.begin(), M_syms.end(),
                                    [=]( GiNaC::symbol const& s ) { return s.get_name() == sname; } );
            if ( it != M_syms.end() )
            {
                return it-M_syms.begin();
            }
            return invalid_uint16_type_value;
        }

    void setParameterValues( vec_type const& p )
    {
        CHECK( M_params.size() == M_syms.size() ) << "Invalid number of parameters " << M_params.size() << " >= symbol size : " << M_syms.size();
        M_params = p;
    }
    virtual void setParameterValues( std::map<std::string, value_type> const& mp )
    {
        CHECK( M_params.size() == M_syms.size() ) << "Invalid number of parameters " << M_params.size() << " >= symbol size : " << M_syms.size();
        for ( auto const& p : mp )
        {
            auto it = std::find_if( M_syms.begin(), M_syms.end(),
                                    [&p]( GiNaC::symbol const& s ) { return s.get_name() == p.first; } );
            if ( it != M_syms.end() )
            {
                M_params[it - M_syms.begin()] = p.second;
                M_symbolNameToValue[p.first] = p.second;
                VLOG( 2 ) << "setting parameter : " << p.first << " with value: " << p.second;
                VLOG( 2 ) << "parameter: \n"
                          << M_params;
            }
            else
            {
                VLOG( 1 ) << "Invalid parameters : " << p.first << " with value: " << p.second;
            }
        }
    }

    void eraseParameterValues( std::set<std::string> const& symbNames )
        {
            for ( std::string const& symbName : symbNames )
                M_symbolNameToValue.erase( symbName );
        }


  protected:
    std::vector<GiNaC::symbol> M_syms;
    vec_type M_params;
    std::set<std::pair<uint16_type, uint16_type>> M_indexSymbolXYZ;
    std::set<std::pair<uint16_type, uint16_type>> M_indexSymbolN;
    std::map<std::string, value_type> M_symbolNameToValue;
    size_type M_context;
    bool M_isNumericExpression;
};

} // namespace Feel::vf

#endif /* __GiNaCBase_H */
