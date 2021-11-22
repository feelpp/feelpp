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

#include <Eigen/Core>
#include <feel/feelvf/exprbase.hpp>
#include <ginac/ginac.h>
namespace Feel::vf
{

//class ExprDynamicBase;

class GiNaCBase : public Feel::vf::ExprDynamicBase
{
  public:
    typedef double value_type;
    typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> vec_type;

    GiNaCBase() = default;
    GiNaCBase( std::vector<GiNaC::symbol> const& syms );
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
    std::set<std::pair<uint16_type, uint16_type>> const& indexSymbolGeom() const { return M_indexSymbolGeom; }
    bool hasAnySymbolXYZ() const { return !M_indexSymbolXYZ.empty(); }
    bool hasAnySymbolN() const { return !M_indexSymbolN.empty(); }
    bool hasAnySymbolGeom() const { return !M_indexSymbolGeom.empty(); }
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
    std::set<std::pair<uint16_type, uint16_type>> M_indexSymbolGeom;
    std::map<std::string, value_type> M_symbolNameToValue;
    size_type M_context;
    bool M_isNumericExpression;
};

} // namespace Feel::vf

#endif /* __GiNaCBase_H */
