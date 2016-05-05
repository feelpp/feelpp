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

namespace Feel
{
namespace vf
{

class GiNaCBase
{
  public:
    typedef double value_type;
    typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> vec_type;

    GiNaCBase() {}
    GiNaCBase( std::vector<GiNaC::symbol> const& syms )
        : M_syms( syms ),
          M_params( vec_type::Zero( M_syms.size() ) ),
          M_indexSymbolXYZ(),
          M_indexSymbolN()
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
            LOG( INFO ) << "index symbol relation  " << is.first << " and " << is.second << "\n";
        for ( auto const& is : M_indexSymbolN )
            LOG( INFO ) << "index symbol relation  " << is.first << " and " << is.second << "\n";

        this->setParameterFromOption();
    }
    GiNaCBase( GiNaCBase const& g )
        : M_syms( g.M_syms ),
          M_params( g.M_params ),
          M_indexSymbolXYZ( g.M_indexSymbolXYZ ),
          M_indexSymbolN( g.M_indexSymbolN )
    {
        this->setParameterFromOption();
    }
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

    std::set<std::pair<uint16_type, uint16_type>> const& indexSymbolXYZ() const { return M_indexSymbolXYZ; }
    std::set<std::pair<uint16_type, uint16_type>> const& indexSymbolN() const { return M_indexSymbolN; }

    void setParameterFromOption()
    {
        using namespace GiNaC;
        std::map<std::string, value_type> m;
        for ( auto const& s : M_syms )
        {
            if ( Environment::vm().count( s.get_name() ) )
            {
                // use try/catch in order to catch casting exception for
                // option that do not return double. Indeed we are only
                // collecting symbols in option database which can be cast
                // to numerical types
                try
                {
                    value_type v = option( _name = s.get_name() ).as<double>();
                    m.insert( std::make_pair( s.get_name(), v ) );
                    LOG( INFO ) << "symbol " << s.get_name() << " found in option with value " << v;
                }
                catch ( ... )
                {
                }

                //                    try
                //                    {
                //                        expression_type e( soption( _name=s.get_name() ), 0 );
                //                        if( is_a<numeric>(e) )
                //                        {
                //                            LOG(INFO) << "symbol " << s.get_name() << " found in option with value " << v;
                //                        }
                //                        else
                //                        {
                //                            ;
                //                        }
                //                    }
                //                    catch(...)
                //                    {}
            }
        }
        this->setParameterValues( m );
    }

    void setParameterValues( vec_type const& p )
    {
        CHECK( M_params.size() == M_syms.size() ) << "Invalid number of parameters " << M_params.size() << " >= symbol size : " << M_syms.size();
        M_params = p;
    }
    void setParameterValues( std::map<std::string, value_type> const& mp )
    {
        CHECK( M_params.size() == M_syms.size() ) << "Invalid number of parameters " << M_params.size() << " >= symbol size : " << M_syms.size();
        for ( auto const& p : mp )
        {
            auto it = std::find_if( M_syms.begin(), M_syms.end(),
                                    [&p]( GiNaC::symbol const& s ) { return s.get_name() == p.first; } );
            if ( it != M_syms.end() )
            {
                M_params[it - M_syms.begin()] = p.second;
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

  protected:
    std::vector<GiNaC::symbol> M_syms;
    vec_type M_params;
    std::set<std::pair<uint16_type, uint16_type>> M_indexSymbolXYZ;
    std::set<std::pair<uint16_type, uint16_type>> M_indexSymbolN;
};
}
} // vf / Feel

#endif /* __GiNaCBase_H */
