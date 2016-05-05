/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 15 Mar 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feel.hpp>
#include <iostream>

#include <feel/feelmodels/modelparameters.hpp>

namespace Feel
{

ModelParameters::ModelParameters( WorldComm const& world )
    : M_worldComm( world )
{
}

ModelParameters::ModelParameters( pt::ptree const& p, WorldComm const& world )
    : M_worldComm( world )
{
}

ModelParameters::~ModelParameters()
{
}

void ModelParameters::setPTree( pt::ptree const& p )
{
    M_p = p;
    setup();
}

void ModelParameters::setup()
{
    for ( auto const& v : M_p )
    {
        LOG( INFO ) << "reading parameter " << v.first << "\n";
        std::string t = v.first; // parameter name
        auto f = v.second;
        {
            try
            {
                auto val = M_p.get<double>( v.first );
                LOG( INFO ) << "adding parameter " << t << " with value " << val;
                this->operator[]( t ) = ModelParameter( t, val, val, val );
            }
            catch ( ... )
            {
                try
                {
                    auto val = M_p.get<std::string>( v.first );
                    LOG( INFO ) << "adding parameter " << t << " with value " << val;
                    this->operator[]( t ) = ModelParameter( t, val, M_directoryLibExpr, M_worldComm );
                }
                catch ( ... )
                {
                    try
                    {
                        auto val = f.get<double>( "value" );
                        LOG( INFO ) << "adding parameter " << t << " with value " << val;
                        this->operator[]( t ) = ModelParameter( t, val, val, val );
                    }
                    catch ( ... )
                    {
                        try
                        {
                            auto val = f.get<std::string>( "value" );
                            LOG( INFO ) << "adding parameter " << t << " with value " << val;
                            this->operator[]( t ) = ModelParameter( t, val, M_directoryLibExpr, M_worldComm );
                        }
                        catch ( ... )
                        {
                            this->operator[]( t ) = ModelParameter( t, 0., 0., 0. );
                        }
                    }
                }
            }
        }
    }
}

void ModelParameters::saveMD( std::ostream& os )
{
    auto myMap = toParameterValues();
    os << "### Parameters\n";
    os << "|Name|Value|Min|Max|\n";
    os << "|---|---|---|---|\n";
    for ( auto it = this->begin(); it != this->end(); it++ )
    {
        os << "|**" << it->first << "**|"
           << it->second.value() << "|"
           << it->second.min() << "|"
           << it->second.max() << "|\n";
    }
    os << "\n";
}

void ModelParameters::updateParameterValues()
{
#if 0
    for( auto const& p : *this )
        if ( p.second.hasExpression() )
            for ( auto const& mysymb : p.second.expression().expression().symbols() )
                std::cout << "p.first " << p.first << " mysymb " << mysymb.get_name() << "\n";
#endif
    std::map<std::string, double> mp;
    for ( auto const& p : *this )
        if ( !p.second.hasExpression() || p.second.expression().expression().symbols().empty() )
            mp[p.first] = p.second.value();
    // update all parameters with expression which not depends of symbols (can be improved)
    this->setParameterValues( mp );
}
void ModelParameters::setParameterValues( std::map<std::string, double> const& mp )
{
    for ( auto& p : *this )
        p.second.setParameterValues( mp );
}

std::map<std::string, double>
ModelParameters::toParameterValues() const
{
    std::map<std::string, double> pv;
    for ( auto const& p : *this )
        pv[p.first] = p.second.value();
    return pv;
}
}
