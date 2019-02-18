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
#include <iostream>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>

#include <feel/feelmodels/modelparameters.hpp>
#include <boost/optional.hpp>

namespace Feel {

ModelParameters::ModelParameters( WorldComm const& world )
    :
    M_worldComm( world )
{}

ModelParameters::ModelParameters( pt::ptree const& p, WorldComm const& world )
    :
    M_worldComm( world )
{}

ModelParameters::~ModelParameters()
{}

void
ModelParameters::setPTree( pt::ptree const& p )
{
    M_p = p;
    setup();
}


void
ModelParameters::setup()
{
    for( auto const& v : M_p )
    {
        LOG(INFO) << "reading parameter " << v.first  << "\n";
        std::string t = v.first; // parameter name
        auto f= v.second;
        if ( f.empty() ) // parameter is not a node but a leaf
        {
            if( boost::optional<double> d = f.get_value_optional<double>() )
            {
                LOG(INFO) << "adding parameter " << t << " with value " << *d;
                this->operator[](t) = ModelParameter( t, *d, 0, 0 );
            }
            else if ( boost::optional<std::string> s = f.get_value_optional<std::string>() )
            {
                LOG(INFO) << "adding parameter " << t << " with expr " << *s;
                this->operator[](t) = ModelParameter( t, *s, M_directoryLibExpr, M_worldComm );
            }
            else
            {
                LOG(INFO) << "adding parameter " << t << " with default value 0";
                this->operator[](t) = ModelParameter( t, 0, 0, 0 );
            }
        }
        else
        {
            std::string name = t;
            if( boost::optional<std::string> n = f.get_optional<std::string>("name") )
                name = *n;
            std::string type = "expression";
            if( boost::optional<std::string> n = f.get_optional<std::string>("type") )
                type = *n;
            if ( type == "expression" || type == "value" )
            {
                double val = 0.;
                double min = 0.;
                double max = 0.;
                if( boost::optional<double> d = f.get_optional<double>("min") )
                    min = *d;
                if( boost::optional<double> d = f.get_optional<double>("max") )
                    max = *d;
                if( boost::optional<double> d = f.get_optional<double>("value") )
                    this->operator[](t) = ModelParameter( name, *d, min, max );
                else if ( boost::optional<std::string> s = f.get_optional<std::string>("value") )
                    this->operator[](t) = ModelParameter( name, *s, M_directoryLibExpr, M_worldComm, min, max );
                else
                    this->operator[](t) = ModelParameter( name, 0., min, max );
            }
            else if ( type == "fit" )
            {
                std::string filename;
                if( boost::optional<std::string> d = f.get_optional<std::string>("filename") )
                    filename = *d;
                else
                    CHECK( false ) << "filename is required with fit type";
                std::string exprStr;
                if( boost::optional<std::string> d = f.get_optional<std::string>("expr") )
                    exprStr = *d;
                else
                    CHECK( false ) << "expr is required with fit type";
                std::string interpType = "Akima";
                if( boost::optional<std::string> d = f.get_optional<std::string>("interpolation") )
                    interpType = *d;

                auto itFindType = InterpolationTypeMap.find( interpType );
                CHECK( itFindType != InterpolationTypeMap.end() ) << "invalid interpolator type " << type;
                InterpolationType interpolatorEnumType = itFindType->second;


                std::shared_ptr<Interpolator> interpolator = Interpolator::New( /*interpType*/interpolatorEnumType, filename, M_worldComm );
                this->operator[](t) = ModelParameter( name, interpolator, exprStr, M_directoryLibExpr, M_worldComm );
            }
        }
    }
}

void
ModelParameters::saveMD(std::ostream &os)
{
  auto myMap = toParameterValues();
  os << "### Parameters\n";
  os << "|Name|Value|Min|Max|\n";
  os << "|---|---|---|---|\n";
  for(auto it =this->begin(); it != this->end(); it++)
  {
    os << "|**" << it->first << "**|" 
      << it->second.value() << "|" 
      << it->second.min() << "|" 
      << it->second.max() << "|\n";
  }
  os << "\n";
}


void
ModelParameters::updateParameterValues()
{
#if 0
    for( auto const& p : *this )
        if ( p.second.hasExpression() )
            for ( auto const& mysymb : p.second.expression().expression().symbols() )
                std::cout << "p.first " << p.first << " mysymb " << mysymb.get_name() << "\n";
#endif
    std::map<std::string,double> mp;
    for( auto const& p : *this )
        if ( !p.second.hasExpression() || p.second.expression().expression().symbols().empty() )
            mp[p.first] = p.second.value();
    // update all parameters with expression which not depends of symbols (can be improved)
    this->setParameterValues( mp );
}
void
ModelParameters::setParameterValues( std::map<std::string,double> const& mp )
{
    for( auto & p : *this )
        p.second.setParameterValues( mp );
}

std::map<std::string,double>
ModelParameters::toParameterValues() const
{
    std::map<std::string,double> pv;
    for( auto const& p : *this )
        if ( p.second.type() != "fit" )
            pv[p.first]=p.second.value();
    return pv;
}



}
