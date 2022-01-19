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

ModelParameters::ModelParameters( worldcomm_ptr_t const& world )
    :
    super( world )
{}

ModelParameters::~ModelParameters()
{}

void
ModelParameters::setPTree( nl::json const& jarg )
{
    M_p = jarg;
    setup();
}


void
ModelParameters::setup()
{
    for ( auto const& [jargkey,jargval] : M_p.items() )
    {
        std::string const& t = jargkey; // parameter name
        LOG(INFO) << "reading parameter " << t;

        if ( jargval.is_primitive() )
        {
            ModelExpression modelexpr;
            modelexpr.setExpr( jargval, this->worldComm(), M_directoryLibExpr );
            if ( modelexpr.hasAtLeastOneExpr() )
                this->emplace( std::make_pair( t, ModelParameter( t, modelexpr ) ) );
        }
        else if ( jargval.is_object() )
        {
            std::string name = t;
            if ( jargval.contains("name") )
            {
                auto const& j_name = jargval.at("name");
                if ( j_name.is_string() )
                    name = j_name.get<std::string>();
            }
            std::string type = "expression";
            if ( jargval.contains("type") )
            {
                auto const& j_type = jargval.at("type");
                if ( j_type.is_string() )
                    type = j_type.get<std::string>();
            }
            std::string desc;
            if ( jargval.contains("description") )
            {
                auto const& j_description = jargval.at("description");
                if ( j_description.is_string() )
                    desc = j_description.get<std::string>();
            }

            if ( type == "expression" || type == "value" )
            {
                double val = 0.;
                double min = 0.;
                double max = 0.;
                if ( jargval.contains("min") )
                {
                    auto const& j_min = jargval.at("min");
                    if ( j_min.is_number() )
                        min = j_min.get<double>();
                }
                if ( jargval.contains("max") )
                {
                    auto const& j_max = jargval.at("max");
                    if ( j_max.is_number() )
                        max = j_max.get<double>();
                }

                ModelExpression modelexpr;
                if ( jargval.contains("value") )
                    modelexpr.setExpr( jargval.at("value"), this->worldComm(), M_directoryLibExpr );
                this->operator[](t) = ModelParameter( t, modelexpr, min, max, desc );
            }
            else if ( type == "fit" )
            {
                std::string filename;
                if ( jargval.contains("filename") )
                    filename = jargval.at("filename").get<std::string>();
                CHECK( !filename.empty() ) << "filename is required with fit type";

                ModelExpression modelexpr;
                if ( jargval.contains("expr") )
                    modelexpr.setExpr( jargval.at("expr"), this->worldComm(), M_directoryLibExpr );
                CHECK( modelexpr.hasAtLeastOneExpr() ) << "expr is required with fit type";

                std::string interpType = "Akima";
                if ( jargval.contains("interpolation") )
                    interpType = jargval.at("interpolation").get<std::string>();

                auto itFindType = InterpolationTypeMap.find( interpType );
                CHECK( itFindType != InterpolationTypeMap.end() ) << "invalid interpolator type " << type;
                InterpolationType interpolatorEnumType = itFindType->second;

                std::string abscissa;
                if ( jargval.contains("abscissa") )
                    abscissa = jargval.at("abscissa").get<std::string>();

                std::string ordinate;
                if ( jargval.contains("ordinate") )
                    ordinate = jargval.at("ordinate").get<std::string>();

                std::shared_ptr<Interpolator> interpolator = Interpolator::New( /*interpType*/interpolatorEnumType, filename, abscissa, ordinate, this->worldComm() );
                this->operator[](t) = ModelParameter( name, interpolator,  modelexpr, desc );
            }
        }
    }
    this->updateParameterValues();
}

void
ModelParameters::saveMD(std::ostream &os)
{
#if 0
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
#endif
}


void
ModelParameters::updateParameterValues()
{
    // erase parameter in expr
    std::set<std::string> allSymbNames;
    for( auto const& p : *this )
        p.second.updateSymbolNames( allSymbNames );
    for( auto & p : *this )
        p.second.eraseParameterValues( allSymbNames );


    int previousParam = 0;
    while ( true )
    {
        std::map<std::string,double> mp = this->toParameterValues();
        if ( mp.size() == previousParam )
            break;
        previousParam = mp.size();
        this->setParameterValues( mp );
    }
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
    {
        auto const& mparam = p.second;
        if ( !mparam.isEvaluable() )
            continue;

        p.second.updateParameterValues( pv );
    }
    return pv;
}

void
ModelParameter::updateInformationObject( std::string const& symbol, nl::json::array_t & ja ) const
{
    if ( this->type() == "expression" )
    {
        auto [exprStr,compInfo] = M_expr.exprInformations();
        ja.push_back( symbolExprInformations( symbol, exprStr, compInfo, this->name() ) );
    }
    else if ( this->type() == "fit" )
    {
        auto compInfo = SymbolExprComponentSuffix( 1, 1 );
        ja.push_back( symbolExprInformations( symbol, "fit(.)", compInfo, this->name() ) );
    }
}

void
ModelParameters::updateInformationObject( nl::json & p ) const
{
    nl::json::array_t ja;

    for( auto const& [symbName, mparam] : *this )
        mparam.updateInformationObject( symbName,ja );

    p = ja;
}

}
