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

ModelParameters::ModelParameters( pt::ptree const& p, worldcomm_ptr_t const& world )
    :
    super( world )
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
        ModelExpression modelexpr;
        if ( f.empty() ) // parameter is not a node but a leaf
        {
            if( boost::optional<double> d = f.get_value_optional<double>() )
            {
                VLOG(2) << "adding parameter " << t << " with value " << *d;
                modelexpr.setExprScalar( Feel::vf::expr( *d ) );
            }
            else if ( boost::optional<std::string> s = f.get_value_optional<std::string>() )
            {
                VLOG(2) << "adding parameter " << t << " with expr " << *s;
                modelexpr.setExpr( *s, this->worldComm(), M_directoryLibExpr );
            }
            if ( modelexpr.hasAtLeastOneExpr() )
                this->operator[](t) = ModelParameter( t, modelexpr );
        }
        else
        {
            std::string name = t;
            if( boost::optional<std::string> n = f.get_optional<std::string>("name") )
                name = *n;
            std::string type = "expression";
            if( boost::optional<std::string> n = f.get_optional<std::string>("type") )
                type = *n;
            std::string desc = "";
            if( boost::optional<std::string> n = f.get_optional<std::string>("description") )
                desc = *n;
            if ( type == "expression" || type == "value" )
            {
                double val = 0.;
                double min = 0.;
                double max = 0.;
                if( boost::optional<double> d = f.get_optional<double>("min") )
                    min = *d;
                if( boost::optional<double> d = f.get_optional<double>("max") )
                    max = *d;

                modelexpr.setExpr( "value", f, this->worldComm(), M_directoryLibExpr );
                this->operator[](t) = ModelParameter( t, modelexpr, min, max, desc );
            }
            else if ( type == "fit" )
            {
                std::string filename;
                if( boost::optional<std::string> d = f.get_optional<std::string>("filename") )
                    filename = *d;
                else
                    CHECK( false ) << "filename is required with fit type";

                modelexpr.setExpr( "expr", f, this->worldComm(), M_directoryLibExpr );
                CHECK( modelexpr.hasAtLeastOneExpr() ) << "expr is required with fit type";

                std::string interpType = "Akima";
                if( boost::optional<std::string> d = f.get_optional<std::string>("interpolation") )
                    interpType = *d;

                auto itFindType = InterpolationTypeMap.find( interpType );
                CHECK( itFindType != InterpolationTypeMap.end() ) << "invalid interpolator type " << type;
                InterpolationType interpolatorEnumType = itFindType->second;

                std::string abscissa;
                if( boost::optional<std::string> d = f.get_optional<std::string>("abscissa") )
                    abscissa = *d;
                std::string ordinate;
                if( boost::optional<std::string> d = f.get_optional<std::string>("ordinate") )
                    ordinate = *d;

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
