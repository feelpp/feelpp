/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 5 September 2019

 Copyright (C) 2019 Feel++ Consortium

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

#include <feel/feelmodels/modelinitialconditions.hpp>

namespace Feel {

void
ModelInitialCondition::setup( pt::ptree const& p, std::string const& typeIC )
{
    M_p = p;
    if ( typeIC == "Expression" )
    {
        M_modelExpr.setExpr( "expr", p, this->worldComm(), M_directoryLibExpr );
        CHECK( M_modelExpr.hasAtLeastOneExpr() ) << "the entry expr is required";
        M_isExpression = true;
        M_isFile = false;
        if ( auto ptmarkers = p.get_child_optional("markers") )
            M_markers.setPTree(*ptmarkers);
    }
    else if ( typeIC == "File" )
    {
        if ( auto itFileName = p.get_optional<std::string>("filename") )
            M_fileName = *itFileName;
        if ( auto itFormat = p.get_optional<std::string>("format") )
            M_fileType = *itFormat;
        if ( auto itDir = p.get_optional<std::string>("directory") )
            M_fileDirectory = *itDir;
        M_isExpression = false;
        M_isFile = true;
    }
    else // default case: we just store the ptree
    {
        M_isExpression = false;
        M_isFile = false;
    }

}

void
ModelInitialCondition::setParameterValues( std::map<std::string,double> const& mp )
{
    if ( this->isExpression() )
        M_modelExpr.setParameterValues( mp );
}



void
ModelInitialConditionTimeSet::setup( pt::ptree const& p )
{
    M_p = p;
    for ( auto const& icAllExprPtree : p )
    {
        std::string const typeIC = icAllExprPtree.first; // typeIC = Expression, File, ...
        for ( auto const& icExprPtree : icAllExprPtree.second )
        {
            std::string icName = icExprPtree.first;
            double time = 0;
            if( boost::optional<double> ittime = icExprPtree.second.get_optional<double>( "time" ) )
                time = *ittime;

            ModelInitialCondition mic( this->worldCommPtr(), M_directoryLibExpr );
            mic.setName( icName );
            mic.setup( icExprPtree.second, typeIC );
            this->operator[]( time )[typeIC].push_back( mic );
        }
    }
}

void
ModelInitialConditionTimeSet::setParameterValues( std::map<std::string,double> const& mp )
{
    for ( auto & [time,icByType] : *this )
        for (auto & [type,vecOfIC] : icByType )
            for (auto & ic : vecOfIC )
                ic.setParameterValues( mp );
}


void
ModelInitialConditions::setPTree( pt::ptree const& p )
{
    for( auto const& subPtree : p )
    {
        std::string name = subPtree.first;
        if ( this->find( name ) != this->end() )
        {
            this->operator[]( name ).setup( subPtree.second );
        }
        else
        {
            ModelInitialConditionTimeSet micts( this->worldCommPtr(), M_directoryLibExpr );
            micts.setup( subPtree.second );
            this->operator[]( name ) = micts;
        }
    }
}

void
ModelInitialConditions::setParameterValues( std::map<std::string,double> const& mp )
{
    for ( auto & [name,ic] : *this )
        ic.setParameterValues( mp );
}

ModelInitialConditionTimeSet const&
ModelInitialConditions::get( std::string const& name ) const
{
    auto itFindIC = this->find( name );
    if ( itFindIC != this->end() )
        return itFindIC->second;
    return M_emptylInitialConditionTimeSet;
}

} // namespace Feel
