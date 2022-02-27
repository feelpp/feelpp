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
ModelInitialCondition::setup( nl::json const& jarg, std::string const& typeIC )
{
    if ( typeIC == "Expression" )
    {
        CHECK( jarg.contains("expr") ) << "an expression is required";
        M_modelExpr.setExpr( jarg.at("expr"), this->worldComm(), M_directoryLibExpr );
        CHECK( M_modelExpr.hasAtLeastOneExpr() ) << "the entry expr is required";
        M_isExpression = true;
        M_isFile = false;
        if ( jarg.contains("markers") )
            M_markers.setup( jarg.at("markers") );
    }
    else if ( typeIC == "File" )
    {
        if ( jarg.contains("filename") )
            M_fileName = jarg.at("filename").get<std::string>();
        if ( jarg.contains("format") )
            M_fileType = jarg.at("format").get<std::string>();
        if ( jarg.contains("directory") )
            M_fileDirectory = jarg.at("directory").get<std::string>();
        M_isExpression = false;
        M_isFile = true;
    }
    else // default case: we just store the ptree
    {
        M_isExpression = false;
        M_isFile = false;
        M_jsonSetup = jarg;
    }
}

void
ModelInitialCondition::setParameterValues( std::map<std::string,double> const& mp )
{
    if ( this->isExpression() )
        M_modelExpr.setParameterValues( mp );
}



void
ModelInitialConditionTimeSet::setup( nl::json const& jarg )
{
    //M_p = p;
    for ( auto const& [jargkey,jargval] : jarg.items() )
    {
        std::string const typeIC = jargkey; // typeIC = Expression, File, ...
        for ( auto const& [j_typekey,j_typeval] : jargval.items() )
        {
            std::string icName = j_typekey;
            double time = 0;
            if ( j_typeval.contains("time") )
            {
                auto const& j_type_time = j_typeval.at("time");
                if ( j_type_time.is_number() )
                    time = j_type_time.get<double>();
                else if ( j_type_time.is_string() )
                    time = std::stod( j_type_time.get<std::string>() );
            }

            ModelInitialCondition mic( this->worldCommPtr(), M_directoryLibExpr );
            mic.setName( icName );
            mic.setup( j_typeval, typeIC );
            this->operator[]( time )[typeIC].push_back( std::move( mic ) );
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
ModelInitialConditions::setPTree( nl::json const& jarg )
{
    bool useToolboxName = true;
    if ( useToolboxName )
    {
        for ( auto const& [jargkey,jargval] : jarg.items() )
            this->setupInternal(jargkey,jargval );
    }
    else
    {
        this->setupInternal( "", jarg );
    }
}
void
ModelInitialConditions::setupInternal( std::string const& tbname, nl::json const& jarg )
{
    for ( auto const& [jargkey,jargval] : jarg.items() )
    {
        std::string const& fieldname = jargkey;
        if ( this->has(tbname,fieldname ) )
        {
            this->operator[](tbname)[fieldname].setup( jargval );
        }
        else
        {
            ModelInitialConditionTimeSet micts( this->worldCommPtr(), M_directoryLibExpr );
            micts.setup( jargval );
            this->operator[]( tbname ).emplace( std::make_pair( fieldname, std::move( micts ) ) );
        }
    }
}
void
ModelInitialConditions::setParameterValues( std::map<std::string,double> const& mp )
{
    for ( auto & [tbname,icdata] : *this )
        for ( auto & [name,ic] : icdata )
            ic.setParameterValues( mp );
}

bool
ModelInitialConditions::has( std::string const& tbname, std::string const& fieldname ) const
{
    auto itFindTB = this->find( tbname );
    if ( itFindTB == this->end() )
        return false;
    auto itFindField = itFindTB->second.find( fieldname );
    if ( itFindField == itFindTB->second.end() )
        return false;
    return true;
}
ModelInitialConditionTimeSet const&
ModelInitialConditions::get( std::string const& tbname, std::string const& fieldname ) const
{
    if ( !this->has( tbname,fieldname ) )
        return M_emptylInitialConditionTimeSet;
    return this->find( tbname )->second.find( fieldname )->second;
}

} // namespace Feel
