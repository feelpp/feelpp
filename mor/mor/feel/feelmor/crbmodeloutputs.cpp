/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Mar 2015
 
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

#include <feel/feelmor/crbmodeloutputs.hpp>

namespace Feel {

CRBModelOutput::CRBModelOutput( worldcomm_ptr_t const& worldComm )
    :
    super( worldComm )
{}

CRBModelOutput::CRBModelOutput( std::string name, nl::json const& jarg, worldcomm_ptr_t const& worldComm, std::string const& directoryLibExpr )
    :
    super( worldComm ),
    M_p( jarg ),
    M_directoryLibExpr( directoryLibExpr ),
    M_name( name ),
    M_markers( name )
{
    if ( jarg.contains("type") )
    {
        auto const& j_type = jarg.at("type");
        if ( j_type.is_string() )
            M_type = j_type.get<std::string>();
    }
    else LOG(WARNING) << "Missing type entry in output " << M_name << "\n";

     if ( jarg.contains("markers") )
         M_markers.setup( jarg.at("markers") );
     //else LOG(WARNING) << "Output " << M_name << " does not have any range\n";

     if ( jarg.contains("topodim") )
     {
         auto const& j_topodim = jarg.at("topodim");
         if ( j_topodim.is_number_unsigned() )
             M_dim = j_topodim.get<int>();
         else if ( j_topodim.is_string() )
             M_dim = std::stoi( j_topodim.get<std::string>() );
     }
     // else LOG(WARNING) << "Output " << M_name << " does not have any dimension\n";

     if ( jarg.contains("coord") && jarg.at("coord").is_array() )
     {
         for( auto const& [key,value] : jarg.at("coord").items() )
         {
             if ( value.is_number() )
                 M_coord.push_back(value.get<double>());
             else if ( value.is_string() )
                 M_coord.push_back(std::stod( value.get<std::string>() ) );
         }
     }

     if ( jarg.contains("radius") )
     {
         auto const& j_radius = jarg.at("radius");
         if ( j_radius.is_number() )
             M_radius = j_radius.get<double>();
         else if ( j_radius.is_string() )
             M_radius = std::stod( j_radius.get<std::string>() );
     }

    if ( jarg.contains("expr") )
        M_expr.setExpr( jarg.at("expr"), this->worldComm(), M_directoryLibExpr );
    else
        M_expr.setExpr( std::string("crb_u:crb_u"), this->worldComm(), M_directoryLibExpr );
}


std::string CRBModelOutput::getString( std::string const& key ) const
{
    if ( !M_p.contains(key) )
    {
        CHECK( false ) << "invalid key " << key;
        return {};
    }

    auto const& j_key =  M_p.at(key);
    if ( !j_key.is_string() )
    {
        CHECK( false ) << "key " << key << " not a string";
        return {};
    }

    return j_key.get<std::string>();
}

std::ostream& operator<<( std::ostream& os, CRBModelOutput const& o )
{
    os << "Output " << o.name()
       << "[ type: " << o.type()
       << ", dim: " << o.dim()
       << ", range: [";
    for( auto const& r : o.markers() )
    {
        os << r << " ";
    }
    os << "]]" << std::endl;
    return os;
}

CRBModelOutputs::CRBModelOutputs( worldcomm_ptr_t const& worldComm )
    :
    super( worldComm )
{}

CRBModelOutputs::CRBModelOutputs( nl::json const& jarg, worldcomm_ptr_t const& worldComm )
    :
    super( worldComm ),
    M_p( jarg )
{
    setup();
}

CRBModelOutputs::~CRBModelOutputs()
{}

void
CRBModelOutputs::setPTree( nl::json const& jarg )
{
    M_p = jarg;
    setup();
}



CRBModelOutput
CRBModelOutputs::loadOutput( std::string const& s, std::string const& name )
{
    nl::json j;
    std::ifstream i(s);
    i >> j;
    return this->getOutput( j, name );
}

void
CRBModelOutputs::setup()
{
    auto const& jarg = M_p;
    for ( auto const& [jargkey,jargval] : jarg.items() )
    {
        LOG(INFO) << "CRBOutput :" << jargkey;
        if ( jargval.contains("filename") )
        {
            auto const& j_filename = jargval.at("filename");
            std::string fname = Environment::expand( j_filename.get<std::string>() );
            LOG(INFO) << "  - filename = " << fname << std::endl;
            this->insert( std::make_pair( jargkey, this->loadOutput( fname, jargkey ) ) );

        }
        else
            this->insert( std::make_pair( jargkey, this->getOutput( jargval, jargkey ) ) );
    }
}

CRBModelOutput
CRBModelOutputs::getOutput( nl::json const& jarg, std::string const& name )
{
    CRBModelOutput o( name, jarg, this->worldCommPtr(), M_directoryLibExpr );
    LOG(INFO) << "adding crboutput " << o;
    return o;
}

void
CRBModelOutputs::setParameterValues( std::map<std::string,double> const& mp )
{
    for( auto & p : *this )
        p.second.setParameterValues( mp );
}

} // namespace Feel
