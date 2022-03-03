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
#include <boost/property_tree/json_parser.hpp>

#include <feel/feelmodels/modeloutputs.hpp>

namespace Feel {

ModelOutput::ModelOutput( worldcomm_ptr_t const& worldComm )
    :
    super( worldComm )
{}

ModelOutput::ModelOutput( std::string name, nl::json const& jarg, worldcomm_ptr_t const& worldComm, std::string const& directoryLibExpr )
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
    //else LOG(WARNING) << "Missing type entry in output " << M_name << "\n";

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

     if ( jarg.contains("coord") && jarg.at("coord").is_arraay() )
     {
         for( auto const& el : jarg.at("coord").items() )
         {
             if ( el.is_number() )
                 M_coord.push_back(el.get<double>());
             else if ( el.is_string() )
                 M_coord.push_back(std::stod( el.get<std::string>() ) );
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
}


std::string ModelOutput::getString( std::string const& key ) const
{
    if ( !M_p.contains("key") )
    {
        CHECK( false ) << "invalid key";
        return {};
    }

    auto const& j_key =  M_p.at("key");
    if ( !j_key.is_string() )
    {
        CHECK( false ) << "invalid key";
        return {};
    }

    return j_key.get<std::string>();
}

std::ostream& operator<<( std::ostream& os, ModelOutput const& o )
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

ModelOutputs::ModelOutputs( worldcomm_ptr_t const& worldComm )
    :
    super( worldComm )
{}

ModelOutputs::ModelOutputs( nl::json const& jarg, worldcomm_ptr_t const& worldComm )
    :
    super( worldComm ),
    M_p( jarg )
{
    setup();
}

ModelOutput
ModelOutputs::loadOutput( std::string const& s, std::string const& name )
{
    nl::json j;
    std::ifstream i(s);
    i >> j;
    return this->getOutput( j, name );
}

void
ModelOutputs::setup()
{
    auto const& jarg = M_p;
    for ( auto const& [jargkey,jargval] : jarg.items() )
    {
        LOG(INFO) << "Output :" << jargkey;
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

ModelOutput
ModelOutputs::getOutput( nl::json const& jarg, std::string const& name )
{
    ModelOutput o( name, jarg, this->worldCommPtr(), M_directoryLibExpr );
    LOG(INFO) << "adding output " << o;
    return o;
}

} // namespace Feel
