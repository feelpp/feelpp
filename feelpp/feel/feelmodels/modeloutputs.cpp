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

ModelOutput::ModelOutput( WorldComm const& worldComm )
    :
    M_worldComm( &worldComm )
{}

ModelOutput::ModelOutput( std::string name, pt::ptree const& p, WorldComm const& worldComm, std::string const& directoryLibExpr )
    :
    M_worldComm( &worldComm ),
    M_p( p ),
    M_directoryLibExpr( directoryLibExpr ),
    M_name( name )
{
    try {
        M_type  = M_p.get<std::string>( "type"  );
    }
    catch ( pt::ptree_bad_path& e )
    {
        LOG(WARNING) << "Missing type entry in output " << M_name << "\n";
    }
    if( auto range = p.get_child_optional("range") )
    {
        for( auto const& item : p.get_child("range") )
            M_range.insert(item.second.template get_value<std::string>());
        if( M_range.empty() )
            M_range.insert(p.get<std::string>("range") );
    }
    else
    {
        LOG(WARNING) << "Output " << M_name << " does not have any range\n";
    }
    if( auto t = p.get_optional<int>("topodim") )
    {
        M_dim = *t;
    }
    else
    {
        LOG(WARNING) << "Output " << M_name << " does not have any dimension\n";
    }
}

std::string ModelOutput::getString( std::string const& key ) const
{
    try { return M_p.get<std::string>( key ); }
    catch( pt::ptree_error const&e ) {
        LOG(ERROR) << "output " << M_name << ": key " << key << ": " << e.what() << std::endl;
        exit(1);
    }
}

std::ostream& operator<<( std::ostream& os, ModelOutput const& o )
{
    os << "Output " << o.name()
       << "[ type: " << o.type()
       << ", dim: " << o.dim()
       << ", range: [";
    for( auto const& r : o.range() )
    {
        os << r << " ";
    }
    os << "]]" << std::endl;
    return os;
}

ModelOutputs::ModelOutputs( WorldComm const& worldComm )
    :
    M_worldComm( &worldComm )
{}

ModelOutputs::ModelOutputs( pt::ptree const& p, WorldComm const& worldComm )
    :
    M_worldComm( &worldComm ),
    M_p( p )
{
    setup();
}

ModelOutput
ModelOutputs::loadOutput( std::string const& s, std::string const& name )
{
    pt::ptree p;
    pt::read_json( s, p );
    return this->getOutput( p, name );
}

void
ModelOutputs::setup()
{
    for( auto const& v : M_p )
    {
        LOG(INFO) << "Output :" << v.first  << "\n";
        if ( auto fname = v.second.get_optional<std::string>("filename") )
        {
            LOG(INFO) << "  - filename = " << Environment::expand( fname.get() ) << std::endl;
            this->insert( std::make_pair( v.first, this->loadOutput( Environment::expand( fname.get() ), v.first ) ) );
        }
        else
        {
            this->insert( std::make_pair( v.first, this->getOutput( v.second, v.first ) ) );
        }
    }
}

ModelOutput
ModelOutputs::getOutput( pt::ptree const& v, std::string const& name )
{
    ModelOutput o( name, v, *M_worldComm, M_directoryLibExpr );
    LOG(INFO) << "adding output " << o;
    return o;
}

} // namespace Feel
