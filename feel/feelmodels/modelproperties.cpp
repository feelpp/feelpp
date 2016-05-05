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
#include <feel/feelcore/removecomments.hpp>
#include <feel/feelcore/utility.hpp>

#include <boost/property_tree/json_parser.hpp>
#include <feel/feelmodels/modelproperties.hpp>

namespace Feel
{

ModelProperties::ModelProperties( std::string const& filename, std::string const& directoryLibExpr, WorldComm const& world )
    : M_worldComm( world ),
      M_params( world ),
      M_bc( world ),
      M_postproc( world )
{
    if ( !fs::exists( filename ) )
    {
        LOG( INFO ) << "Could not find " << filename << std::endl;
        return;
    }
    else
    {
        LOG( INFO ) << "Loading " << filename << std::endl;
    }

    auto json_str_wo_comments = removeComments( readFromFile( filename ) );
    LOG( INFO ) << "json file without comment:" << json_str_wo_comments;

    std::istringstream istr( json_str_wo_comments );
    pt::read_json( istr, M_p );
    try
    {
        M_name = M_p.get<std::string>( "Name" );
    }
    catch ( pt::ptree_bad_path& e )
    {
        if ( Environment::isMasterRank() )
            std::cout << "Missing Name entry in model properties\n";
    }
    try
    {
        M_model = M_p.get<std::string>( "Model" );
    }
    catch ( pt::ptree_bad_path& e )
    {
        if ( Environment::isMasterRank() )
            std::cout << "Missing Model entry in model properties\n";
    }
    try
    {
        M_shortname = M_p.get( "ShortName", M_name );
    }
    catch ( pt::ptree_bad_path& e )
    {
        if ( Environment::isMasterRank() )
            std::cout << "Missing ShortName entry in model properties - set it to the Name entry : " << M_name << "\n";
    }
    auto par = M_p.get_child_optional( "Parameters" );
    if ( par )
    {
        LOG( INFO ) << "Model with parameters\n";
        if ( !directoryLibExpr.empty() )
            M_params.setDirectoryLibExpr( directoryLibExpr );
        M_params.setPTree( *par );
    }
    auto func = M_p.get_child_optional( "Functions" );
    if ( func )
    {
        LOG( INFO ) << "Model with functions\n";
        if ( !directoryLibExpr.empty() )
            M_functions.setDirectoryLibExpr( directoryLibExpr );
        M_functions.setPTree( *func );
    }
    auto bc = M_p.get_child_optional( "BoundaryConditions" );
    if ( bc )
    {
        LOG( INFO ) << "Model with boundary conditions\n";
        if ( !directoryLibExpr.empty() )
            M_bc.setDirectoryLibExpr( directoryLibExpr );
        M_bc.setPTree( *bc );
    }
    else
    {
        LOG( WARNING ) << "Model does not have any boundary conditions\n";
    }
    auto mat = M_p.get_child_optional( "Materials" );
    if ( mat )
    {
        LOG( INFO ) << "Model with materials\n";
        M_mat.setPTree( *mat );
    }
    else
    {
        LOG( WARNING ) << "Model does not have any materials\n";
    }
    auto pp = M_p.get_child_optional( "PostProcess" );
    if ( pp )
    {
        LOG( INFO ) << "Model with PostProcess\n";
        if ( !directoryLibExpr.empty() )
            M_postproc.setDirectoryLibExpr( directoryLibExpr );
        M_postproc.setPTree( *pp );
    }
    else
    {
        LOG( WARNING ) << "Model does not have any materials\n";
    }
}

ModelProperties::~ModelProperties()
{
}

std::string ModelProperties::getEntry( std::string& s )
{
    std::string res;
    try
    {
        res = M_p.get<std::string>( s );
    }
    catch ( pt::ptree_bad_path& e )
    {
        if ( Environment::isMasterRank() )
            std::cout << s << " is not an entry in the tree\n";
    }
    return res;
}

void ModelProperties::saveMD( std::ostream& os )
{
    os << "## Model \n";
    os << " - Name **" << name() << "**\n";
    os << " - shortName **" << shortName() << "**\n";
    os << " - description **" << description() << "**\n";
    os << " - model **" << model() << "**\n\n";
    M_params.saveMD( os );
    M_mat.saveMD( os );
    M_bc.saveMD( os );
    M_postproc.saveMD( os );
}

void ModelProperties::put( std::string const& key, std::string const& entry )
{
    M_p.put( key, entry );
}

void ModelProperties::write( std::string const& f )
{
    pt::write_json( f, M_p );
}
}
