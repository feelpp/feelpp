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
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>

#include <feel/feelmodels/modelproperties.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace Feel {


nl::json
ModelProperties::read_json( std::string const& _filename, worldcomm_ptr_t const& world )
{
    std::string filename = Environment::expand( _filename );
    if ( !fs::exists( filename ) )
    {
        if ( world->isMasterRank() )
        {
            std::cout << "[modelProperties]  Could not find :\"" << filename << "\"" <<std::endl;
        }
        LOG(INFO) << "Could not find " << filename << std::endl;
        return {};
    }
    else
    {
        if ( world->isMasterRank() )
        {
            std::cout << "[modelProperties] Loading Model Properties : \"" << filename << "\"" << std::endl;
        }
        LOG(INFO) << "Loading " << filename << std::endl;
    }
#if 0
    auto json_str_wo_comments = removeComments(readFromFile(filename));
    LOG(INFO) << "json file without comment:" << json_str_wo_comments;
    std::istringstream istr( json_str_wo_comments );
    pt::read_json(istr, M_p);
#endif
    std::ifstream ifs(filename);
    return nl::json::parse(ifs,nullptr,true,true);
}

ModelProperties::ModelProperties( std::string const& directoryLibExpr, worldcomm_ptr_t const& world, std::string const& prefix, po::variables_map const& vm )
    :
    super( world ),
    M_prefix( prefix ),
    M_directoryLibExpr( directoryLibExpr ),
    M_params( world ),
    M_mat( world ),
    M_ic( world ),
    M_bc2( world ),
    M_postproc( world ),
    M_outputs( world )
{
    if ( countoption( _name="json.merge_patch",_prefix=M_prefix,_vm=vm ) > 0 )
    {
        for ( std::string const& patch : vsoption( _name="json.merge_patch",_prefix=M_prefix,_vm=vm ) )
            M_json_merge_patch.merge_patch( json::parse( patch ) );
    }
    if ( countoption( _name="json.patch",_prefix=M_prefix,_vm=vm ) > 0 )
    {
        for ( std::string const& patch : vsoption( _name="json.patch", _prefix=M_prefix,_vm=vm ) )
            M_json_patch.push_back( json::parse( patch ) );
    }
}

ModelProperties::ModelProperties( std::string const& filename, std::string const& directoryLibExpr, worldcomm_ptr_t const& world, std::string const& prefix, po::variables_map const& vm )
    :
    ModelProperties( ModelProperties::read_json( filename, world ),directoryLibExpr,world,prefix,vm )
{}

void
ModelProperties::setup( std::vector<std::string> const& filenames )
{
    nl::json j;
    for ( std::string const& filename : filenames )
        j.merge_patch( ModelProperties::read_json( filename, this->worldCommPtr() ) );
    this->setup( std::move( j ) );
}

void
ModelProperties::setupFromFilenameOption( po::variables_map const& vm )
{
    if ( countoption( _name="json.filename",_prefix=M_prefix,_vm=vm ) > 0 )
        this->setup( vsoption( _name="json.filename",_prefix=M_prefix,_vm=vm ) );
}

void
ModelProperties::setupImpl()
{
    if ( !M_json_merge_patch.is_null() )
        M_jsonData.merge_patch( M_json_merge_patch );
    if ( !M_json_patch.empty() )
        M_jsonData = M_jsonData.patch( M_json_patch );
    nl::json & jarg = M_jsonData;
    //std::cout << jarg.dump(1) << std::endl;
    if ( jarg.contains("Name") )
    {
        auto const& j_name = jarg.at("Name");
        if ( j_name.is_string() )
            M_name = j_name.get<std::string>();
    }
    if ( jarg.contains("ShortName") )
    {
        auto const& j_shortname = jarg.at("ShortName");
        if ( j_shortname.is_string() )
            M_shortname = j_shortname.get<std::string>();
    }

    M_unit = "m";
    if ( jarg.contains("Unit") )
    {
        auto const& j_unit = jarg.at("Unit");
        if ( j_unit.is_string() )
            M_unit = j_unit.get<std::string>();
    }

    if ( jarg.contains("Models") )
    {
        LOG(INFO) << "Model with model\n";
        M_models.setPTree( jarg.at("Models") );
    }

    if ( jarg.contains("Parameters") )
    {
        LOG(INFO) << "Model with parameters\n";
        if ( !M_directoryLibExpr.empty() )
            M_params.setDirectoryLibExpr( M_directoryLibExpr );
        M_params.setPTree( jarg.at("Parameters") );
    }
    if ( jarg.contains("BoundaryConditions") )
    {
        LOG(INFO) << "Model with boundary conditions\n";
        if ( M_bc2_enabled )
        {
            pt::ptree pt;
            std::istringstream istr( jarg.at("BoundaryConditions").dump() );
            pt::read_json(istr, pt);
            M_bc2.setPTree( pt );
        }
        else
            M_bc.setup( jarg.at("BoundaryConditions") );
    }

    if ( jarg.contains("InitialConditions") )
    {
        LOG(INFO) << "Model with initial conditions\n";
        if ( !M_directoryLibExpr.empty() )
            M_ic.setDirectoryLibExpr( M_directoryLibExpr );
        M_ic.setPTree( jarg.at("InitialConditions") );
    }

    if ( jarg.contains("Materials") )
    {
        LOG(INFO) << "Model with materials\n";
        if ( !M_directoryLibExpr.empty() )
            M_mat.setDirectoryLibExpr( M_directoryLibExpr );
        M_mat.setPTree( jarg.at("Materials") );
    }

    if ( jarg.contains("PostProcess") )
    {
        LOG(INFO) << "Model with PostProcess\n";
        if ( !M_directoryLibExpr.empty() )
            M_postproc.setDirectoryLibExpr( M_directoryLibExpr );
        M_postproc.setPTree( jarg.at("PostProcess") );
    }

    if ( jarg.contains("Outputs") )
    {
        LOG(INFO) << "Model with outputs\n";
        if ( !M_directoryLibExpr.empty() )
            M_outputs.setDirectoryLibExpr( M_directoryLibExpr );
        M_outputs.setPTree( jarg.at("Outputs") );
    }
}

}
