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



ModelProperties::ModelProperties( std::string const& filename, std::string const& directoryLibExpr, worldcomm_ptr_t const& world, std::string const& prefix )
    :
    super( world ),
    M_prefix( prefix ),
    M_directoryLibExpr( directoryLibExpr ),
    M_params( world ),
    M_mat( world ),
    M_bc( world, false ),
    M_ic( world ),
    M_bc2( world ),
    M_postproc( world ),
    M_outputs( world )
{
    if ( !fs::exists( filename ) )
    {
        if ( Environment::isMasterRank() )
        {
            std::cout << "[modelProperties]  Could not find :\"" << filename << "\"" <<std::endl;
        }
        LOG(INFO) << "Could not find " << filename << std::endl;
        return;
    }
    else
    {
        if ( Environment::isMasterRank() )
        {
            std::cout << "[modelProperties] Loading Model Properties : \"" << filename << "\"" << std::endl;
        }
        LOG(INFO) << "Loading " << filename << std::endl;
    }


    auto json_str_wo_comments = removeComments(readFromFile(filename));
    LOG(INFO) << "json file without comment:" << json_str_wo_comments;


    std::istringstream istr( json_str_wo_comments );
    pt::read_json(istr, M_p);

    this->setup();
}
ModelProperties::ModelProperties( pt::ptree const& pt, std::string const& directoryLibExpr, worldcomm_ptr_t const& world, std::string const& prefix )
    :
    super( world ),
    M_prefix( prefix ),
    M_directoryLibExpr( directoryLibExpr ),
    M_p( pt )
{
    this->setup();
}
ModelProperties::ModelProperties( nl::json const& j, std::string const& directoryLibExpr, worldcomm_ptr_t const& world, std::string const& prefix )
    :
    super( world ),
    M_prefix( prefix ),
    M_directoryLibExpr( directoryLibExpr ),
    M_p()
{
    auto json_str_wo_comments = removeComments(j.dump(1));
    LOG(INFO) << "json file without comment:" << json_str_wo_comments;
    std::istringstream istr( json_str_wo_comments );
    pt::read_json(istr, M_p);
    this->setup();
}

void
ModelProperties::setup()
{
    editPtreeFromOptions( M_p, M_prefix );

    std::ostringstream pt_ostr;
    write_json( pt_ostr, M_p );
    std::istringstream pt_istream( pt_ostr.str() );
    nl::json jarg;
    pt_istream >> jarg;

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

    auto bc = M_p.get_child_optional("BoundaryConditions");
    if ( bc )
    {
        LOG(INFO) << "Model with boundary conditions\n";
        if ( !M_directoryLibExpr.empty() )
            M_bc.setDirectoryLibExpr( M_directoryLibExpr );
        M_bc.setPTree( *bc );
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

ModelProperties::~ModelProperties()
{}

std::string ModelProperties::getEntry(std::string &s)
{
  std::string res;
  try{
  res = M_p.get<std::string>(s);
  }
  catch( pt::ptree_bad_path& e)
  {
        if ( Environment::isMasterRank() )
            std::cout << s << " is not an entry in the tree\n";
  }
  return res;
}

void ModelProperties::saveMD(std::ostream &os)
{
  os << "## Model \n";
  os << " - Name **" << name()<< "**\n";
  os << " - shortName **" << shortName()<< "**\n";
  os << " - description **" << description()<< "**\n";
  M_params.saveMD(os);
  M_mat.saveMD(os);
  M_bc.saveMD(os);
  M_postproc.saveMD(os);
}

void ModelProperties::put(std::string const &key, std::string const &entry)
{
  M_p.put(key,entry);
}

void ModelProperties::write(std::string const &f)
{
    pt::write_json(f,M_p);
}

ModelProperties&
ModelProperties::enableBoundaryConditions2()
{
    M_bc2_enabled = true; 
    auto bc = M_p.get_child_optional("BoundaryConditions");
    if ( bc )
    {
        VLOG(1) << "Model with boundary conditions\n";
        M_bc2.setPTree( *bc );
    }
    return *this; 
}
}
