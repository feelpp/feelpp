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
#include <feel/feelmodels/modelproperties.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace Feel {


ModelProperties::ModelProperties( std::string const& filename )
{
    pt::read_json(filename, M_p);
    try {
        M_name = M_p.get<std::string>( "Name" );
        M_model = M_p.get<std::string>( "Model" );
    }
    catch ( pt::ptree_bad_path& e )
    {
        if ( Environment::isMasterRank() )
            std::cout << "Missing Name or Model entry in model properties\n";
            //std::cout << "Missing path " << e.path() << "\n";
    }
    auto par = M_p.get_child_optional("Parameters");
    if ( par )
    {
        LOG(INFO) << "Model with parameters\n";
        M_params.setPTree( *par );
        
    }
    auto bc = M_p.get_child_optional("BoundaryConditions");
    if ( bc )
    {
        LOG(INFO) << "Model with boundary conditions\n";
        M_bc.setPTree( *bc );
        
    }
    else
    {
        LOG(WARNING) << "Model does not have any boundary conditions\n";
    }
    auto mat = M_p.get_child_optional("Materials");
    if ( mat )
    {
        LOG(INFO) << "Model with materials\n";
        M_mat.setPTree( *mat );
        
    }
    else
    {
        LOG(WARNING) << "Model does not have any materials\n";
    }
    auto pp = M_p.get_child_optional("PostProcess");
    if ( pp )
    {
        LOG(INFO) << "Model with PostProcess\n";
        M_postproc.setPTree( *pp );
        
    }
    else
    {
        LOG(WARNING) << "Model does not have any materials\n";
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

}
