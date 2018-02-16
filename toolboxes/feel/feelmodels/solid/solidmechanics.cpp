/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2014-06-04

  Copyright (C) 2014 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file solidmechanics.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-06-04
 */

#include <feel/feelmodels/solid/solidmechanics.hpp>

#include <feel/feelmodels/modelvf/solidmecgeomapeulerian.hpp>

namespace Feel
{
namespace FeelModels
{

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::SolidMechanics( std::string const& prefix,
                                                    bool buildMesh,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    ModelBaseRepository const& modelRep )
    :
    super_type( prefix, buildMesh, worldComm, subPrefix, modelRep )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());

    //-----------------------------------------------------------------------------//
    // load info from .json file
    this->loadConfigBCFile();
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // set of worldComm for the function spaces
    this->createWorldsComm();
    //-----------------------------------------------------------------------------//
    // build  mesh, space,exporter,...
    if (buildMesh) this->build();
    //-----------------------------------------------------------------------------//

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::self_ptrtype
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::New( std::string const& prefix,
                                         bool buildMesh,
                                         WorldComm const& worldComm,
                                         std::string const& subPrefix,
                                         ModelBaseRepository const& modelRep )
{
    return boost::make_shared<self_type>( prefix, buildMesh, worldComm, subPrefix, modelRep );
}








} // namespace FeelModels
} // namespace Feel
