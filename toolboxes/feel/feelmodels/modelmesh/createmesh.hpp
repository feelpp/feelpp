/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-08-11

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file createmesh.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-08-11
 */

#ifndef FEELPP_MODELS_CREATEMESH_H
#define FEELPP_MODELS_CREATEMESH_H 1

#include <feel/feelmodels/modelcore/modelnumerical.hpp>

namespace Feel
{
namespace FeelModels
{
    template <typename MeshType>
    std::shared_ptr<MeshType>
    reloadMesh(std::string const& nameFile, WorldComm const& worldComm, int straighten=1 );

    template <typename MeshType>
    void
    createMeshModel( ModelNumerical & model, std::shared_ptr<MeshType> & mesh, std::string const& modelMeshRestartFile );

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_CREATEMESH_H
