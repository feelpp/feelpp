/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-06-30

  Copyright (C) 2014 Feel++ Consortium

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
/**
   \file ddframework.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-06-30
 */
#include <feel/feel.hpp>

int main(int argc, char const *argv[])
{
  using namespace Feel;
  auto global_mesh=loadMesh(_mesh = new Mesh<Hypercube<2>>, _respect_partition=false );
  // this mesh is for the coarse preconditioner
  auto coarsemesh=loadMesh(_mesh = new Mesh<Hypercube<2>> );

  auto E = global_mesh->element( Environment::processId() );
  auto Emesh = mesh( E );  // TODO using an in-memory gmsh interface to be done
  auto S = Subdomain( E, Emesh ); // TODO

  /// then do something

  return 0;
}
