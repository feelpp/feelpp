/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-06-14

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file bench1_run2d.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-06-14
 */
#include <bench1_impl.hpp>

namespace Feel
{
void
Bench1::run2d()
{
    const int Dim = 2;
    typedef Mesh<Simplex<Dim> > mesh_type;
    boost::shared_ptr<mesh_type> aMesh;

    std::string shape = vm()["shape"].as<std::string>();


    aMesh = createGMSHMesh( _mesh=new mesh_type,
                            _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                          _usenames=true,
                                          _shape=shape,
                                          _dim=Dim,
                                          _h=meshSize ),
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );

    LOG(INFO) << "run2d starts" << "\n";
    bench1<mesh_type, 1>( aMesh );
    bench1<mesh_type, 2>( aMesh );
    bench1<mesh_type, 4>( aMesh );
    bench1<mesh_type, 6>( aMesh );
    bench1<mesh_type, 8>( aMesh );
    bench1<mesh_type,10>( aMesh );
    LOG(INFO) << "run2d ends" << "\n";

}
}
