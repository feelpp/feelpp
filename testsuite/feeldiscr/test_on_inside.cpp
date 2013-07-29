/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2008-02-07

Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
  \file generecrash.cpp
  \author Vincent HUBER <vincent.huber@cemosis.fr>
  \date 2013-08-22
  */
#include <feel/feel.hpp>
using namespace Feel;

  int
main( int argc, char** argv )
{
  // Initialize Feel++ Environment
  Environment env( _argc=argc, _argv=argv,
      _desc=feel_options(),
      _about=about( _name="test_on_inside" ,
        _author="Feel++ Consortium",
        _email="feelpp-devel@feelpp.org" ) );

  // create the mesh (specify the dimension of geometric entity)
  typedef Mesh<Simplex<2>> mesh_type;

  auto mesh = loadMesh(_mesh = new mesh_type );

  double nFaces_glob = nelements( markedfaces(mesh,"toto"), true );
  double length = integrate( _range = markedfaces(mesh,"toto"),_expr = cst(1.0) ).evaluate()( 0,0 );

  LOG(INFO) << "global number of faces : " << nFaces_glob ;
  LOG(INFO) << "length : " << length << " expected value : " << math::sqrt(2*.8*.8) << " error " << std::fabs(length -  math::sqrt(2*.8*.8))/std::fabs(math::sqrt(2*.8*.8));

  CHECK( nFaces_glob > 1 )
    << "Invalid number of faces marked toto: "
    << nFaces_glob ;

  CHECK( math::abs( length - math::sqrt(2*.8*.8) ) < 1e-10 )
    << "wrong line interface length between (.1,.1) amd (.9,.9)"
    << " length : " << length
    << " expected value : " << math::sqrt(2*.8*.8);

}
