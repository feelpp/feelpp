/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013 Feel++ Consortium

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
   \file unitcircle.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_UNITCIRCLE_HPP)
#define FEELPP_UNITCIRCLE_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>

namespace Feel {

/**
 * build a mesh of the unit circle using triangles
 */
template<int Ngeo=1>
inline
boost::shared_ptr<Mesh<Simplex<2,Ngeo> > >
unitCircle( double h = doption(_name="gmsh.hsize") )
{
    return createGMSHMesh(_mesh=new Mesh<Simplex<2,Ngeo> >,
                          _desc=domain( _name="circle",
                                        _shape="ellipsoid",
                                        _dim=2,
                                        _xmin=-1,
                                        _ymin=-1,
                                        _h=h ) );
}

}

#endif /* FEELPP_UNITCIRCLE_HPP */
