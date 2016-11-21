/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013-2016 Feel++ Consortium

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
   \file unitsquare.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#ifdef FEELPP_HAS_GMSH
#include <feel/feelfilters/creategmshmesh.hpp>
#endif
#include <feel/feelfilters/domain.hpp>

#ifdef FEELPP_HAS_GMSH
namespace Feel {

boost::shared_ptr<Mesh<Simplex<2> > >
unitSquare( double h, PeriodicEntities pe )
{
    return createGMSHMesh(_mesh=new Mesh<Simplex<2> >,
                          _desc=domain( _name="square",
                                        _shape="hypercube",
                                        _dim=2,
                                        _h=h),
                          _periodic = pe );
}

}
#endif

