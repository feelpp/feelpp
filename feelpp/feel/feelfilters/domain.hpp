/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

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
   \file domain.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#ifndef FEELPP_FILTERS_DOMAIN_H
#define FEELPP_FILTERS_DOMAIN_H 1

#include <feel/feelconfig.h>

#ifdef FEELPP_HAS_GMSH
#include <feel/feelfilters/gmsh.hpp>

namespace Feel {
/**
 * \brief generate a simple geometrical domain from required and optional parameters
 *
 * List of required parameters:
 *  - \param _name name of the file that will ge generated without extension
 *  - \param _shape shape of the domain to be generated (simplex or hypercube)
 * List of optional parameters:
 *  - \param _dim dimension of the domain (default: 2)
 *  - \param _order order of the geometry (default: 1)
 *  - \param _h characteristic size of the mesh (default: 0.1)
 *  - \param _convex type of convex used to mesh the domain (default: simplex) (simplex or hypercube)
 *  - \param _addmidpoint add middle point (default: true )
 *  - \param _xmin minimum x coordinate (default: 0)
 *  - \param _xmax maximum x coordinate (default: 1)
 *  - \param _ymin minimum y coordinate (default: 0)
 *  - \param _ymax maximum y coordinate (default: 1)
 *  - \param _zmin minimum z coordinate (default: 0)
 *  - \param _zmax maximum z coordinate (default: 1)
 *
 * \attention this function uses the Boost.Parameter library that allows to
 * enter the parameter in any order.
 *
 */
template <typename ... Ts>
gmsh_ptrtype domain( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    std::string const& name = args.get(_name );
    std::string const& prefix = args.get_else( _prefix, "" );
    worldcomm_ptr_t worldcomm = args.get_else( _worldcomm, Environment::worldCommPtr() );
    std::string const& shape = args.get_else_invocable( _shape, [&prefix](){ return soption(_prefix=prefix,_name="gmsh.domain.shape"); } );
    double shear = args.get_else_invocable( _shear, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.domain.shear"); } );
    bool recombine = args.get_else_invocable( _recombine, [&prefix](){ return boption(_prefix=prefix,_name="gmsh.domain.recombine"); } );
    int dim = args.template get_else<std::is_integral>( _dim, 3);
    int order = args.template get_else<std::is_integral>( _order, 1);
    auto && geo_parameters = args.get_else( _geo_parameters, Gmsh::gpstr2map("") );
    double h = args.get_else_invocable(_h, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.hsize"); } );
    std::string const& convex = args.get_else_invocable(_convex, [&prefix](){ return soption(_prefix=prefix,_name="gmsh.domain.convex"); } );
    bool addmidpoint = args.get_else_invocable( _addmidpoint, [&prefix](){ return boption(_prefix=prefix,_name="gmsh.domain.addmidpoint"); } );
    bool usenames = args.get_else_invocable( _usenames, [&prefix](){ return boption(_prefix=prefix,_name="gmsh.domain.usenames"); } );
    double xmin = args.get_else_invocable( _xmin, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.domain.xmin"); } );
    double xmax = args.get_else_invocable( _xmax, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.domain.xmax"); } );
    double ymin = args.get_else_invocable( _ymin, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.domain.ymin"); } );
    double ymax = args.get_else_invocable( _ymax, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.domain.ymax"); } );
    double zmin = args.get_else_invocable( _zmin, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.domain.zmin"); } );
    double zmax = args.get_else_invocable( _zmax, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.domain.zmax"); } );
    double nx = args.get_else_invocable( _nx, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.domain.nx"); } );
    double ny = args.get_else_invocable( _ny, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.domain.ny"); } );
    double nz = args.get_else_invocable( _nz, [&prefix](){ return doption(_prefix=prefix,_name="gmsh.domain.ny"); } );
    bool substructuring = args.get_else_invocable( _substructuring, [&prefix](){ return boption(_prefix=prefix,_name="gmsh.domain.substructuring"); } );

    gmsh_ptrtype gmsh_ptr = Gmsh::New( shape, 3, 1, convex, worldcomm );
    gmsh_ptr->setPrefix( name );
    //gmsh_ptr->addGeoParameters( gmsh_ptr->retrieveGeoParameters( gmsh_ptr->description() ) );
    gmsh_ptr->addGeoParameters( geo_parameters );
    gmsh_ptr->setCharacteristicLength( h );
    gmsh_ptr->setAddMidPoint( addmidpoint );
    gmsh_ptr->usePhysicalNames( usenames );
    gmsh_ptr->setShear( shear );
    gmsh_ptr->setRecombine( recombine );
    gmsh_ptr->setX( std::make_pair( xmin, xmax ) );
    gmsh_ptr->setY( std::make_pair( ymin, ymax ) );
    gmsh_ptr->setZ( std::make_pair( zmin, zmax ) );
    gmsh_ptr->setNx( nx );
    gmsh_ptr->setNy( ny );
    gmsh_ptr->setNz( nz );
    gmsh_ptr->setSubStructuring( substructuring );
    return gmsh_ptr;
}

}

#endif

#endif /* FEELPP_DOMAIN_HPP */
