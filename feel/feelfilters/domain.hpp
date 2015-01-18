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
   \file domain.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_DOMAIN_HPP)
#define FEELPP_DOMAIN_HPP 1

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
BOOST_PARAMETER_FUNCTION(
    ( gmsh_ptrtype ), // return type
    domain,    // 2. function name
    tag,           // 3. namespace of tag types
    ( required
      ( name,           *( boost::is_convertible<mpl::_,std::string> ) )
      )
    ( optional
      ( prefix,(std::string), "" )
      ( worldcomm,      *, Environment::worldComm() )
      ( shape,          *( boost::is_convertible<mpl::_,std::string> ),  soption(_prefix=prefix,_name="gmsh.domain.shape") )
      ( shear,          *( boost::is_arithmetic<mpl::_> )    ,  doption(_prefix=prefix,_name="gmsh.domain.shear") )
      ( recombine,      *( boost::is_integral<mpl::_> )    , boption(_prefix=prefix,_name="gmsh.domain.recombine") )
      ( dim,              *( boost::is_integral<mpl::_> ), 3 )
      ( order,              *( boost::is_integral<mpl::_> ), 1 )
      ( geo_parameters,  *( boost::icl::is_map<mpl::_> ), Gmsh::gpstr2map("") )
      ( h,              *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.hsize") )
      ( convex,         *( boost::is_convertible<mpl::_,std::string> ), soption(_prefix=prefix,_name="gmsh.domain.convex") )
      ( addmidpoint,    *( boost::is_integral<mpl::_> ), boption(_prefix=prefix,_name="gmsh.domain.addmidpoint") )
      ( usenames,       *( boost::is_integral<mpl::_> ), boption(_prefix=prefix,_name="gmsh.domain.usenames") )
      ( xmin,           *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.domain.xmin") )
      ( xmax,           *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.domain.xmax"))
      ( ymin,           *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.domain.ymin") )
      ( ymax,           *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.domain.ymax") )
      ( zmin,           *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.domain.zmin") )
      ( zmax,           *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.domain.zmax") )
      ( nx,             *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.domain.nx") )
      ( ny,             *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.domain.ny") )
      ( nz,             *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.domain.ny") )
      ( substructuring, *( boost::is_integral<mpl::_> ), boption(_prefix=prefix,_name="gmsh.domain.substructuring") ) ) )
{
    gmsh_ptrtype gmsh_ptr = Gmsh::New( shape, 3, 1, convex, worldcomm );
    gmsh_ptr->setPrefix( name );
    gmsh_ptr->setGeoParameters( gmsh_ptr->retrieveGeoParameters( gmsh_ptr->description() ), 0 );
    gmsh_ptr->setGeoParameters( geo_parameters );
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

#endif /* FEELPP_DOMAIN_HPP */
