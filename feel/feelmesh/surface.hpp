/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-27

  Copyright (C) 2005,2006 EPFL

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
   \file surface.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-27
 */
#ifndef __Surface_H
#define __Surface_H 1

#include <feel/feelmesh/geo0d.hpp>

namespace Feel
{

/**
 * This class defines a surface.  A surface is a two-dimensional
 * object living in three-dimensional space.  Examples of surfaces
 * are planes, hollow spheres, hollow cylinders, etc...  This is
 * a generic base class that describes the useful functionality
 * a surface will provide.  Specific derived classes actually implement
 * the functionality, so this class has pure virtual members.
 *
 * @author Benjamin S. Kirk, 2002
 * @author Christophe Prud'homme, 2005
 */
class Surface
{
public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Constructor.  Does nothing at the moment.
     */
    Surface () {}

    /**
     * Copy-constructor.
     */
    Surface ( const Surface& ) {}

    /**
     * Destructor.
     */
    virtual ~Surface () {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * @returns true if the point p is above the surface,
     * false otherwise.
     */
    virtual bool aboveSurface ( const Point& p ) const = 0;

    /**
     * @returns true if the point p is below the surface,
     * false otherwise.
     */
    virtual bool belowSurface ( const Point& p ) const = 0;

    /**
     * @returns true if the point p is on the surface,
     * false otherwise.  Note that the definition of on
     * the surface really means "very close" to account
     * for roundoff error.
     */
    virtual bool onSurface ( const Point& p ) const = 0;

    /**
     * @returns the closest point on the surface to point p.
     */
    virtual Point closestPoint ( const Point& p ) const = 0;

    /**
     * @returns a unit vector normal to the surface at
     * point p.
     */
    virtual Point unitNormal ( const Point& p ) const = 0;

    /**
     * @returns the \p Point \p world_coords in the
     * surface's coordinate system.  \p world_coords
     * is in the world coordinate system.  This method
     * is not purely virtual, because there may be surfaces
     * that do not have an own coordinate system.  These
     * simply do not have to overload this method.
     */
    virtual Point surfaceCoords ( const Point& world_coords ) const
    {
        Point p ( world_coords );
        return p;
    }



    /**
     * @returns the world (cartesian) coordinates for the
     * surface coordinates \p surf_coords.  This method
     * is not purely virtual, because there may be surfaces
     * that do not have an own coordinate system.  These
     * simply do not have to overload this method.
     */
    virtual Point worldCoords ( const Point& surf_coords ) const
    {
        Point p ( surf_coords );
        return p;
    }

    //@}



protected:

private:

};
} // Feel
#endif /* __Surface_H */
