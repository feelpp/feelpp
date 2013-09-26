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
   \file sphere.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-27
 */
#ifndef __Sphere_H
#define __Sphere_H 1

#include <feel/feelmesh/surface.hpp>

namespace Feel
{
/**
 * This class defines a sphere.  It also computes coordinate
 * transformations between cartesian  \f$ (x, y, z) \f$
 * and spherical  \f$ (r, \theta, \phi) \f$ coordinates.
 * The spherical coordinates are valid in the ranges:
 *
 * - \f$ 0 \le r      < \infty \f$
 * - \f$ 0 \le \theta < \pi \f$
 * - \f$ 0 \le \phi   < 2\pi \f$
 *
 * The coordinates are related as follows:
 * \f$ \phi \f$ is the angle in the xy plane
 * starting with 0. from the positive x axis,
 * \f$ \theta \f$ is measured against the positive
 * z axis.
 \verbatim

 \      | Z
 \theta|
 \    |    .
 \   |   .
 \  |  .
 \ | .
 \|.
 ---------------+---------.---------
 /|\       .          Y
 /phi\     .
 /  |  \   .
 /   |   \ .
 /.........\
 /     |
 X /
 \endverbatim
 *
 * \author Benjamin S. Kirk, Daniel Dreyer, 2002-2003
 * \author Christophe Prud'homme, 2005
 * \date 2002-2003
 */
class Sphere : public Surface
{
    typedef Surface super;

public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    Sphere()
        :
        super(),
        M_radius( -1 )
    {}

    Sphere( Sphere const & s )
        :
        super(),
        M_center( s.M_center ),
        M_radius( s.M_radius )
    {
    }

    Sphere ( const Point& c,
             const double   r )
        :
        super(),
        M_center( c ),
        M_radius( r )
    {
        assert ( r > 0. );

    }



    ~Sphere()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * Returns the radius of the sphere.
     */
    double radius() const
    {
        return M_radius;
    }


    /**
     * @returns the center of the sphere.
     */
    const Point& center() const
    {
        return M_center;
    }




    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the center.
     */
    void setCenter( Point const& p )
    {
        M_center = p;
    }

    /**
     * set the radius
     */
    void setRadius( double r )
    {
        M_radius = r ;
    }


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Defines a sphere of radius r centered at c.
     */
    void createFromCenterRadius ( const Point& c, const double r )
    {
        FEELPP_ASSERT( r > 0 )( r ).error( "radius negative" );
        this->setCenter( c );
        this->setRadius( r );

    }

    /**
     * @returns true if other_sphere intersects this sphere,
     * false otherwise.
     */
    bool intersects ( const Sphere& other_sphere ) const
    {
        FEELPP_ASSERT( M_radius > 0 )( M_radius ).error( "radius negative" );
        FEELPP_ASSERT( other_sphere.radius() > 0 )( other_sphere.radius() ).error( "radius negative" );

        if ( Feel::distance( this->center(), other_sphere.center() ) < ( this->radius() + other_sphere.radius() ) )
            return true;

        return false;
    }


    /**
     * @returns true if the point p is above the surface,
     * false otherwise.
     */
    bool aboveSurface ( const Point& p ) const
    {
        FEELPP_ASSERT( M_radius > 0 )( M_radius ).error( "radius negative" );

        if ( Feel::distance( p, this->center() ) > this->radius() )
            return true;

        return false;
    }


    /**
     * @returns true if the point p is below the surface,
     * false otherwise.
     */
    bool belowSurface ( const Point& p ) const
    {
        return ( !this->aboveSurface ( p ) );
    }


    /**
     * @returns true if the point p is on the surface,
     * false otherwise.  Note that the definition of on
     * the surface really means "very close" to account
     * for roundoff error.
     */
    bool onSurface ( const Point& p ) const
    {
        FEELPP_ASSERT( M_radius > 0 )( M_radius ).error( "radius negative" );


        if ( std::abs( Feel::distance( p, this->center() ) - this->radius() ) < 1.e-10 )
            return true;

        return false;
    }


    /**
     * @return the closest point on the surface to point p.
     */
    Point closestPoint ( const Point& p ) const
    {
        FEELPP_ASSERT( M_radius > 0 )( M_radius ).error( "radius negative" );

        // get the normal from the surface in the direction
        // of p
        Point normal = this->unitNormal ( p );

        // The closest point on the sphere is in the direction
        // of the normal a distance r from the center.
        const Point cp = this->center().node() + normal.node()*this->radius();

        return cp;
    }


    /**
     * @return a unit vector normal to the surface at
     * point p.
     */
    Point unitNormal ( const Point& p ) const
    {
        FEELPP_ASSERT( M_radius > 0 )( M_radius ).error( "radius negative" );

        //assert ( !(p == this->center()) );

        // Create a vector from the center to the point
        Point n = p.node() - this->center().node();
        Point unit_n( n.node()/ ublas::norm_2( n.node() ) );

        return unit_n;
    }



    /**
     * @returns the spherical coordinates for the
     * cartesian coordinates \p cart.
     */
    Point surfaceCoords ( const Point& cart ) const
    {
        // constant translation in the origin
        const Point c ( cart.node() - this->center().node() );

        // phi: special care, so that it gives 0..2pi results
        const double phi = std::atan2( c( 1 ), c( 0 ) );

        return Point( /* radius */ ublas::norm_2( c.node() ),
                                   /* theta  */ std::atan2( std::sqrt( c( 0 )*c( 0 ) + c( 1 )*c( 1 ) ), c( 2 ) ),
                                   /* phi    */ ( ( phi < 0 )  ?  2.*M_PI+phi  :  phi ) );
    }


    /**
     * @returns the cartesian coordinates for the
     * spherical coordinates \p sph.
     */
    Point worldCoords ( const Point& sph ) const
    {
        const double r     = sph( 0 );
        const double theta = sph( 1 );
        const double phi   = sph( 2 );

        // constant translation out of the origin
        return Point ( /* x */ r*std::sin( theta )*std::cos( phi ) + this->center()( 0 ),
                               /* y */ r*std::sin( theta )*std::sin( phi ) + this->center()( 1 ),
                               /* z */ r*std::cos( theta )               + this->center()( 2 ) );
    }




    //@}




private:


    /**
     * The center of the sphere.
     */
    Point M_center;

    /**
     * The radius of the sphere.
     */
    Real  M_radius;
};
} // Feel
#endif /* __Sphere_H */
