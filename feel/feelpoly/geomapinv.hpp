/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-06-20

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file geomapinv.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-06-20
 */
#ifndef __GeoMapInverse_H
#define __GeoMapInverse_H 1

#include <feel/feelmesh/kdtree.hpp>
#include <feel/feelpoly/geomap.hpp>

namespace Feel
{
/**
 * \class GeoMapInverse
 * \brief handles the geometric inversion for a given (supposedly quite large)
     set of points
 *
 *  @author Christophe Prud'homme
 */
template<int Dim,
         int Order = 1,
         int RealDim = Dim,
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Entity = Simplex>
class GeoMapInverse
{
public:


    /** @name Typedefs
     */
    //@{

    typedef GeoMap<Dim,Order,RealDim,T,Entity> gm_type;
    typedef typename gm_type::Inverse gic_type;
    typedef T value_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    GeoMapInverse( value_type eps = 1e-12 )
        :
        M_eps( eps ),
        M_tree()
    {}
    GeoMapInverse( GeoMapInverse const & gmi )
        :
        M_eps( gmi.M_eps ),
        M_tree( gmi.M_tree )
    {}
    virtual ~GeoMapInverse()
    {}

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
     * clear the kd-tree
     */
    void clear( void )
    {
        M_tree.clear();
    }

    /**
     *  Add the points contained in c to the list of points.
     */
    template<class CONT>
    void
    addPoints( const CONT &c )
    {
        M_tree.reserve( std::distance( c.begin(),c.end() ) );
        typename CONT::const_iterator it = c.begin(), ite = c.end();

        for ( ; it != ite; ++it )
            M_tree.addPointWithId( it->node(), it->id(), 0 );
    }

    /**
     * @return Number of points.
     */
    size_type nPoints( void ) const
    {
        return M_tree.nPoints();
    }

    /**
     *  Add point p to the list of points.
     */
    size_type addPoint( typename node<value_type>::type const& p )
    {
        return M_tree.addPoint( p );
    }

    /**
     * add new points in the kd-tree
     *
     * if lid is equal to invalid_uint16_type_value then it means that it is of no use and will be ignore
     */
    void addPointWithId( typename node<value_type>::type const& p, size_type id, uint16_type comp )
    {
        M_tree.addPointWithId( p, id, comp );
    }

    /**
     * add new points in the kd-tree
     */
    void addPointWithId( boost::tuple<typename node<value_type>::type, size_type, uint16_type > const& p )
    {
        M_tree.addPointWithId( boost::get<0>( p ), boost::get<1>( p ), boost::get<2>( p ) );
    }

    /**
     * Find all the points present in the box between min and max.
     */
    size_type
    pointsInBox( KDTree::points_type &ipts,
                 typename node<value_type>::type const& min,
                 typename node<value_type>::type const& max ) const
    {
        M_tree.pointsInBox( ipts, min, max );
        return ipts.size();
    }


    //@}

protected:

    value_type M_eps;
    mutable KDTree M_tree;
    //gic_type M_gic;

};
} // Feel
#endif /* __GeoMapInverse_H */
