/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-15

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
   \file boundingbox.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-15
 */
#ifndef __BoundingBox_H
#define __BoundingBox_H 1

#include <feel/feelalg/glas.hpp>

namespace Feel
{
/*!
  \class BoundingBox
  \brief bounding box for a matrix of points

  a matrix of points is a matrix of size (dimension x number of
  points).

  the rows contains the coordinates while the columns are the points

  @author Christophe Prud'homme
  @see
*/
template<typename T = double>
struct BoundingBox
{
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef typename node<T>::type node_type;
    typedef typename matrix_node<T>::type matrix_node_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    BoundingBox( bool is_lin = true )
        :
        is_linear( is_lin )
    {}
    BoundingBox( BoundingBox const & bb )
        :
        is_linear( bb.is_linear ),
        min( bb.min ),
        max( bb.max )
    {}
    ~BoundingBox()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{
    void
    make( matrix_node_type const& __ptab )
    {
        typename matrix_node_type::const_iterator2 it = __ptab.begin2();
        typename matrix_node_type::const_iterator2 en = __ptab.end2();

        size_type P = __ptab.size1();

        min.resize( __ptab.size1() );
        max.resize( __ptab.size1() );

        min = ublas::column( __ptab, 0 );
        max = ublas::column( __ptab, 0 );

        ++it;

        typename node_type::iterator itmin = min.begin();
        typename node_type::iterator itmax = max.begin();

        while ( it != en )
        {
            typename matrix_node_type::const_iterator1 it1 = it.begin();

            for ( size_type i = 0; i < P; ++i, ++it1 )
            {
                min[i] = std::min( min[i], *it1 );
                max[i] = std::max( max[i], *it1 );
            }

#if 0
            std::for_each( it.begin(), it.end(),
                           ( lambda::var( min[i] ) = std::min( lambda::var( min[i] ), lambda::_1 ),
                             lambda::var( max[i] ) = std::max( lambda::var( max[i] ), lambda::_1 ),
                             std::cout << "min: " << lambda::var( min ) << "\n",
                             std::cout << "max: " << lambda::var( max ) << "\n" ) );
#endif

            // enlarge the box for non-linear transformations
            if ( !is_linear )
            {
                for ( size_type i = 0; i < P; ++i )
                {
                    value_type e = ( max[i]-min[i] ) * 0.2;
                    min[i] -= e;
                    max[i] += e;
                }
            }

            ++it;
        }

    }


    //@}

    /** @name Accessors
     */
    //@{

    bool isLinear() const
    {
        return is_linear;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}


    node_type min;
    node_type max;
    bool is_linear;

};
}
#endif /* __BoundingBox_H */
