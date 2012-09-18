/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-03-03

  Copyright (C) 2006 EPFL

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
   \file electrostatic.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2008-05-05
 */

#ifndef __electrostatic_H
#define __electrostatic_H 1

#include <feel/feelpoly/warpblend.hpp>

namespace Feel
{
template< class Convex,
          uint16_type Order,
          typename T = double >
class PointSetElectrostatic : public  PointSetInterpolation<Convex::nDim, Order, T, Simplex>
{

public :

    typedef PointSetWarpBlend<Convex, Order, T> super;

    typedef T value_type;

    static const uint32_type Dim = Convex::nDim;
    static const uint16_type nPoints1D = Order+1;

    typedef typename super::return_type return_type;

    typedef ublas::vector<value_type> vector_type;

    static const uint32_type topological_dimension = Convex::topological_dimension;
    static const uint32_type numPoints = super::numPoints;

    typedef Reference<Convex, Dim, Convex::nOrder, Dim, value_type> reference_convex_type;

    typedef typename reference_convex_type::points_type points_type;

    static const uint32_type nbPtsPerFace = super::nbPtsPerFace;

    typedef std::vector<uint16_type> orbits_type;

    reference_convex_type RefConv;

    PointSetElectrostatic( int interior = 0 )
    {
        PointSetWarpBlend<Convex, Order, T> G( interior );

        points_type final_pts = G.points();

        //Copies information about WarpBlend points
        this->setEid( G.getEid() );
        this->setPtE( G.getPtE() );

        if ( Order > 2 )
        {
            entities = std::make_pair( RefConv.vertices(), equiVertices() );

            defineOrbits();

            points_type Gt( Dim+1, nbPtsPerFace );

            uint16_type p=0;

            for ( uint16_type i = 0; i < getPrimalNumber(); i++ )
            {
                points_type orbit = generateOrbit( getOrbitId( i ),
                                                   ublas::subrange( getPrimalPts(), 0, 3, i, i+1 ) );

                ublas::subrange( Gt, 0, 3, p, p+orbit.size2() ) = orbit;

                p += orbit.size2();
            }

            points_type pts = toEquilateral( toCartesian( Gt ) );

            pts = orderInteriorFace( pts );

            if ( interior )
                final_pts = pts;

            else
                final_pts = putInPointset ( final_pts, pts, G.interiorRangeById( 2, 0 ) );

        }

        this->setPoints( final_pts );

        this->setName( "fekete", Order );
    }

    ~PointSetElectrostatic() {}

private :

    std::pair<points_type, points_type> entities;

    orbits_type orbits;

    void addOrbit( uint16_type orbit_id, uint16_type num )
    {
        for ( uint16_type i=0; i<num; i++ )
            orbits.push_back( orbit_id );
    }

    //only define orbits for interior points; all the other we know that are the gausslobatto points in the edges
    void defineOrbits()
    {
        if ( Order == 3 )
            addOrbit( 1, 1 );

        else if ( Order == 4 )
            addOrbit( 3, 1 );

        else if ( Order == 5 )
            addOrbit( 3, 2 );

        else if ( Order == 6 )
        {
            addOrbit( 1, 1 );
            addOrbit( 3, 1 );
            addOrbit( 6, 1 );
        }

        else if ( Order == 7 )
        {
            addOrbit( 3, 3 );
            addOrbit( 6, 1 );
        }

        else if ( Order == 8 )
        {
            addOrbit( 3, 3 );
            addOrbit( 6, 2 );
        }

        else if ( Order == 9 )
        {
            addOrbit( 1, 1 );
            addOrbit( 3, 3 );
            addOrbit( 6, 3 );
        }

        else if ( Order == 10 )
        {
            addOrbit( 3, 4 );
            addOrbit( 6, 4 );
        }

        else if ( Order == 11 )
        {
            addOrbit( 3, 5 );
            addOrbit( 6, 5 );
        }

        else if ( Order == 12 )
        {
            addOrbit( 1, 1 );
            addOrbit( 3, 4 );
            addOrbit( 6, 7 );
        }

        else if ( Order == 13 )
        {
            addOrbit( 3, 6 );
            addOrbit( 6, 8 );
        }
    }

    uint16_type getOrbitId( uint16_type i )
    {
        return orbits[i];
    }

    // returns the number of points defined in the triangle with which we calculate the orbits
    uint16_type getPrimalNumber()
    {
        uint16_type __num[ 16 ] = { 1,
                                    1,
                                    1,
                                    1,
                                    2,
                                    3,
                                    4,
                                    5
                                  };

        return __num[Order-1];
    }

    // defines the points with which we calculate the orbits
    points_type getPrimalPts()
    {
        points_type bar_coord( Dim+1, getPrimalNumber() );

        if ( Order == 3 )
        {
            bar_coord( 0,0 ) = value_type( 1.0/3.0 );
            bar_coord( 1,0 ) = value_type( 1.0/3.0 );
        }

        else if ( Order == 4 )
        {
            bar_coord( 0,0 ) = value_type( 0.2371200168 );
            bar_coord( 1,0 ) = value_type( 0.2371200168 );
        }

        else if ( Order == 5 )
        {
            bar_coord( 0,0 ) = value_type( 0.410515151 );
            bar_coord( 1,0 ) = value_type( 0.410515151 );

            bar_coord( 0,1 ) = value_type( 0.1575181512 );
            bar_coord( 1,1 ) = value_type( 0.1575181512 );
        }

        else if ( Order == 6 )
        {
            bar_coord( 0,0 ) = value_type( 1.0/3 );
            bar_coord( 1,0 ) = value_type( 1.0/3 );

            bar_coord( 0,1 ) = value_type( 0.1061169285 );
            bar_coord( 1,1 ) = value_type( 0.1061169285 );

            bar_coord( 0,2 ) = value_type( 0.3097982151 );
            bar_coord( 1,2 ) = value_type( 0.5569099204 );
        }

        else if ( Order == 7 )
        {
            bar_coord( 0,0 ) = value_type( 0.4477725053 );
            bar_coord( 1,0 ) = value_type( 0.4477725053 );

            bar_coord( 0,1 ) = value_type( 0.2604038024 );
            bar_coord( 1,1 ) = value_type( 0.2604038024 );

            bar_coord( 0,2 ) = value_type( 0.0660520784 );
            bar_coord( 1,2 ) = value_type( 0.0660520784 );

            bar_coord( 0,3 ) = value_type( 0.2325524777 );
            bar_coord( 1,3 ) = value_type( 0.6759625951 );
        }

#if 0

        else if ( Order == 8 )
        {
            bar_coord( 0,0 ) = value_type( 0.8444999609 );
            bar_coord( 1,0 ) = value_type( 0.0777500195 );

            bar_coord( 0,1 ) = value_type( 0.4683305115 );
            bar_coord( 1,1 ) = value_type( 0.4683305115 );

            bar_coord( 0,2 ) = value_type( 0.3853668203 );
            bar_coord( 1,2 ) = value_type( 0.3853668204 );

            bar_coord( 0,3 ) = value_type( 0.7172965409 );
            bar_coord( 1,3 ) = value_type( 0.2335581033 );

            bar_coord( 0,4 ) = value_type( 0.5853134902 );
            bar_coord( 1,4 ) = value_type( 0.2667701010 );
        }

        else if ( Order == 9 )
        {
            bar_coord( 0,0 ) = value_type( 1.0/3 );
            bar_coord( 1,0 ) = value_type( 1.0/3 );

            bar_coord( 0,1 ) = value_type( 0.1704318201 );
            bar_coord( 1,1 ) = value_type( 0.1704318201 );

            bar_coord( 0,2 ) = value_type( 0.0600824712 );
            bar_coord( 1,2 ) = value_type( 0.4699587644 );

            bar_coord( 0,3 ) = value_type( 0.0489345696 );
            bar_coord( 1,3 ) = value_type( 0.0489345696 );

            bar_coord( 0,4 ) = value_type( 0.1784337588 );
            bar_coord( 1,4 ) = value_type( 0.3252434900 );

            bar_coord( 0,5 ) = value_type( 0.0588564879 );
            bar_coord( 1,5 ) = value_type( 0.3010242110 );

            bar_coord( 0,6 ) = value_type( 0.0551758079 );
            bar_coord( 1,6 ) = value_type( 0.1543901944 );
        }

        else if ( Order == 10 )
        {
            bar_coord( 0,0 ) = value_type( 0.9145236987 );
            bar_coord( 1,0 ) = value_type( 0.0427381507 );

            bar_coord( 0,1 ) = value_type( 0.5331019411 );
            bar_coord( 1,1 ) = value_type( 0.2334490294 );

            bar_coord( 0,2 ) = value_type( 0.4814795342 );
            bar_coord( 1,2 ) = value_type( 0.4814795342 );

            bar_coord( 0,3 ) = value_type( 0.3800851251 );
            bar_coord( 1,3 ) = value_type( 0.3800851251 );

            bar_coord( 0,4 ) = value_type( 0.8150971991 );
            bar_coord( 1,4 ) = value_type( 0.1351329831 );

            bar_coord( 0,5 ) = value_type( 0.6778669104 );
            bar_coord( 1,5 ) = value_type( 0.2844305545 );

            bar_coord( 0,6 ) = value_type( 0.6759450113 );
            bar_coord( 1,6 ) = value_type( 0.2079572403 );

            bar_coord( 0,7 ) = value_type( 0.5222323306 );
            bar_coord( 1,7 ) = value_type( 0.3633472465 );
        }

        else if ( Order == 11 )
        {
            bar_coord( 0,0 ) = value_type( 0.9201760661 );
            bar_coord( 1,0 ) = value_type( 0.0399119670 );

            bar_coord( 0,1 ) = value_type( 0.8097416696 );
            bar_coord( 1,1 ) = value_type( 0.0951291652 );

            bar_coord( 0,2 ) = value_type( 0.4216558161 );
            bar_coord( 1,2 ) = value_type( 0.2891720920 );

            bar_coord( 0,3 ) = value_type( 0.4200100315 );
            bar_coord( 1,3 ) = value_type( 0.4200100315 );

            bar_coord( 0,4 ) = value_type( 0.4832770031 );
            bar_coord( 1,4 ) = value_type( 0.4832770031 );

            bar_coord( 0,5 ) = value_type( 0.8236881237 );
            bar_coord( 1,5 ) = value_type( 0.1452587341 );

            bar_coord( 0,6 ) = value_type( 0.7030268141 );
            bar_coord( 1,6 ) = value_type( 0.2021386640 );

            bar_coord( 0,7 ) = value_type( 0.6642752329 );
            bar_coord( 1,7 ) = value_type( 0.3066778199 );

            bar_coord( 0,8 ) = value_type( 0.5605605456 );
            bar_coord( 1,8 ) = value_type( 0.3510551601 );

            bar_coord( 0,9 ) = value_type( 0.5584153138 );
            bar_coord( 1,9 ) = value_type( 0.2661283688 );

        }

        else if ( Order == 12 )
        {
            bar_coord( 0,0 ) = value_type( 1.0/3 );
            bar_coord( 1,0 ) = value_type( 1.0/3 );

            bar_coord( 0,1 ) = value_type( 0.1988883477 );
            bar_coord( 1,1 ) = value_type( 0.4005558262 );

            bar_coord( 0,2 ) = value_type( 0.2618405201 );
            bar_coord( 1,2 ) = value_type( 0.2618405201 );

            bar_coord( 0,3 ) = value_type( 0.0807386775 );
            bar_coord( 1,3 ) = value_type( 0.0807386775 );

            bar_coord( 0,4 ) = value_type( 0.0336975736 );
            bar_coord( 1,4 ) = value_type( 0.0336975736 );

            bar_coord( 0,5 ) = value_type( 0.1089969290 );
            bar_coord( 1,5 ) = value_type( 0.3837518758 );

            bar_coord( 0,6 ) = value_type( 0.1590834479 );
            bar_coord( 1,6 ) = value_type( 0.2454317980 );

            bar_coord( 0,7 ) = value_type( 0.0887037176 );
            bar_coord( 1,7 ) = value_type( 0.1697134458 );

            bar_coord( 0,8 ) = value_type( 0.0302317829 );
            bar_coord( 1,8 ) = value_type( 0.4071849276 );

            bar_coord( 0,9 ) = value_type( 0.0748751152 );
            bar_coord( 1,9 ) = value_type( 0.2874821712 );

            bar_coord( 0,10 ) = value_type( 0.0250122615 );
            bar_coord( 1,10 ) = value_type( 0.2489279690 );

            bar_coord( 0,11 ) = value_type( 0.0262645218 );
            bar_coord( 1,11 ) = value_type( 0.1206826354 );
        }

#endif

        for ( uint16_type i=0; i < getPrimalNumber(); i++ )
            bar_coord( 2,i ) = 1 - bar_coord( 0,i ) - bar_coord( 1,i );

        return bar_coord;
    }

    points_type cyclic( points_type pts, uint16_type permOrder )
    {
        if ( permOrder > 0 )
        {
            for ( uint16_type i=0; i<pts.size2(); i++ )
            {
                value_type a = pts( 0,i );
                value_type b = pts( 1,i );
                value_type c = pts( 2,i );

                pts( 0,i ) = c;
                pts( 1,i ) = a;
                pts( 2,i ) = b;
            }
        }

        if ( permOrder > 1 )
            pts = cyclic( pts, permOrder-1 );

        return pts;
    }


    points_type symmetries( points_type pts, uint16_type permOrder )
    {
        if ( permOrder < 3 )
            pts = cyclic( pts, permOrder );

        else
        {
            for ( uint16_type i=0; i<pts.size2(); i++ )
            {
                value_type a = pts( 0,i );
                value_type b = pts( 1,i );
                value_type c = pts( 2,i );

                pts( 0,i ) = a;
                pts( 1,i ) = c;
                pts( 2,i ) = b;

                pts = cyclic( pts, permOrder-3 );
            }
        }

        return pts;
    }


    points_type generateOrbit( uint16_type numOrbit, points_type orbitPoint )
    {
        points_type orbit( Dim+1, numOrbit );

        if ( numOrbit == 1 )
            orbit = orbitPoint;

        else
        {
            for ( uint16_type i=0; i<numOrbit; i++ )
                ublas::subrange( orbit, 0, 3, i, i+1 ) = symmetries( orbitPoint, i );
        }

        return orbit;
    }


    points_type equiVertices ()
    {
        points_type V ( ublas::scalar_matrix<value_type>( Dim, Dim+1, value_type( 0 ) ) );

        for ( uint16_type i=0; i < 3; i++ )
        {
            value_type angle = M_PI*( value_type( 7 ) + value_type( 4 )*value_type( i ) )/value_type( 6 );

            V( 0,i ) = math::cos( angle );
            V( 1,i ) = math::sin( angle );
        }

        V *= value_type( 2 )/math::sqrt( value_type( 3 ) );

        return V;
    }

    vector_type getVertex ( uint16_type element, uint16_type i )
    {
        vector_type p;

        if ( element == 0 )
            p = ublas::column( entities.first, i );

        else
            p = ublas::column( entities.second, i );

        return p;
    }

    points_type toCartesian ( points_type pts )
    {
        points_type C( Dim, Dim );

        for ( uint16_type i = 0; i < Dim; i++ )
            ublas::column( C, i ) = getVertex( 1, ( Dim - Dim%2 + i )%Dim ) - getVertex( 1,Dim );

        points_type bar_coord = ublas::subrange( pts, 1, Dim+1, 0, pts.size2() );

        points_type Gt = ublas::prod( C, bar_coord );

        for ( uint16_type i=0; i < pts.size2(); i++ )
            ublas::column( Gt, i ) += getVertex( 1,Dim );

        return Gt;
    }

    points_type toEquilateral ( points_type pts )
    {
        points_type coord_ref_elem ( Dim, Dim );
        points_type coord_equi_elem ( Dim, Dim );

        for ( uint16_type i = 0; i < Dim; i++ )
        {
            ublas::column( coord_ref_elem, i ) = getVertex( 0,Dim ) - getVertex( 0,i );
            ublas::column( coord_equi_elem, i ) = getVertex( 1,Dim ) - getVertex( 1,i );
        }

        points_type A ( Dim, Dim );
        points_type b ( Dim, 1 );

        LU< points_type > lu( coord_ref_elem );
        lu.inverse( A );

        A = ublas::prod( coord_equi_elem , A );

        ublas::column( b,0 ) = getVertex( 1,0 ) - ublas::prod( A, getVertex( 0,0 ) );

        for ( uint16_type i=0; i < pts.size2(); i++ )
            ublas::column( pts, i ) -= ublas::column( b,0 );

        LU< points_type > lu2( A );
        pts = lu2.solve( pts ) ;

        return pts;
    }

    points_type putInPointset ( points_type final_pts, points_type pts, std::pair<uint16_type, uint16_type> position )
    {
        ublas::subrange( final_pts, 0, 2, position.first, position.second ) = pts;

        return final_pts;
    }

    template<int n>
    static bool order( vector_type a, vector_type b )
    {
        return a( n ) < b( n );
    }

    points_type orderInteriorFace( points_type pts )
    {
        std::vector<vector_type> aux( pts.size2() );

        for ( uint16_type i=0; i<pts.size2(); i++ )
            aux[i] = ublas::column( pts, i );

        std::sort( aux.begin(), aux.end(), &order<1> );

        for ( uint16_type p=0, i=Order-3; i>0; i-- )
        {
            std::sort( aux.begin()+p, aux.begin()+p+i+1, &order<0> );
            p += i+1;
        }

        for ( uint16_type i=0; i<pts.size2(); i++ )
            ublas::column( pts, i ) = aux[i];

        return pts;
    }
};

} // Feel
#endif /* __Electrostatic_H */
