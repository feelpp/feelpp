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
   \file warpblend.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2006-10-02
 */
#ifndef __WarpBlend_H
#define __WarpBlend_H 1

#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/pointsetinterpolation.hpp>
#include <feel/feelpoly/jacobi.hpp>
#include <feel/feelpoly/geomap.hpp>
#include <feel/feelmesh/geond.hpp>
#include <feel/feelmesh/geo0d.hpp>
#include <feel/feelpoly/gausslobatto.hpp>

#include <feel/feelpoly/equispaced.hpp>


namespace Feel
{

template< class Convex,
          uint16_type Order,
          typename T = double >
class PointSetWarpBlend : public  PointSetInterpolation<Convex::nDim, Order, T, Simplex>
{

public :

    typedef PointSetEquiSpaced<Convex, Order, T> super;

    typedef T value_type;

    static const uint32_type Dim = Convex::nDim;
    static const uint16_type nPoints1D = Order+1;

    typedef typename super::return_type return_type;

    typedef ublas::vector<value_type> vector_type;

    static const uint32_type topological_dimension = Convex::topological_dimension;
    static const uint32_type numPoints = super::numPoints;

    typedef Reference<Convex, Dim, Convex::nOrder, /*Dim*/Convex::nRealDim, value_type> reference_convex_type;

    typedef typename reference_convex_type::points_type points_type;

    reference_convex_type RefConv;

    PointSetWarpBlend( int interior = 0 )
    {
        PointSetEquiSpaced<Convex, Order, T> G( interior );

        points_type final_pts = G.points();

        //Copies information about EquiSpaced points
        this->setEid( G.getEid() );
        this->setPtE( G.getPtE() );

        if ( Dim == 1 )
        {
            PointSetGaussLobatto<Hypercube<1,1>, Order, value_type> gausslobatto( interior );

            final_pts = gausslobatto.points();
        }

        else if ( Order > 2 )
        {
            entities = std::make_pair( RefConv.vertices(), equiVertices() );

            if ( Dim == 2 )
                final_pts = transformPoints<2>( final_pts );

            else
            {
                uint16_type max_entity_dim = 1;

                if ( Order > 3 )
                    max_entity_dim = 2;

                blendAxis();

                for ( uint16_type d = 1; d < max_entity_dim +1 ; ++d )
                {
                    for ( int e = RefConv.entityRange( d ).begin();
                            e < RefConv.entityRange( d ).end();
                            ++e )
                    {
                        points_type pts = toEquilateral( G.pointsBySubEntity( d, e, 0 ), true );

                        points_type coord_bar = toBarycentric( pts );

                        pts += calculateFaceDeformation( coord_bar, entityMap( d,e ), mpl::int_<3>() );

                        final_pts = putInPointset ( final_pts,
                                                    toEquilateral( pts, false ),
                                                    G.interiorRangeById( d, e ) );
                    }
                }

                if ( Order > 3 )
                {
                    final_pts = putInPointset ( final_pts,
                                                transformPoints<3>( G.pointsBySubEntity( 3, 0, 0 ) ),
                                                G.interiorRangeById( 3, 0 ) );
                }
            }
        }

        this->setPoints( final_pts );

        this->setName( "warpblend", Order );
    }

    ~PointSetWarpBlend() {}

private :

    std::pair<points_type, points_type> entities;

    std::map<uint16_type, std::pair<vector_type, vector_type > > axis;

    points_type equiVertices ()
    {
        points_type V ( ublas::scalar_matrix<value_type>( Dim, Dim+1, value_type( 0 ) ) );

        for ( uint16_type i=0; i < 3; i++ )
        {
            value_type angle = M_PI*( value_type( 7 ) + value_type( 4 )*value_type( i ) )/value_type( 6 );

            V( 0,i ) = math::cos( angle );
            V( 1,i ) = math::sin( angle );

            if ( Dim == 3 )
                V( 2,i ) = - math::sqrt( value_type( 2 ) )/value_type( 4 );
        }

        V *= value_type( 2 )/math::sqrt( value_type( 3 ) );

        if ( Dim == 3 )
            V( 2,3 ) = math::sqrt( value_type( 6 ) )/value_type( 2 );

        return V;
    }

    void blendAxis()
    {
        std::vector<uint16_type> indices;

        for ( uint16_type i=0; i<4; i++ )
        {
            for ( uint16_type j=0; j<4; j++ )
            {
                if ( j != ( 18 +i*( -6 + ( i-1 )*( -3 + 4*( i-2 ) ) ) )/6 )
                    indices.push_back( j );
            }

            axis[i].first  = getVertex( 1,indices[1] ) - getVertex( 1,indices[0] );
            axis[i].second = getVertex( 1,indices[2] ) - ( getVertex( 1,indices[1] ) + getVertex( 1,indices[0] ) )/value_type( 2 );

            axis[i].first  = axis[i].first/ublas::norm_2( axis[i].first );
            axis[i].second = axis[i].second/ublas::norm_2( axis[i].second );

            indices.clear();
        }
    }

    //correspondence between the faces and edges of the reference and the equilateral simplex used
    uint16_type entityMap ( uint16_type top_dim, uint16_type id )
    {
        if ( top_dim == 2 )
            id = ( 12 + id*( 6 + ( id-1 )*( -9 + 4*( id-2 ) ) ) )/6;

        else
        {
            if ( id == 5 )
                id = 2;

            else if ( id > 2 )
                id = 1;

            else
                id = 0;
        }

        return id;
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

    points_type toEquilateral ( points_type pts, bool toEqui )
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

        if ( toEqui )
        {
            pts = ublas::prod( A, pts );

            for ( uint16_type i=0; i < pts.size2(); i++ )
                ublas::column( pts, i ) += ublas::column( b,0 );
        }

        else
        {
            for ( uint16_type i=0; i < pts.size2(); i++ )
                ublas::column( pts, i ) -= ublas::column( b,0 );

            LU< points_type > lu( A );
            pts = lu.solve( pts ) ;
        }

        return pts;
    }

    points_type toBarycentric ( points_type pts )
    {
        points_type C( Dim, Dim );

        for ( uint16_type i = 0; i < Dim; i++ )
            ublas::column( C, i ) = getVertex( 1, ( Dim - Dim%2 + i )%Dim ) - getVertex( 1,Dim );

        for ( uint16_type i=0; i < pts.size2(); i++ )
            ublas::column( pts, i ) -= getVertex( 1,Dim );

        LU< points_type > lu( C );

        points_type solution = lu.solve( pts );

        pts.resize( entities.second.size2(), pts.size2() );

        ublas::subrange( pts, 1, Dim+1, 0, pts.size2() ) = solution;

        ublas::row( pts, 0 ) = ublas::scalar_vector<value_type>( pts.size2(), value_type( 1 ) );

        for ( uint16_type i = 1; i < Dim+1; i++ )
            ublas::row( pts, 0 ) -= ublas::row( pts, i );

        return pts;
    }

    vector_type gll()
    {
        vector_type gx( Order+1 );
        vector_type gw( Order+1 );

        details::dyna::gausslobattojacobi<value_type, vector_type, vector_type>( Order+1, gw, gx );

        return gx;
    }

    value_type alpha( uint16_type d ) const
    {
        value_type p;

        double __alpha_2d[ 16 ] = { 0,
                                    0,
                                    0,
                                    1.4152,
                                    0.1001,
                                    0.2751,
                                    0.9808,
                                    1.0999,
                                    1.2832,
                                    1.3648,
                                    1.4773,
                                    1.4959,
                                    1.5743,
                                    1.5770,
                                    1.6223,
                                    1.6258
                                  };

        double __alpha_3d[ 16 ] = { 0,
                                    0,
                                    0,
                                    0,
                                    0.1002,
                                    1.1332,
                                    1.5608,
                                    1.3413,
                                    1.2577,
                                    1.1603,
                                    1.0153,
                                    0.6080,
                                    0.4523,
                                    0.8856,
                                    0.8717,
                                    0.9655
                                  };

        if ( d == 2 )
            p = ( Order<=15 )?__alpha_2d[ Order ]:( 5./3. );

        else
            p = ( Order<=15 )?__alpha_3d[ Order ]:1.0;

        return p;
    }

    template<typename AE>
    vector_type warpFactor( vector_type const& x, ublas::vector_expression<AE> const& xout ) const
    {
        vector_type warp( xout().size() );
        warp = ublas::scalar_vector<value_type>( xout().size(), value_type( 0 ) );
        vector_type xeq( glas::linspace( value_type( -1 ), value_type( 1 ), nPoints1D, 0 ) );

        vector_type d( xout().size() );

        for ( int i = 0; i < nPoints1D; ++i )
        {
            d = ublas::scalar_vector<value_type>( xout().size(), x( i )-xeq( i ) );

            for ( int j = 1; j < nPoints1D-1; ++j )
                if ( i != j )
                {
                    d = ( ublas::element_prod( d, xout ) - d *xeq( j ) )/( xeq( i ) - xeq( j ) );
                }

            // deflate end roots
            if ( i != 0 )
                d *= value_type( -1 )/( xeq( i ) - xeq( 0 ) );

            if ( i != nPoints1D-1 )
                d *= value_type( 1 )/( xeq( i ) - xeq( nPoints1D-1 ) );

            warp += d;
        }

        return warp;
    }

    //retrieves the right barycentric coordinates for each warp function
    points_type getCoordinates( points_type const& pts, uint16_type face_id )
    {
        std::map<uint16_type, std::vector<uint16_type> > faces;

        for ( uint16_type i=0; i<4; i++ )
            for ( uint16_type j=0; j<4; j++ )
            {
                if ( j != i )
                    faces[i].push_back( j );
            }

        for ( uint16_type i=2; i<4; i++ )
        {
            faces[i][1] = faces[i][2];
            faces[i][2] = 1;
        }

        points_type coord ( 3, pts.size2() );

        for ( uint16_type i=0; i<3; i++ )
            ublas::row( coord, i ) = ublas::row( pts, faces[face_id][i] );

        return coord;
    }

    points_type calculateFaceDeformation( points_type const& coord_bar, mpl::int_<2> )
    {
        points_type blend ( 3, coord_bar.size2() );

        for ( uint16_type i = 0; i < 3; i++ )
            ublas::row( blend, i ) = ublas::element_prod( ublas::row( coord_bar, ( i+1 )%3 ), ublas::row( coord_bar, ( i+2 )%3 ) );


        blend *= value_type( 4 );

        points_type warpfactor( 3, coord_bar.size2() );

        for ( uint16_type i = 0; i < 3; i++ )
            ublas::row( warpfactor, i ) = warpFactor( gll(),  ublas::row( coord_bar, ( i+2 )%3 ) - ublas::row( coord_bar, ( i+1 )%3 ) );

        points_type warp( 3, coord_bar.size2() );

        warp = ublas::element_prod( ublas::element_prod( blend, warpfactor ),
                                    ublas::scalar_matrix<value_type>( 3, coord_bar.size2(), value_type( 1 ) ) +
                                    alpha( 2 )*alpha( 2 )*ublas::element_prod( coord_bar, coord_bar ) );

        points_type deformation( ublas::scalar_matrix<value_type>( 2, warp.size2(), value_type( 0 ) ) );

        for ( uint16_type i=0; i < warp.size1(); i++ )
        {
            value_type angle = value_type( 2*i )*M_PI/value_type( 3 );

            ublas::row( deformation, 0 ) += math::cos( angle )*ublas::row( warp, i );
            ublas::row( deformation, 1 ) += math::sin( angle )*ublas::row( warp, i );
        }

        return deformation;
    }

    //calculates the face deformation for one of the faces of the tetrahedra
    points_type calculateFaceDeformation( points_type const& pts, uint16_type face_id, mpl::int_<3> )
    {
        points_type coord_bar;

        coord_bar = getCoordinates( pts, face_id );

        points_type def = calculateFaceDeformation( coord_bar, mpl::int_<2>() );

        points_type w ( 3, pts.size2() );

        for ( uint16_type i=0; i<3; i++ )
            ublas::row( w, i ) = axis[face_id].first( i )*ublas::row( def,0 ) + axis[face_id].second( i )*ublas::row( def,1 );

        return w;
    }

    points_type blendFunction( points_type const& pts )
    {
        points_type blend ( 4, pts.size2() );

        for ( uint16_type i=0; i<4; i++ )
        {
            vector_type aux1 ( ublas::scalar_vector<value_type>( pts.size2(), value_type( 1 ) ) );
            vector_type aux2 = aux1;

            for ( uint16_type j=0; j<4; j++ )
            {
                if ( j != i )
                {
                    aux1 = ublas::element_prod( aux1, ublas::row( pts, j ) );
                    aux2 = ublas::element_prod( aux2, value_type( 2 )*ublas::row( pts, j ) + ublas::row( pts, i ) );
                }
            }

            ublas::row( blend, i ) = ublas::element_div( aux1, aux2 );
        }

        blend *= value_type( 8 );

        //add beta deformation
        blend = ublas::element_prod( blend, ublas::scalar_matrix<value_type>( 4, pts.size2(), value_type( 1 ) ) +
                                     alpha( 3 )*alpha( 3 )*ublas::element_prod( pts, pts ) );

        return blend;
    }

    points_type calculateFaceDeformation( points_type const& coord_bar, mpl::int_<3> )
    {
        points_type blend = blendFunction( coord_bar );

        std::map<uint16_type, points_type > warp;

        for ( uint16_type face_id = 0; face_id < 4; face_id++ )
        {
            warp[face_id].resize( 3, coord_bar.size2() );

            for ( uint16_type i=0; i<3; i++ )
            {
                vector_type aux = ublas::row( calculateFaceDeformation( coord_bar,
                                              entityMap( 2, face_id ),
                                              mpl::int_<3>() ), i );

                ublas::row( warp[face_id], i ) = ublas::element_prod( ublas::row( blend, face_id ), aux );
            }
        }

        for ( uint16_type face_id = 1; face_id < 4; face_id++ )
            warp[0] += warp[face_id];

        return warp[0];
    }

    template<int N>
    points_type transformPoints( points_type pts )
    {
        pts = toEquilateral( pts, true );

        points_type coord_bar = toBarycentric( pts );

        pts += calculateFaceDeformation( coord_bar, mpl::int_<N>() );

        return toEquilateral( pts, false );
    }

    points_type putInPointset ( points_type final_pts, points_type pts, std::pair<uint16_type, uint16_type> position )
    {
        ublas::subrange( final_pts, 0, 3, position.first, position.second ) = pts;

        return final_pts;
    }
};

} // Feel
#endif /* __WarpBlend_H */
