/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

This file is part of the Feel library

Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2006-03-06

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
   \file mapped.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2006-09-26
*/
#ifndef __PointSetMapped_H
#define __PointSetMapped_H 1

#include <stdexcept>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>


#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/visitor.hpp>
#include <feel/feelcore/traits.hpp>

#include <feel/feelalg/glas.hpp>

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/geoelement.hpp>

#include <feel/feelpoly/equispaced.hpp>





namespace Feel
{
namespace ublas = boost::numeric::ublas;

template< typename element_type,
          class Convex,
          uint16_type Order,
          typename T = double ,
          template<class, uint16_type, class> class PointSetType = PointSetEquiSpaced>
class PointSetMapped : public PointSetType<Convex, Order, T>
{

public :
    typedef PointSetType<Convex, Order, T> pointset_type;

    typedef T value_type;

    static const uint32_type Dim = Convex::nDim;
    static const uint32_type convexOrder = Convex::nOrder;

    static const bool is_simplex = Convex::is_simplex;

    typedef typename pointset_type::nodes_type nodes_type;
    typedef typename matrix_node<value_type>::type points_type;

    typedef typename element_type::gm_type gm_type;
    typedef typename element_type::gm_ptrtype gm_ptrtype;

    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;

    typedef ublas::vector<uint16_type> permutation_vector_type;
    typedef ublas::mapped_matrix<uint16_type> permutation_matrix_type;

    typedef mpl::if_< mpl::bool_< is_simplex >,
            Simplex<Dim, Order, Dim> ,
            Hypercube<Dim, Order, Dim> > conv_order_type;

    static const uint32_type nbPtsPerVertex = conv_order_type::type::nbPtsPerVertex;
    static const uint32_type nbPtsPerEdge = conv_order_type::type::nbPtsPerEdge;
    static const uint32_type nbPtsPerFace = conv_order_type::type::nbPtsPerFace;

    typedef Reference<Convex, Dim, convexOrder, Dim, value_type> RefElem;

    typedef typename pointset_type::range_type range_type;
    typedef typename pointset_type::index_map_type index_map_type;

    RefElem RefConv;

    PointSetMapped( element_type const& _elt )
        :
        M_elt( _elt )
    {
        pointset_type pts;

        points_type Gt = updatePoints<2>( updatePoints<1>( pts.points(), pts, mpl::bool_< ( Order > 2 ) >() ),
                                          pts,
                                          mpl::bool_< ( Dim == 3 && Order > 2+uint16_type( is_simplex ) ) >() );

        this->setPoints( Gt );
    }

    ~PointSetMapped() {}

    permutation_vector_type getVectorPermutation ( face_permutation_type P )
    {
        return vector_permutation[P];
    }

    permutation_matrix_type getMatrixPermutation ( face_permutation_type P )
    {
        return matrix_permutation[P];
    }

    points_type pointsBySubEntity( uint16_type top_dim, uint16_type local_id, bool boundary = 0, bool real = 0 )
    {
        index_map_type index_list = this->entityToLocal( top_dim, local_id, boundary );

        //Rearrange the order of the vertices, edges and faces

        uint16_type matrix_size = 0;

        if ( index_list[0].size() != 0 )
            matrix_size +=index_list[0].size()*nbPtsPerVertex;

        if ( ( top_dim>=1 ) && ( index_list[1].size() != 0 ) )
            matrix_size +=index_list[1].size()*nbPtsPerEdge;

        if ( ( top_dim>=2 ) && ( index_list[2].size() != 0 ) )
            matrix_size +=index_list[2].size()*nbPtsPerFace;

        points_type G ( Dim, matrix_size );

        for ( uint16_type i=0, p=0; i < top_dim+1; i++ )
        {
            if ( index_list[i].size() )
            {
                for ( uint16_type j=0; j < index_list[i].size(); j++ )
                {
                    //Something else needs to be done here: depending on the entity, and in the case
                    //we have boundary=1, the points in the boundary have to be reordered

                    points_type aux = this->interiorPointsById( i, index_list[i][j] );

                    ublas::subrange( G, 0, Dim, p, p+aux.size2() ) = aux;

                    p+=aux.size2();
                }
            }
        }

        if ( real )
            G = transformToReal( G );

        return G;
    }

private:

    element_type M_elt;

    std::map<face_permutation_type, permutation_vector_type> vector_permutation;
    std::map<face_permutation_type, permutation_matrix_type> matrix_permutation;

    template <int top_dim>
    nodes_type updatePoints( nodes_type const& pts, pointset_type const& /*G*/, mpl::bool_<false> )
    {
        return pts;
    }

    template <int top_dim>
    nodes_type updatePoints( nodes_type const& thepts, pointset_type const& G, mpl::bool_<true> )
    {
        nodes_type pts( thepts );

        if ( top_dim == 2 )
        {
            generateFacePermutations( Order - 1-uint16_type( is_simplex ),
                                      mpl::int_<Dim>(), mpl::bool_< is_simplex >() );
        }

        for ( size_type e = RefConv.entityRange( top_dim ).begin();
                e < RefConv.entityRange( top_dim ).end();
                ++e )
        {
            points_type Gt = G.pointsBySubEntity( top_dim, e );

            Gt = permutatePoints ( Gt, e, mpl::int_<top_dim>() );

            uint16_type pos = G.interiorRangeById( top_dim, e ).first;

            ublas::subrange( pts, 0, Dim, pos, pos+Gt.size2() ) = Gt;
        }

        return pts;
    }

    // Permutates points in the reference element that, when mapped to the
    // real one, match others generated in the same real edge
    points_type permutatePoints ( nodes_type const& theGt, size_type entity_local_id, mpl::int_<1> )
    {
        nodes_type Gt( theGt );
        edge_permutation_type permutation = M_elt.edgePermutation( entity_local_id );

        if ( permutation != edge_permutation_type( edge_permutation_type::IDENTITY ) )
        {
            for ( uint16_type i=0; i <= ( Gt.size2()-1 )/2 - Gt.size2()%2; i++ )
                ublas::column( Gt, i ).swap( ublas::column( Gt, Gt.size2()-1-i ) );
        }

        return Gt;
    }

    // Permutates points in the reference element that, when mapped to the
    // real one, match others generated in the same real face
    points_type permutatePoints ( nodes_type const& Gt, size_type entity_local_id, mpl::int_<2> )
    {
        face_permutation_type permutation = M_elt.facePermutation( entity_local_id );
        nodes_type res( Gt );

        if ( permutation != face_permutation_type( face_permutation_type::IDENTITY ) )
            res = prod( Gt, getMatrixPermutation( permutation ) );

        return res;
    }

    permutation_matrix_type vectorToMatrixPermutation ( permutation_vector_type const& v )
    {
        permutation_matrix_type P ( v.size(), v.size(), v.size()*v.size() );

        for ( uint16_type i = 0; i < v.size(); ++i )
            P( i,v( i ) ) = 1;

        return P;
    }

    permutation_vector_type matrixToVectorPermutation ( permutation_matrix_type const& P )
    {
        FEELPP_ASSERT( P.size1() == P.size2() ).error( "invalid permutation" );

        permutation_vector_type v ( P.size1() );

        for ( uint16_type i = 0; i < v.size(); ++i )
            for ( uint16_type j = 0; j < v.size(); ++j )
                if ( P( i,j ) == 1 )
                    v( i ) = j;

        return v;
    }

    void setPermutation( int out, int first, int second )
    {
        matrix_permutation[face_permutation_type( out )] = ublas::prod( matrix_permutation[face_permutation_type( first )],
                matrix_permutation[face_permutation_type( second )] );

        vector_permutation[face_permutation_type( out )] = matrixToVectorPermutation( matrix_permutation[face_permutation_type( out )] );
    }

    void generateFacePermutations( uint16_type /*n_side_points*/, mpl::int_<2>, mpl::bool_<true> ) {}
    void generateFacePermutations( uint16_type /*n_side_points*/, mpl::int_<2>, mpl::bool_<false> ) {}

    //Permutations for tetrahedra
    void generateFacePermutations( uint16_type n_side_points, mpl::int_<3>, mpl::bool_<true> )
    {
        uint16_type npoints = n_side_points*( n_side_points+1 )/2;

        permutation_vector_type _vec( npoints );

        //define hypotenuse reverse permutation
        for ( uint16_type i = 0; i < n_side_points; ++i )
        {
            for ( uint16_type j = 0; j <= i; j++ )
            {
                uint16_type _first = i + j*( j-1 )/2 + j*( n_side_points - j );
                uint16_type _last = i + ( i-j )*( i-j-1 )/2 + ( i-j )*( n_side_points - ( i-j ) );

                _vec( _first ) = _last;
            }
        }

        vector_permutation[face_permutation_type( face_permutation_type::REVERSE_HYPOTENUSE )] = _vec;
        matrix_permutation[face_permutation_type( face_permutation_type::REVERSE_HYPOTENUSE )] = vectorToMatrixPermutation( _vec );

        //define base reverse permutation
        for ( uint16_type i = 0; i < n_side_points; ++i )
        {
            uint16_type _begin = i*n_side_points - i*( i-1 )/2;
            uint16_type _end = ( i+1 )*n_side_points - ( i+1 )*i/2 - 1;

            for ( uint16_type j = 0; j <= n_side_points - i - 1; j++ )
                _vec( _begin + j ) = _end - j;
        }

        vector_permutation[face_permutation_type( face_permutation_type::REVERSE_BASE )] = _vec;
        matrix_permutation[face_permutation_type( face_permutation_type::REVERSE_BASE )] = vectorToMatrixPermutation( _vec );

        setPermutation( face_permutation_type::ROTATION_ANTICLOCK,
                        face_permutation_type::REVERSE_BASE,
                        face_permutation_type::REVERSE_HYPOTENUSE );

        setPermutation( face_permutation_type::ROTATION_CLOCKWISE,
                        face_permutation_type::REVERSE_HYPOTENUSE,
                        face_permutation_type::REVERSE_BASE );

        setPermutation( face_permutation_type::REVERSE_HEIGHT,
                        face_permutation_type::REVERSE_BASE,
                        face_permutation_type::ROTATION_CLOCKWISE );
    }

    //Permutations for hexahedra
    void generateFacePermutations( uint16_type n_side_points, mpl::int_<3>, mpl::bool_<false> )
    {
        uint16_type npoints = n_side_points*( n_side_points );

        permutation_vector_type _vec( npoints );

        //define base permutation (tau_2)
        uint16_type p=0;

        for ( int16_type i = n_side_points-1; i >= 0; --i )
        {
            uint16_type _first = i*n_side_points;

            for ( uint16_type j = 0; j <= n_side_points-1; j++ )
            {
                _vec( p ) = _first + j;
                p++;
            }
        }

        vector_permutation[face_permutation_type( face_permutation_type::REVERSE_BASE )] = _vec;
        matrix_permutation[face_permutation_type( face_permutation_type::REVERSE_BASE )] = vectorToMatrixPermutation( _vec );

        //define once anticlockwise permutation (tau_3)
        p=0;

        for ( int16_type i = n_side_points-1; i >= 0; --i )
        {
            for ( uint16_type j = 0; j <= n_side_points-1; j++ )
            {
                _vec( p ) = i + n_side_points*j;
                p++;
            }
        }

        vector_permutation[face_permutation_type( face_permutation_type::ROTATION_ANTICLOCK )] = _vec;
        matrix_permutation[face_permutation_type( face_permutation_type::ROTATION_ANTICLOCK )] = vectorToMatrixPermutation( _vec );

        setPermutation( face_permutation_type::SECOND_DIAGONAL,
                        face_permutation_type::REVERSE_BASE,
                        face_permutation_type::ROTATION_ANTICLOCK );

        setPermutation( face_permutation_type::REVERSE_HEIGHT,
                        face_permutation_type::SECOND_DIAGONAL,
                        face_permutation_type::ROTATION_ANTICLOCK );

        setPermutation( face_permutation_type::ROTATION_TWICE_CLOCKWISE,
                        face_permutation_type::REVERSE_HEIGHT,
                        face_permutation_type::REVERSE_BASE );

        setPermutation( face_permutation_type::PRINCIPAL_DIAGONAL,
                        face_permutation_type::REVERSE_HEIGHT,
                        face_permutation_type::ROTATION_ANTICLOCK );

        setPermutation( face_permutation_type::ROTATION_CLOCKWISE,
                        face_permutation_type::ROTATION_ANTICLOCK,
                        face_permutation_type::ROTATION_TWICE_CLOCKWISE );
    }

    points_type transformToReal( points_type const& Gt )
    {
        gm_ptrtype _gm_ptr( new gm_type );

        typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( _gm_ptr, Gt ) );

        typename gm_type::template Context<vm::POINT, element_type> gmc( _gm_ptr, M_elt, __geopc );

        return gmc.xReal();
    }
};

} // Feel
#endif /* __PointSetMapped_H */
