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
   \file quadmapped.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2006-09-26
*/
#ifndef __Quadmapped_H
#define __Quadmapped_H 1

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelpoly/equispaced.hpp>

#include <feel/feelcore/visitor.hpp>

#include <stdexcept>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/glas.hpp>

#include <feel/feelmesh/geoelement.hpp>
#include <feel/feelpoly/im.hpp>



namespace Feel
{
namespace ublas = boost::numeric::ublas;

/*template< class Convex,
          uint16_type Order,
          typename T = double ,
          template<class, uint16_type, class> class PointSetType = Gauss>*/

template<typename PsetType>
class QuadMapped : public PsetType
{
    typedef PsetType super;
public :
    typedef super pointset_type;

    typedef typename super::value_type value_type;

    typedef typename super::convex_type convex_type;
    static const uint32_type Dim = convex_type::nDim;
    static const uint32_type convexOrder = convex_type::nOrder;

    static const bool is_simplex = convex_type::is_simplex;

    typedef typename pointset_type::nodes_type nodes_type;
    typedef typename matrix_node<value_type>::type points_type;



    typedef ublas::vector<uint16_type> permutation_vector_type;
    typedef ublas::mapped_matrix<uint16_type> permutation_matrix_type;

    //typedef mpl::if_< mpl::bool_< is_simplex >, Simplex<Dim, Order, Dim> , Hypercube<Dim, Order, Dim> > conv_order_type;
    typedef convex_type conv_order_type;

#if 0
    static const uint32_type nbPtsPerVertex = conv_order_type::type::nbPtsPerVertex;
    static const uint32_type nbPtsPerEdge = conv_order_type::type::nbPtsPerEdge;
    static const uint32_type nbPtsPerFace = conv_order_type::type::nbPtsPerFace;
#else
    static const uint32_type nbPtsPerVertex = conv_order_type::nbPtsPerVertex;
    static const uint32_type nbPtsPerEdge = conv_order_type::nbPtsPerEdge;
    static const uint32_type nbPtsPerFace = conv_order_type::nbPtsPerFace;
#endif
    typedef typename convex_type::vertex_permutation_type vertex_permutation_type;
    typedef typename convex_type::edge_permutation_type edge_permutation_type;
    typedef typename convex_type::face_permutation_type face_permutation_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<Dim>, mpl::int_<2> >, mpl::identity<edge_permutation_type>, mpl::identity<face_permutation_type> >::type::type permutation_type;
    typedef std::vector<std::map<permutation_type, points_type> > permutation_points_type;

    QuadMapped()
        :
        super()
    {
    }

    ~QuadMapped() {}

    permutation_points_type
    operator()( pointset_type const& /*pts*/ )
    {
        permutation_points_type ppts( convex_type::numTopologicalFaces );

        // loop over all permutations on the reference convex and
        // store all permutations of the pointset
        for ( int f = 0; f < convex_type::numTopologicalFaces; ++f )
        {
            __typeof__(  typename pointset_type::Face( *this, f ) ) facequad = typename pointset_type::Face( *this, f );

            generate( ppts, facequad, f, mpl::int_<Dim>() );
        }

        return ppts;
    }

private:

    permutation_vector_type getVectorPermutation ( face_permutation_type P )
    {
        return vector_permutation[P];
    }

    permutation_matrix_type getMatrixPermutation ( face_permutation_type P )
    {
        return matrix_permutation[P];
    }

    void
    generate( permutation_points_type& ppts, typename pointset_type::Face const& pts, int f, mpl::int_<1> )
    {
        vertex_permutation_type pbegin( vertex_permutation_type::IDENTITY );
        vertex_permutation_type pend( vertex_permutation_type::N_PERMUTATIONS );

        for ( vertex_permutation_type p = pbegin; p < pend; ++p )
        {
            ppts[f].insert( std::make_pair( p, pts.points() ) );
            DVLOG(2) << "[quadmapped] ppts[" << f << "]["
                          << p.value() << "]="
                          << ppts[f].find( p )->second << "\n";
        }
    }
    void
    generate( permutation_points_type& ppts, typename pointset_type::Face const& pts, int f, mpl::int_<2> )
    {
        edge_permutation_type pbegin( edge_permutation_type::IDENTITY );
        edge_permutation_type pend( edge_permutation_type::N_PERMUTATIONS );

        for ( edge_permutation_type p = pbegin; p < pend; ++p )
        {
            ppts[f].insert( std::make_pair( p, permutatePoints( pts.points(), p ) ) );
            DVLOG(2) << "[quadmapped] ppts[" << f << "][" << p.value() << "]=" << ppts[f].find( p )->second << "\n";
        }
    }
    void
    generate( permutation_points_type& ppts, typename pointset_type::Face const& pts, int f, mpl::int_<3> )
    {
        uint16_type npts = pts.nPoints();

        // generate all permutations
        generateFacePermutations( npts, mpl::bool_<is_simplex>() );
        //generateFacePermutations( npts, mpl::bool_<false>() );
        face_permutation_type pbegin( face_permutation_type::IDENTITY );
        face_permutation_type pend( face_permutation_type::N_PERMUTATIONS );

        for ( face_permutation_type p = pbegin; p < pend; ++p )
        {
            ppts[f].insert( std::make_pair( p, permutatePoints( pts.points(), p ) ) );
            DVLOG(2) << "[quadmapped] ppts[" << f << "][" << p.value() << "]=" << ppts[f].find( p )->second << "\n";
        }
    }

    // Permutates points in the reference element that, when mapped to the
    // real one, match others generated in the same real edge
    points_type permutatePoints ( points_type Gt,
                                  edge_permutation_type edge_perm )
    {
        if ( edge_perm != edge_permutation_type( edge_permutation_type::IDENTITY ) )
        {
            // this formula works with 1-point quadratures
            for ( int i=0; i <= ( int( Gt.size2() )-1 )/2 - int( Gt.size2() )%2; i++ )
                ublas::column( Gt, i ).swap( ublas::column( Gt, Gt.size2()-1-i ) );
        }

        return Gt;
    }

    // Permutates points in the reference element that, when mapped to the
    // real one, match others generated in the same real face
    points_type permutatePoints ( nodes_type const& Gt,
                                  face_permutation_type face_perm )
    {
        FEELPP_ASSERT( face_perm != face_permutation_type( face_permutation_type::NO_PERMUTATION ) )
        ( face_perm ).error( "invalid permutation" );

        nodes_type res( Gt );

        if ( face_perm != face_permutation_type( face_permutation_type::IDENTITY ) )
            res = prod( Gt, getMatrixPermutation( face_perm ) );

        FEELPP_ASSERT( res.size1() == Gt.size1() &&
                       res.size2() == Gt.size2() ) ( res )( Gt ).error ( "invalid permutation operation" );
        return res;
    }

    permutation_matrix_type vectorToMatrixPermutation ( permutation_vector_type const& v )
    {
        permutation_matrix_type P ( v.size(), v.size(), v.size()*v.size() );
        P.clear();

        for ( uint16_type i = 0; i < v.size(); ++i )
            P( i,v( i ) ) = 1;

        return P;
    }

    permutation_vector_type matrixToVectorPermutation ( permutation_matrix_type const& P )
    {
        FEELPP_ASSERT( P.size1() == P.size2() ).error( "invalid permutation" );

        permutation_vector_type v ( P.size1() );
        v.clear();

        for ( uint16_type i = 0; i < v.size(); ++i )
            for ( uint16_type j = 0; j < v.size(); ++j )
                if ( P( i,j ) == 1 )
                    v( i ) = j;

        return v;
    }

    void setPermutation( face_permutation_type out, face_permutation_type first, face_permutation_type second )
    {
        matrix_permutation[out] = ublas::prod( matrix_permutation[first],
                                               matrix_permutation[second] );

        vector_permutation[out] = matrixToVectorPermutation( matrix_permutation[out] );
    }


    //Permutations for tetrahedra
    void generateFacePermutations( uint16_type npoints, mpl::bool_<true> )
    {
        //uint16_type npoints = n_side_points*(n_side_points+1)/2;
        //uint16_type n_side_points = (uint16_type) (-1 + math::sqrt( (double)8*npoints + 1 ))/2;
        uint16_type n_side_points = ( uint16_type )math::sqrt( ( double )npoints );
        //FEELPP_ASSERT( npoints == (n_side_points-1)*(n_side_points-2)/2 )
        //FEELPP_ASSERT( npoints == (n_side_points)*(n_side_points))
        //( npoints )( n_side_points ).error( "invalid number of points" );
        permutation_vector_type _vec( npoints );
        _vec.clear();

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

        setPermutation( ( face_permutation_type )face_permutation_type::ROTATION_ANTICLOCK,
                        ( face_permutation_type )face_permutation_type::REVERSE_BASE,
                        ( face_permutation_type )face_permutation_type::REVERSE_HYPOTENUSE );

        setPermutation( ( face_permutation_type )face_permutation_type::ROTATION_CLOCKWISE,
                        ( face_permutation_type )face_permutation_type::REVERSE_HYPOTENUSE,
                        ( face_permutation_type )face_permutation_type::REVERSE_BASE );

        setPermutation( ( face_permutation_type )face_permutation_type::REVERSE_HEIGHT,
                        ( face_permutation_type )face_permutation_type::REVERSE_BASE,
                        ( face_permutation_type )face_permutation_type::ROTATION_CLOCKWISE );
    }

    //Permutations for hexahedra
    void generateFacePermutations( uint16_type npoints, mpl::bool_<false> )
    {
        //uint16_type npoints = n_side_points*(n_side_points);
        uint16_type n_side_points = ( uint16_type )math::sqrt( ( double )npoints );
        FEELPP_ASSERT( npoints == n_side_points*n_side_points )
        ( npoints )( n_side_points ).error( "invalid number of points" );
        permutation_vector_type _vec( npoints );
        _vec.clear();
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

        vector_permutation[( face_permutation_type )face_permutation_type::REVERSE_BASE] = _vec;
        matrix_permutation[( face_permutation_type )face_permutation_type::REVERSE_BASE] = vectorToMatrixPermutation( _vec );

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

        vector_permutation[( face_permutation_type )face_permutation_type::ROTATION_ANTICLOCK] = _vec;
        matrix_permutation[( face_permutation_type )face_permutation_type::ROTATION_ANTICLOCK] = vectorToMatrixPermutation( _vec );

        setPermutation( ( face_permutation_type )face_permutation_type::SECOND_DIAGONAL,
                        ( face_permutation_type )face_permutation_type::REVERSE_BASE,
                        ( face_permutation_type )face_permutation_type::ROTATION_ANTICLOCK );

        setPermutation( ( face_permutation_type )face_permutation_type::REVERSE_HEIGHT,
                        ( face_permutation_type )face_permutation_type::SECOND_DIAGONAL,
                        ( face_permutation_type )face_permutation_type::ROTATION_ANTICLOCK );

        setPermutation( ( face_permutation_type )face_permutation_type::ROTATION_TWICE_CLOCKWISE,
                        ( face_permutation_type )face_permutation_type::REVERSE_HEIGHT,
                        ( face_permutation_type )face_permutation_type::REVERSE_BASE );

        setPermutation( ( face_permutation_type )face_permutation_type::PRINCIPAL_DIAGONAL,
                        ( face_permutation_type )face_permutation_type::REVERSE_HEIGHT,
                        ( face_permutation_type )face_permutation_type::ROTATION_ANTICLOCK );

        setPermutation( ( face_permutation_type )face_permutation_type::ROTATION_CLOCKWISE,
                        ( face_permutation_type )face_permutation_type::ROTATION_ANTICLOCK,
                        ( face_permutation_type )face_permutation_type::ROTATION_TWICE_CLOCKWISE );
    }

private:

    std::map<face_permutation_type, permutation_vector_type> vector_permutation;
    std::map<face_permutation_type, permutation_matrix_type> matrix_permutation;



};

} // Feel
#endif /* __QuadMapped_H */
