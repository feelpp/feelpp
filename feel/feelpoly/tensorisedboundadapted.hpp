/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
       Date: 2005-12-13

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
   \file tensorizedboundadapted.hpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2005-12-13
 */
#ifndef __TensorisedBoundadapted_H
#define __TensorisedBoundadapted_H 1


#include <boost/lambda/if.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feelalg/glas.hpp>

#include <feel/feelpoly/principal.hpp>

namespace Feel
{

/**
 * \class TensorisedBoundaryAdapted
 * \brief TensorisedBoundaryAdapted Basis on simplex products
 *
 * This class represents the Boundary Adapted Basis made from
 * Jacobi polynomials up to degree \c
 * Degree on a simplex in dimension \c Dim.
 *
 *
 * The Boundary adapted basis is constructed to preserve a part
 * of the Jacobi polynomials' orthogonality. However we need
 * to modify the basis in order to manage easily the boundary condtions.
 *
 * \ingroup Polynomial
 * @author Gilles Steiner
 * @see dubiner.hpp
 */

template<uint16_type Dim,
         uint16_type Degree,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
class TensorisedBoundaryAdapted
{
public:

    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Degree;
    static const uint16_type nConvexOrder = nDim+nOrder+1;

    /** @name Typedefs
     */
    //@{

    typedef TensorisedBoundaryAdapted<Dim, Degree, T, StoragePolicy> self_type;

    /*
     * for now can be used only as a basis but we might be interested
     * to have then expressed in other basis like in the moment basis
     */
    typedef self_type basis_type;

    /*
     * numerical type
     */
    typedef T value_type;
    typedef int16_type int_type;
    typedef uint16_type uint_type;

    /*
     * Geometry where the polynomials are defined and constructed
     */
    typedef Hypercube<nDim, nConvexOrder, nDim> convex_type;
    typedef Reference<convex_type, nDim, nConvexOrder, nDim, T> reference_convex_type;
    typedef typename reference_convex_type::points_type points_type;

    static const uint16_type numVertices = reference_convex_type::numVertices;
    static const uint16_type numFaces = reference_convex_type::numFaces;

    template<uint16_type order>
    struct Convex
    {
        typedef Hypercube<nDim, order, nDim> type;
    };

    /*
     * storage policy
     */
    typedef StoragePolicy<value_type> storage_policy;
    typedef typename storage_policy::vector_type vector_type;
    typedef typename storage_policy::matrix_type matrix_type;
    typedef typename storage_policy::vector_matrix_type vector_matrix_type;
    typedef typename storage_policy::vector_vector_matrix_type  vector_vector_matrix_type;
    typedef typename storage_policy::matrix_node_type matrix_node_type;
    typedef typename storage_policy::node_type node_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    TensorisedBoundaryAdapted()
        :
        M_refconvex(),
        M_pts( nDim, numVertices ),
        M_pts_face( numVertices )
    {
        PointSetEquiSpaced<convex_type, nOrder, value_type> pts;

        // only the points associated with the vertices
        M_pts = pts.pointsByEntity( 0 );

        // get the points for each faces
        for ( uint16_type e = M_refconvex.entityRange( nDim-1 ).begin();
                e < M_refconvex.entityRange( nDim-1 ).end();
                ++e )
        {
            M_pts_face[e] = pts.pointsBySubEntity( nDim-1, e, 1 );
        }

        matrix_type A( ublas::trans( evaluate( pts.points() ) ) );
        matrix_type D = ublas::identity_matrix<value_type>( A.size1(), A.size2()  );
        LU<matrix_type> lu( A );
        matrix_type C = lu.solve( D );
        vector_matrix_type d ( derivate( pts.points() ) );
        M_D.resize( d.size() );

        for ( size_type i = 0; i < d.size(); ++i )
        {
            M_D[i] = ublas::prod( d[i], C );
            glas::clean( M_D[i] );
        }
    }

    TensorisedBoundaryAdapted( TensorisedBoundaryAdapted const & d )
        :
        M_refconvex( d.M_refconvex ),
        M_pts( d.M_pts ),
        M_pts_face( d.M_pts_face )
    {
    }

    ~TensorisedBoundaryAdapted()
    {
    }

    //@}

    /**
     * Access to the points of the reference convex associated
     **/

    points_type points()
    {
        return M_pts;
    }

    /**
     * Access to the points associated with the face \c f
     **/
    points_type const& points( int f ) const
    {
        return M_pts_face[f];
    }

    /** @name Operator overloads
     */
    //@{

    self_type const& operator=( self_type const& d )
    {
        if ( this != &d )
        {
            M_pts = d.M_pts;
            M_D = d.M_D;
        }

        return *this;
    }

    matrix_type operator()( node_type const& pt ) const
    {
        points_type pts( pt.size(), 1 );
        ublas::column( pts, 0 ) = pt;
        return evaluate( pts );
    }

    matrix_type operator()( points_type const& pts ) const
    {
        return evaluate( pts );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the maximum degree of the TensorisedBoundaryAdapted polynomial to be
     * constructed
     */
    uint_type degree() const
    {
        return nOrder;
    }

    /**
     * \return self as a basis
     */
    self_type const& basis() const
    {
        return *this;
    }

    /**
     * \return the family name of the finite element
     */
    std::string familyName() const
    {
        return "tensorizedboundaryadapted";
    }


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    matrix_type coeff() const
    {
        return ublas::identity_matrix<value_type>( reference_convex_type::polyDims( nOrder ), M_pts.size2() );
    }


    /**
     * evaluate the TensorisedBoundaryAdapted polynomials at a set of points \p __pts
     *
     * \arg __pts is a set of points
     */
    matrix_type evaluate( points_type const& __pts ) const
    {
        return evaluate( __pts, mpl::int_<nDim>() );
    }

    template<typename AE>
    vector_matrix_type derivate( ublas::matrix_expression<AE>  const& __pts ) const
    {
        return derivate( __pts, mpl::int_<nDim>() );
    }

    /**
     * \brief derivatives of Dubiner polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    matrix_type const& d( uint16_type i )
    {
        return M_D[i];
    }

    /**
     * \brief derivatives of Dubiner polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    matrix_type const& derivate( uint16_type i )
    {
        return M_D[i];
    }

    //@}

private:

    /**
     * Evaluation at a set of points of the expansion basis in 1D
     */

    /**
     * Here we compute the modified basis in one dimension following
     * Sherwin-Karniadakis description (p.111), i.e.
     * \f[ \psi_i(x) = \left\{\begin{array}{ll} \frac{1-x}{2}, & i=0, \\ \frac{1-x}{2}\frac{1+x}{2} P^{(1,1)}_{i-1}(x), & 1 \leq i < N, \\ \frac{1+x}{2}, & i=N, \\ \end{array}\right.   \f]
     *
     */
    matrix_type
    evaluate( points_type const& __pts, mpl::int_<1> ) const
    {
        return M_pfunc.evaluate_1( ublas::row( __pts,0 ) );
    }

    /**
     * derivation at a set of points of the expansion basis in 1D
     */
    template<typename AE>
    vector_matrix_type
    derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<1> ) const
    {
        FEELPP_ASSERT( __pts().size1() == 1 )( __pts().size1() )( __pts().size2() ).error( "invalid points" );

        vector_matrix_type D( 1 );
        D[0].resize( nOrder+1, __pts().size2() );
        D[0] = M_pfunc.derivate_1( ublas::row( __pts(),0 ) );

        return D;
    }

    /**
     * Evaluation at a set of points of the expansion basis in 2D on
     * the square
     */
    matrix_type evaluate( points_type const& __pts, mpl::int_<2> ) const;

    /**
     * derivation at a set of points of the expansion basis in 2D on
     * the square
     */
    template<typename AE>
    vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<2> ) const;

    /**
     * Evaluation at a set of points of the expansion basis in 3D on
     * the cube
     */
    matrix_type evaluate( points_type const& __pts, mpl::int_<3> ) const;

    /**
     * derivation at a set of points of the expansion basis in 3D on
     * the cube
     */
    template<typename AE>
    vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<3> ) const;

private:
    reference_convex_type M_refconvex;
    points_type M_pts;
    std::vector<points_type> M_pts_face;

    Principal<Degree, T, StoragePolicy> M_pfunc;

    vector_matrix_type M_D;
};

template<uint16_type Dim,
         uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename TensorisedBoundaryAdapted<Dim, Degree,  T, StoragePolicy>::matrix_type
TensorisedBoundaryAdapted<Dim, Degree,  T, StoragePolicy>::evaluate( points_type const& __pts, mpl::int_<2> ) const
{
    matrix_type res( convex_type::polyDims( nOrder ), __pts.size2() );

    vector_type eta1s = ublas::row( __pts, 0 );
    vector_type eta2s = ublas::row( __pts, 1 );

    matrix_type psi_1( M_pfunc.evaluate_1( eta1s ) );
    matrix_type psi_2( M_pfunc.evaluate_1( eta2s ) );

    /*

    C ------- D    y
    |         |    ^
    |         |    |
    |         |    |
    |         |    |
    |         |    |-----> x
    A ------- B

    */

    /* 4 Vertex */

    /* A */
    ublas::row( res,0 ) = ublas::element_prod( ublas::row( psi_1,0 )      , ublas::row( psi_2,0 ) );
    /* B */
    ublas::row( res,1 ) = ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2,0 ) );
    /* C */
    ublas::row( res,2 ) = ublas::element_prod( ublas::row( psi_1,0 )      , ublas::row( psi_2,nOrder ) );
    /* D */
    ublas::row( res,3 ) = ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2,nOrder ) );

    /* Global index */
    uint_type G_i = 3;

    /* Edges */
    //@{
    /* AB = (nOrder-1) */
    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) = ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,0 ) );
    }

    /* AC = (nOrder-1)*/
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( res,G_i ) = ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2,q ) );

        if ( q%2==0 )
            ublas::row( res,G_i )*= value_type( -1.0 );
    }

    /* CD = (nOrder-1) */
    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) = ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,nOrder ) );

        if ( p%2==0 )
            ublas::row( res,G_i )*= value_type( -1.0 );
    }

    /* BD = (nOrder-1)*/
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( res,G_i ) = ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2,q ) );
    }

    //@}

    /* Interior : (N-1)*(N-1) */
    //@{
    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res,G_i ) = ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,q ) );
        }

    //@}

    return res;
}


template<uint16_type Dim,
         uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename TensorisedBoundaryAdapted<Dim, Degree,  T, StoragePolicy>::vector_matrix_type
TensorisedBoundaryAdapted<Dim, Degree,  T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<2> ) const
{
    vector_matrix_type res( 2 );
    res[0].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[1].resize( convex_type::polyDims( nOrder ), __pts().size2() );

    vector_type eta1s = ublas::row( __pts(), 0 );
    vector_type eta2s = ublas::row( __pts(), 1 );

    matrix_type psi_1( M_pfunc.evaluate_1( eta1s ) );
    matrix_type psi_2( M_pfunc.evaluate_1( eta2s ) );

    matrix_type dpsi_1( M_pfunc.derivate_1( eta1s ) );
    matrix_type dpsi_2( M_pfunc.derivate_1( eta2s ) );

    /* 4 Vertex */

    /* A */
    ublas::row( res[0],0 ) = ublas::element_prod( ublas::row( dpsi_1,0 )      , ublas::row( psi_2,0 ) );
    ublas::row( res[1],0 ) = ublas::element_prod( ublas::row( psi_1,0 )      , ublas::row( dpsi_2,0 ) );
    /* B */
    ublas::row( res[0],1 ) = ublas::element_prod( ublas::row( dpsi_1,nOrder ) , ublas::row( psi_2,0 ) );
    ublas::row( res[1],1 ) = ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( dpsi_2,0 ) );
    /* C */
    ublas::row( res[0],2 ) = ublas::element_prod( ublas::row( dpsi_1,0 )      , ublas::row( psi_2,nOrder ) );
    ublas::row( res[1],2 ) = ublas::element_prod( ublas::row( psi_1,0 )      , ublas::row( dpsi_2,nOrder ) );
    /* D */
    ublas::row( res[0],3 ) = ublas::element_prod( ublas::row( dpsi_1,nOrder ) , ublas::row( psi_2,nOrder ) );
    ublas::row( res[1],3 ) = ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( dpsi_2,nOrder ) );

    /* Global index */
    uint_type G_i = 3;

    /* Edges */
    //@{
    /* AB = (nOrder-1) */
    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2,0 ) );
        ublas::row( res[1],G_i ) = ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2,0 ) );
    }

    /* AC = (nOrder-1)*/
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( dpsi_1,0 ) , ublas::row( psi_2,q ) );
        ublas::row( res[1],G_i ) = ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( dpsi_2,q ) );
    }

    /* CD = (nOrder-1) */
    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2,nOrder ) );
        ublas::row( res[1],G_i ) = ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2,nOrder ) );
    }

    /* BD = (nOrder-1)*/
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( dpsi_1,nOrder ) , ublas::row( psi_2,q ) );
        ublas::row( res[1],G_i ) = ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( dpsi_2,q ) );
    }

    //@}

    /* Interior : (N-1)*(N-1) */
    //@{
    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2,q ) );
            ublas::row( res[1],G_i ) = ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2,q ) );
        }

    //@}

    return res;
}

template<uint16_type Dim,
         uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename TensorisedBoundaryAdapted<Dim, Degree,  T, StoragePolicy>::matrix_type
TensorisedBoundaryAdapted<Dim, Degree,  T, StoragePolicy>::evaluate( points_type const& __pts, mpl::int_<3> ) const
{
    matrix_type res( convex_type::polyDims( nOrder ), __pts.size2() );

    ublas::vector<value_type> eta1s = ublas::row( __pts, 0 );
    ublas::vector<value_type> eta2s = ublas::row( __pts, 1 );
    ublas::vector<value_type> eta3s = ublas::row( __pts, 2 );

    matrix_type psi_1( M_pfunc.evaluate_1( eta1s ) );
    matrix_type psi_2( M_pfunc.evaluate_1( eta2s ) );
    matrix_type psi_3( M_pfunc.evaluate_1( eta3s ) );

    /*

    H---------G
    /|        /|
    / |       / |
    /  |      /  |
    E---------F   |        z    y
    |   |     |   |        ^   ^
    |   D-----|---C        |  /
    |  /      |  /         | /
    | /       | /          |/
    |/        |/           |-----> x
    A---------B

    */


    /* 8 Vertex */

    /* A */
    ublas::row( res,0 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,0 ) );
    /* B */
    ublas::row( res,1 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,0 ) );
    /* C */
    ublas::row( res,2 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,0 ) );
    /* D */
    ublas::row( res,3 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,0 ) );
    /* E */
    ublas::row( res,4 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,nOrder ) );
    /* F */
    ublas::row( res,5 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,nOrder ) );
    /* G */
    ublas::row( res,6 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,nOrder ) );
    /* H */
    ublas::row( res,7 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,nOrder ) );


    /* Global index */
    uint_type G_i = 7;

    /* 12 Edges */
    //@{
    /* AB : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,0 ) );
    }

    /* BC : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,p ) ), ublas::row( psi_3,0 ) );
    }

    /* DC : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,0 ) );
    }

    /* AD : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,p ) ), ublas::row( psi_3,0 ) );
    }

    /* AE : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,p ) );
    }

    /* BF : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,p ) );
    }

    /* CG : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,p ) );
    }

    /* DH : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,p ) );
    }

    /* EF : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,nOrder ) );
    }

    /* FG : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,p ) ), ublas::row( psi_3,nOrder ) );
    }

    /* HG : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,nOrder ) );
    }

    /* EH : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,p ) ), ublas::row( psi_3,nOrder ) );
    }

    //@}


    /* 6 Faces : (N-1)(N-1) modes */
    //@{

    /* ABCD */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res,G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,q ) ), ublas::row( psi_3,0 ) );
        }

    /* ABFE */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res,G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,0 ) ), ublas::row( psi_3,q ) );
        }

    /* BCGF */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res,G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2,p ) ), ublas::row( psi_3,q ) );
        }

    /* CGHD */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res,G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,q ) );
        }

    /* ADHE */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res,G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2,p ) ), ublas::row( psi_3,q ) );
        }

    /* EFGH */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res,G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,q ) ), ublas::row( psi_3,nOrder ) );
        }

    //@}


    /* Interior : (N-1)^3 */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
            for ( uint_type r = 1; r < nOrder; ++r )
            {
                ++G_i;
                ublas::row( res,G_i ) =
                    ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,q ) ), ublas::row( psi_3,r ) );
            }

    return res;
}

template<uint16_type Dim,
         uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename TensorisedBoundaryAdapted<Dim, Degree,  T, StoragePolicy>::vector_matrix_type
TensorisedBoundaryAdapted<Dim, Degree,  T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<3> ) const
{
    vector_matrix_type res( 3 );
    res[0].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[1].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[2].resize( convex_type::polyDims( nOrder ), __pts().size2() );


    vector_type eta1s = ublas::row( __pts(), 0 );
    vector_type eta2s = ublas::row( __pts(), 1 );
    vector_type eta3s = ublas::row( __pts(), 2 );

    matrix_type psi_1( M_pfunc.evaluate_1( eta1s ) );
    matrix_type dpsi_1( M_pfunc.derivate_1( eta1s ) );

    matrix_type psi_2( M_pfunc.evaluate_1( eta2s ) );
    matrix_type dpsi_2( M_pfunc.derivate_1( eta2s ) );

    matrix_type psi_3( M_pfunc.evaluate_1( eta3s ) );
    matrix_type dpsi_3( M_pfunc.derivate_1( eta3s ) );

    /* 8 Vertex */

    ublas::row( res[0],0 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,0 ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,0 ) );
    ublas::row( res[1],0 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( dpsi_2,0 ) ), ublas::row( psi_3,0 ) );
    ublas::row( res[2],0 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,0 ) ), ublas::row( dpsi_3,0 ) );

    ublas::row( res[0],1 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,nOrder ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,0 ) );
    ublas::row( res[1],1 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( dpsi_2,0 ) ), ublas::row( psi_3,0 ) );
    ublas::row( res[2],1 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,0 ) ), ublas::row( dpsi_3,0 ) );

    ublas::row( res[0],2 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,nOrder ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,0 ) );
    ublas::row( res[1],2 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( dpsi_2,nOrder ) ), ublas::row( psi_3,0 ) );
    ublas::row( res[2],2 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,nOrder ) ), ublas::row( dpsi_3,0 ) );

    ublas::row( res[0],3 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,0 ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,0 ) );
    ublas::row( res[1],3 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( dpsi_2,nOrder ) ), ublas::row( psi_3,0 ) );
    ublas::row( res[2],3 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,nOrder ) ), ublas::row( dpsi_3,0 ) );

    ublas::row( res[0],4 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,0 ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,nOrder ) );
    ublas::row( res[1],4 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( dpsi_2,0 ) ), ublas::row( psi_3,nOrder ) );
    ublas::row( res[2],4 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,0 ) ), ublas::row( dpsi_3,nOrder ) );

    ublas::row( res[0],5 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,nOrder ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,nOrder ) );
    ublas::row( res[1],5 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( dpsi_2,0 ) ), ublas::row( psi_3,nOrder ) );
    ublas::row( res[2],5 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,0 ) ), ublas::row( dpsi_3,nOrder ) );

    ublas::row( res[0],6 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,nOrder ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,nOrder ) );
    ublas::row( res[1],6 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( dpsi_2,nOrder ) ), ublas::row( psi_3,nOrder ) );
    ublas::row( res[2],6 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,nOrder ) ), ublas::row( dpsi_3,nOrder ) );

    ublas::row( res[0],7 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,0 ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,nOrder ) );
    ublas::row( res[1],7 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( dpsi_2,nOrder ) ), ublas::row( psi_3,nOrder ) );
    ublas::row( res[2],7 ) =
        ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,nOrder ) ), ublas::row( dpsi_3,nOrder ) );


    /* Global index */
    uint_type G_i = 7;

    /* 12 Edges */
    //@{
    /* AB : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,p ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,0 ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( dpsi_2,0 ) ), ublas::row( psi_3,0 ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( psi_2,0 ) ), ublas::row( dpsi_3,0 ) );
    }

    /* BC : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,nOrder ), ublas::row( psi_2,p ) ), ublas::row( psi_3,0 ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( dpsi_2,p ) ), ublas::row( psi_3,0 ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,p ) ), ublas::row( dpsi_3,0 ) );
    }

    /* DC : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,p ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,0 ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( dpsi_2,nOrder ) ), ublas::row( psi_3,0 ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( psi_2,nOrder ) ), ublas::row( dpsi_3,0 ) );
    }

    /* AD : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,0 ), ublas::row( psi_2,p ) ), ublas::row( psi_3,0 ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( dpsi_2,p ) ), ublas::row( psi_3,0 ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,p ) ), ublas::row( dpsi_3,0 ) );
    }

    /* AE : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,0 ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,p ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( dpsi_2,0 ) ), ublas::row( psi_3,p ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,0 ) ), ublas::row( dpsi_3,p ) );
    }

    /* BF : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,nOrder ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,p ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( dpsi_2,0 ) ), ublas::row( psi_3,p ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,0 ) ), ublas::row( dpsi_3,p ) );
    }

    /* CG : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,nOrder ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,p ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( dpsi_2,nOrder ) ), ublas::row( psi_3,p ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,nOrder ) ), ublas::row( dpsi_3,p ) );
    }

    /* DH : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,0 ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,p ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( dpsi_2,nOrder ) ), ublas::row( psi_3,p ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,nOrder ) ), ublas::row( dpsi_3,p ) );
    }

    /* EF : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,p ), ublas::row( psi_2,0 ) ), ublas::row( psi_3,nOrder ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( dpsi_2,0 ) ), ublas::row( psi_3,nOrder ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( psi_2,0 ) ), ublas::row( dpsi_3,nOrder ) );
    }

    /* FG : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,nOrder ), ublas::row( psi_2,p ) ), ublas::row( psi_3,nOrder ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( dpsi_2,p ) ), ublas::row( psi_3,nOrder ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ), ublas::row( psi_2,p ) ), ublas::row( dpsi_3,nOrder ) );
    }

    /* HG : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,p ), ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,nOrder ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( dpsi_2,nOrder ) ), ublas::row( psi_3,nOrder ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ), ublas::row( psi_2,nOrder ) ), ublas::row( dpsi_3,nOrder ) );
    }

    /* EH : (nOrder-1) */

    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,0 ), ublas::row( psi_2,p ) ), ublas::row( psi_3,nOrder ) );
        ublas::row( res[1],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( dpsi_2,p ) ), ublas::row( psi_3,nOrder ) );
        ublas::row( res[2],G_i ) =
            ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ), ublas::row( psi_2,p ) ), ublas::row( dpsi_3,nOrder ) );
    }

    //@}


    /* 6 Faces : (N-1)(N-1) modes */
    //@{

    /* ABCD */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res[0],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2,q ) ), ublas::row( psi_3,0 ) );
            ublas::row( res[1],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2,q ) ), ublas::row( psi_3,0 ) );
            ublas::row( res[2],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,q ) ), ublas::row( dpsi_3,0 ) );
        }

    /* ABFE */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res[0],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2,0 ) ), ublas::row( psi_3,q ) );
            ublas::row( res[1],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2,0 ) ), ublas::row( psi_3,q ) );
            ublas::row( res[2],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,0 ) ), ublas::row( dpsi_3,q ) );
        }

    /* BCGF */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res[0],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,nOrder ) , ublas::row( psi_2,p ) ), ublas::row( psi_3,q ) );
            ublas::row( res[1],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( dpsi_2,p ) ), ublas::row( psi_3,q ) );
            ublas::row( res[2],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2,p ) ), ublas::row( dpsi_3,q ) );
        }

    /* CGHD */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res[0],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2,nOrder ) ), ublas::row( psi_3,q ) );
            ublas::row( res[1],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2,nOrder ) ), ublas::row( psi_3,q ) );
            ublas::row( res[2],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,nOrder ) ), ublas::row( dpsi_3,q ) );
        }

    /* ADHE */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res[0],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,0 ) , ublas::row( psi_2,p ) ), ublas::row( psi_3,q ) );
            ublas::row( res[1],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( dpsi_2,p ) ), ublas::row( psi_3,q ) );
            ublas::row( res[2],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2,p ) ), ublas::row( dpsi_3,q ) );
        }

    /* EFGH */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
        {
            ++G_i;
            ublas::row( res[0],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2,q ) ), ublas::row( psi_3,nOrder ) );
            ublas::row( res[1],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2,q ) ), ublas::row( psi_3,nOrder ) );
            ublas::row( res[2],G_i ) =
                ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,q ) ), ublas::row( dpsi_3,nOrder ) );
        }

    //@}


    /* Interior : (N-1)^3 */

    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder; ++q )
            for ( uint_type r = 1; r < nOrder; ++r )
            {
                ++G_i;
                ublas::row( res[0],G_i ) =
                    ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2,q ) ), ublas::row( psi_3,r ) );
                ublas::row( res[1],G_i ) =
                    ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2,q ) ), ublas::row( psi_3,r ) );
                ublas::row( res[2],G_i ) =
                    ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2,q ) ), ublas::row( dpsi_3,r ) );
            }


    return res;
}




}
#endif /* __BoundAdapted_H */
