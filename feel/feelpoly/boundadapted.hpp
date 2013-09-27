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
   \file boundadapted.hpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2005-12-13
 */
#ifndef __Boundadapted_H
#define __Boundadapted_H 1


#include <boost/lambda/if.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feelalg/glas.hpp>

#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/principal.hpp>
#include <feel/feelpoly/basis.hpp>

namespace Feel
{
template<class Convex,uint16_type O,typename T> class PointSetWarpBlend;
template<uint16_type Dim,uint16_type RealDim,uint16_type Degree,typename NormalizationPolicy,typename T,template<class> class StoragePolicy> class Dubiner;

template<uint16_type Dim,
         uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
class BoundaryAdapted;

template<uint16_type Dim,
         uint16_type Degree,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
struct BoundaryAdaptedTraits
{
    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Degree;
    static const uint16_type nConvexOrderDiff = nDim+nOrder+1;
    static const bool is_normalized = false;


    typedef T value_type;

    template<uint16_type order, typename V = value_type>
    struct Convex
    {
        typedef Simplex<nDim, order, nDim> type;
        typedef Reference<Simplex<nDim, order, nDim>, nDim, order, nDim, V>  reference_type;
    };

    template<typename NewT>
    struct ChangeValueType
    {
        typedef BoundaryAdapted<Dim, Degree, NewT, StoragePolicy> type;
        typedef BoundaryAdaptedTraits<Dim, Degree, NewT, StoragePolicy> traits_type;

    };

    template<uint16_type NewOrder>
    struct ChangeOrder
    {
        typedef BoundaryAdapted<Dim, NewOrder, T, StoragePolicy> type;
        typedef BoundaryAdaptedTraits<Dim, NewOrder, T, StoragePolicy> traits_type;
    };

    /*
     * Geometry where the polynomials are defined and constructed
     */

    typedef typename Convex<nOrder>::type convex_type;
    typedef typename Convex<nOrder>::reference_type reference_convex_type;
    typedef typename Convex<nConvexOrderDiff>::type diff_convex_type;
    typedef typename Convex<nConvexOrderDiff>::reference_type diff_reference_convex_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
            mpl::identity<PointSetWarpBlend<diff_convex_type,nConvexOrderDiff,value_type> >,
            mpl::identity<PointSetEquiSpaced<diff_convex_type, nConvexOrderDiff,value_type> > >::type::type diff_pointset_type;

    static const uint16_type numVertices = reference_convex_type::numVertices;
    static const uint16_type numFaces = reference_convex_type::numFaces;

    /*
     * storage policy
     */

    typedef StoragePolicy<value_type> storage_policy;
    typedef typename storage_policy::vector_type vector_type;
    typedef typename storage_policy::matrix_type matrix_type;
    typedef typename storage_policy::vector_matrix_type vector_matrix_type;
    typedef typename storage_policy::vector_vector_matrix_type vector_vector_matrix_type;
    typedef typename storage_policy::matrix_node_type matrix_node_type;
    typedef typename storage_policy::points_type points_type;
    typedef typename storage_policy::node_type node_type;
};

template<int D, int O>
struct BoundaryAdaptedTag
{
    static const int Dim = D;
    static const int Order = O;
};

/**
 * \class BoundaryAdapted
 * \brief BoundaryAdapted Basis on simplex
 *
 * This class represents the Boundary Adapted Basis made from
 * Dubiner polynomials up to degree \c
 * Degree on a simplex in dimension \c Dim.
 *
 *
 * The Boundary adapted basis is construct to preserve a part
 * of the Dubiner polynomials' orthogonality. However we need
 * to modify the basis in order to manage easily the boundary condtions.
 *
 * \ingroup Polynomial
 * @author Gilles Steiner
 * @see boundadapted.hpp
 */

template<uint16_type Dim,
         uint16_type Degree,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
class BoundaryAdapted : Basis<BoundaryAdaptedTag<Dim, Degree>, T>
{
    typedef Basis<BoundaryAdaptedTag<Dim, Degree>, T > super;

public:
    typedef BoundaryAdaptedTraits<Dim, Degree, T, StoragePolicy> traits_type;

    static const uint16_type nDim = traits_type::nDim;
    static const uint16_type nOrder = traits_type::nOrder;
    static const uint16_type nConvexOrderDiff = traits_type::nConvexOrderDiff;
    static const bool is_normalized = traits_type::is_normalized;
    static const uint16_type numVertices = traits_type::numVertices;
    static const uint16_type numFaces = traits_type::numFaces;

    /** @name Typedefs
     */
    //@{

    typedef BoundaryAdapted<Dim, Degree, T, StoragePolicy> self_type;

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
    typedef typename traits_type::convex_type convex_type;
    typedef typename traits_type::reference_convex_type reference_convex_type;
    typedef typename traits_type::diff_pointset_type diff_pointset_type;

    /*
     * storage policy
     */
    typedef typename traits_type::storage_policy storage_policy;
    typedef typename traits_type::vector_type vector_type;
    typedef typename traits_type::matrix_type matrix_type;
    typedef typename traits_type::vector_matrix_type vector_matrix_type;
    typedef typename traits_type::vector_vector_matrix_type vector_vector_matrix_type;
    typedef typename traits_type::matrix_node_type matrix_node_type;
    typedef typename traits_type::points_type points_type;
    typedef typename traits_type::node_type node_type;

    typedef Principal<Degree, T,StoragePolicy> principal_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    BoundaryAdapted()
        :
        super( *this ),
        M_refconvex(),
        M_pts( nDim, numVertices ),
        M_pts_face( reference_convex_type::numVertices )
    {
        PointSetEquiSpaced<convex_type,nOrder,value_type> pts;

        // only the points associated with the vertices
        M_pts = pts.pointsByEntity( 0 );
        DVLOG(2) << "[boundaryadapted] pts= " <<  M_pts << "\n";

        // get the points for each faces
        for ( uint16_type e = M_refconvex.entityRange( nDim-1 ).begin();
                e < M_refconvex.entityRange( nDim-1 ).end();
                ++e )
        {
            M_pts_face[e] = pts.pointsBySubEntity( nDim-1, e, 1 );
            DVLOG(2) << "[boundaryadapted] face " << e << " pts " <<  M_pts_face[e] << "\n";
        }
    }

    BoundaryAdapted( BoundaryAdapted const & d )
        :
        super( *this ),
        M_refconvex( d.M_refconvex ),
        M_pts( d.M_pts ),
        M_pts_face( d.M_pts_face )
    {}

    ~BoundaryAdapted()
    {}

    //@}

    /**
     * Access to the points of the reference convex associated
     **/
    points_type const& points() const
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
            M_refconvex = d.M_refconvex;
            M_pts = d.M_pts;
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
     * \return the maximum degree of the BoundaryAdapted polynomial to be
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
     * \return the \c familyName()
     */
    std::string familyName() const
    {
        return "dubiner.boundary.adapted";
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
     * evaluate the BoundaryAdapted polynomials at a set of points \p __pts
     *
     * \arg __pts is a set of points
     */
    static matrix_type evaluate( points_type const& __pts )
    {
        return evaluate( __pts, mpl::int_<nDim>() );
    }

    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE>  const& __pts )
    {
        return derivate( __pts, mpl::int_<nDim>() );
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
    static evaluate( points_type const& __pts, mpl::int_<1> )
    {
        matrix_type E = M_pfunc.evaluate_1( ublas::row( __pts,0 ) );
        matrix_type D;
        D.resize( E.size1(), E.size2() );

        ublas::row( D, 0 ) = ublas::row( E, 0 );
        ublas::row( D, 1 ) = ublas::row( E, E.size2()-1 );

        for ( unsigned int i=1; i < E.size2()-1 ; ++i )
        {
            ublas::row( D, i+1 ) = ublas::row( E, i );
        }

        return D;
    }

    /**
     * derivation at a set of points of the expansion basis in 1D
     */
    template<typename AE>
    vector_matrix_type
    static derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<1> )
    {
        FEELPP_ASSERT( __pts().size1() == 1 )( __pts().size1() )( __pts().size2() ).error( "invalid points" );

        vector_matrix_type E( 1 );
        E[0].resize( nOrder+1, __pts().size2() );
        E[0] = M_pfunc.derivate_1( ublas::row( __pts(),0 ) );

        vector_matrix_type D( 1 );
        D[0].resize( nOrder+1, __pts().size2() );

        ublas::row( D[0], 0 ) = ublas::row( E[0], 0 );
        ublas::row( D[0], 1 ) = ublas::row( E[0], E[0].size2()-1 );

        for ( unsigned int i=1; i < E[0].size2()-1 ; ++i )
        {
            ublas::row( D[0], i+1 ) = ublas::row( E[0], i );
        }

        return D;
    }

    /**
     * Evaluation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    static matrix_type evaluate( points_type const& __pts, mpl::int_<2> );

    /**
     * derivation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<2> );

    /**
     * Evaluation at a set of points of the expansion basis in 3D on
     * the tetrahedron
     */
    static matrix_type evaluate( points_type const& __pts, mpl::int_<3> );

    /**
     * derivation at a set of points of the expansion basis in 3D on
     * the tetrahedron
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<3> );

private:
    reference_convex_type M_refconvex;
    points_type M_pts;
    std::vector<points_type> M_pts_face;
    static principal_type M_pfunc;

};
template<uint16_type Dim,
         uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename BoundaryAdapted<Dim, Degree,  T, StoragePolicy>::principal_type
BoundaryAdapted<Dim, Degree,  T, StoragePolicy>::M_pfunc;

template<uint16_type Dim,
         uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename BoundaryAdapted<Dim, Degree,  T, StoragePolicy>::matrix_type
BoundaryAdapted<Dim, Degree,  T, StoragePolicy>::evaluate( points_type const& __pts, mpl::int_<2> )
{
    matrix_type res( convex_type::polyDims( nOrder ), __pts.size2() );

    details::etas<TRIANGLE, value_type> etas( __pts );
    vector_type eta1s = ublas::row( etas(), 0 );
    vector_type eta2s = ublas::row( etas(), 1 );

    matrix_type psi_1( M_pfunc.evaluate_1( eta1s ) );
    vector_matrix_type psi_2( M_pfunc.evaluate_2( eta2s ) );

    /* 3 Vertex */

    ublas::row( res,0 ) = ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2[0],0 ) );
    ublas::row( res,1 ) = ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2[nOrder],0 ) );

    vector_type ones( ublas::scalar_vector<value_type>( __pts.size2(), value_type( 1.0 ) ) );

    ublas::row( res,2 ) = ( ones + eta2s ) / 2.0;

    /* Global index */
    uint_type G_i = 2;

    /* Edges */
    //@{

    /* BD = \Gamma_0 : (nOrder-1)*/
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( res,G_i ) = ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2[nOrder],q ) );
    }

    /* AC = \Gamma_1 : (nOrder-1)*/
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( res,G_i ) = ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2[0],q ) );

        if ( q%2==0 )
            ublas::row( res,G_i )*= value_type( -1.0 );
    }

    /* AB = \Gamma_2 : (nOrder-1) */
    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) = ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2[p],0 ) );
    }

    //@}


    /* Interior : (nOrder-1)(nOrder-2)/2 */
    //@{
    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder-p; ++q )
        {
            ++G_i;
            ublas::row( res,G_i ) = ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2[p],q ) );
        }

    //@}

    return res;
}


template<uint16_type Dim,
         uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename BoundaryAdapted<Dim, Degree,  T, StoragePolicy>::vector_matrix_type
BoundaryAdapted<Dim, Degree,  T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<2> )
{
    vector_matrix_type res( 2 );
    res[0].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[1].resize( convex_type::polyDims( nOrder ), __pts().size2() );

    details::etas<TRIANGLE, value_type> etas( __pts );
    vector_type eta1s = ublas::row( etas(), 0 );
    vector_type eta2s = ublas::row( etas(), 1 );

    matrix_type psi_1( M_pfunc.evaluate_1( eta1s ) );
    vector_matrix_type psi_2( M_pfunc.evaluate_2( eta2s ) );

    matrix_type dpsi_1( M_pfunc.derivate_1( eta1s ) );
    vector_matrix_type dpsi_2( M_pfunc.derivate_2( eta2s ) );

    vector_type ones( ublas::scalar_vector<value_type>( __pts().size2(), value_type( 1.0 ) ) );
    vector_type d1( value_type( 2.0 )*ublas::element_div( ones,ones-eta2s ) ); // 2 / (1 - xi_2)
    vector_type d2( ublas::element_div( ones+eta1s ,ones-eta2s ) ); // 2(1+xi_1) / (1 - xi_2)^2

    /* 3 Vertex */

    ublas::row( res[0],0 ) = ublas::element_prod( ublas::row( dpsi_1,0 ) , ublas::row( psi_2[0],0 ) );
    ublas::row( res[1],0 ) = ublas::element_prod( ublas::row( res[0],0 ) , d2 );

    ublas::row( res[0],0 ) = ublas::element_prod( ublas::row( res[0],0 ) , d1 );
    ublas::row( res[1],0 ) += ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( dpsi_2[0],0 ) );

    ublas::row( res[0],1 ) = ublas::element_prod( ublas::row( dpsi_1,nOrder ) , ublas::row( psi_2[nOrder],0 ) );
    ublas::row( res[1],1 ) = ublas::element_prod( ublas::row( res[0],1 ) , d2 );

    ublas::row( res[0],1 ) = ublas::element_prod( ublas::row( res[0],1 ) , d1 );
    ublas::row( res[1],1 ) += ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( dpsi_2[nOrder],0 ) );

    ublas::row( res[0],2 ) = ublas::scalar_vector<value_type>( __pts().size2(),value_type( 0.0 ) );
    ublas::row( res[1],2 ) = value_type( 0.5 ) * ones;

    /* Global index */
    uint_type G_i = 2;

    /* Edges */
    //@{

    /* BD : (nOrder-1)*/
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( dpsi_1,nOrder ) , ublas::row( psi_2[nOrder],q ) );
        ublas::row( res[1],G_i ) = ublas::element_prod( ublas::row( res[0],G_i ) , d2 );

        ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( res[0],G_i ) , d1 );
        ublas::row( res[1],G_i ) += ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( dpsi_2[nOrder],q ) );
    }

    /* AC : (nOrder-1)*/
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( dpsi_1,0 ) , ublas::row( psi_2[0],q ) );
        ublas::row( res[1],G_i ) = ublas::element_prod( ublas::row( res[0],G_i ) , d2 );

        ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( res[0],G_i ) , d1 );
        ublas::row( res[1],G_i ) += ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( dpsi_2[0],q ) );

        if ( q%2==0 )
        {
            ublas::row( res[0],G_i )*= value_type( -1.0 );
            ublas::row( res[1],G_i )*= value_type( -1.0 );
        }
    }

    /* AB : (nOrder-1) */
    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2[p],0 ) );
        ublas::row( res[1],G_i ) = ublas::element_prod( ublas::row( res[0],G_i ) , d2 );

        ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( res[0],G_i ) , d1 );
        ublas::row( res[1],G_i ) += ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2[p],0 ) );
    }

    /* Interior : (nOrder-1)(nOrder-2)/2 */
    //@{
    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder-p; ++q )
        {
            ++G_i;
            ublas::row( res[0],G_i ) =  ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2[p],q ) );
            ublas::row( res[1],G_i ) =  ublas::element_prod( ublas::row( res[0],G_i ) , d2 );

            ublas::row( res[0],G_i ) = ublas::element_prod( ublas::row( res[0],G_i ) , d1 );
            ublas::row( res[1],G_i ) +=  ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2[p],q ) );
        }

    //@}


    return res;
}

template<uint16_type Dim,
         uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename BoundaryAdapted<Dim, Degree,  T, StoragePolicy>::matrix_type
BoundaryAdapted<Dim, Degree,  T, StoragePolicy>::evaluate( points_type const& __pts, mpl::int_<3> )
{
    matrix_type res( convex_type::polyDims( nOrder ), __pts.size2() );

    FEELPP_ASSERT( __pts.size1() == 3 )( __pts.size1() ).error( "invalid space dimension" );

    details::etas<TETRAHEDRON, value_type> etas( __pts );
    vector_type eta1s = ublas::row( etas(), 0 );
    vector_type eta2s = ublas::row( etas(), 1 );
    vector_type eta3s = ublas::row( etas(), 2 );

    matrix_type psi_1( M_pfunc.evaluate_1( eta1s ) );
    vector_matrix_type psi_2( M_pfunc.evaluate_2( eta2s ) );
    vector_vector_matrix_type psi_3( M_pfunc.evaluate_3( eta3s ) );

    ublas::scalar_vector<value_type> ones( eta3s.size(), value_type( 1.0 ) );

    /* 4 Vertex */

    ublas::row( res,0 ) =
        ublas::element_prod(
            ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2[0],0 ) ),
            ublas::row( psi_3[0][0],0 )
        );
    ublas::row( res,1 ) =
        ublas::element_prod(
            ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2[nOrder],0 ) ),
            ublas::row( psi_3[nOrder][0],0 )
        );

    ublas::row( res,2 ) =
        ublas::element_prod(
            ublas::row( psi_2[0],nOrder ), ublas::row( psi_3[0][nOrder],0 )
        );

    ublas::row( res,3 ) = ( ones + eta3s ) / 2.0;


    /* Global index */
    uint_type G_i = 3;

    /* Edges */
    //@{
    /* AB : (nOrder-1) */
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2[nOrder],q ) ),
                ublas::row( psi_3[nOrder][q],0 )
            );
    }

    /* AC : (nOrder-1)*/
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2[0],q ) ),
                ublas::row( psi_3[0][q],0 )
            );
        /* Ensure Consistent anticlockwise orientation of the
           edges */

        if ( q%2==0 )
        {
            ublas::row( res,G_i )*= value_type( -1.0 );
        }
    }

    /* BC : (nOrder-1)*/
    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2[p],0 ) ),
                ublas::row( psi_3[p][0],0 )
            );
    }

    /* AD : (nOrder-1)*/
    for ( uint_type r = 1; r < nOrder; ++r )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2[0],0 ) ),
                ublas::row( psi_3[0][0],r )
            );
    }

    /* BD : (nOrder-1)*/
    for ( uint_type r = 1; r < nOrder; ++r )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2[nOrder],0 ) ),
                ublas::row( psi_3[nOrder][0],r )
            );
    }

    /* CD : (nOrder-1)*/
    for ( uint_type r = 1; r < nOrder; ++r )
    {
        ++G_i;
        ublas::row( res,G_i ) =
            ublas::element_prod( ublas::row( psi_2[0],nOrder ),ublas::row( psi_3[0][nOrder],r ) );
    }

    //@}


    /* Faces  */
    //@{

    /* BCE : (N-1)(N-2)/2 */
    for ( uint_type q = 1; q < nOrder; ++q )
        for ( uint_type r = 1; r < nOrder-q; ++r )
        {
            ++G_i;
            ublas::row( res,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2[nOrder],q ) ),
                    ublas::row( psi_3[nOrder][q],r )
                );
        }

    /* ACE : (N-1)(N-2)/2 */
    for ( uint_type q = 1; q < nOrder; ++q )
        for ( uint_type r = 1; r < nOrder-q; ++r )
        {
            ++G_i;
            ublas::row( res,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2[0],q ) ),
                    ublas::row( psi_3[0][q],r )
                );
        }

    /* ABE : (N-1)(N-2)/2 */
    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type r = 1; r < nOrder-p; ++r )
        {
            ++G_i;
            ublas::row( res,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2[p],0 ) ),
                    ublas::row( psi_3[p][0],r )
                );
        }


    /* ABC : (N-1)(N-2)/2 */
    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder-p; ++q )
        {
            ++G_i;
            ublas::row( res,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2[p],q ) ),
                    ublas::row( psi_3[p][q],0 )
                );
        }

    //@}


    /* Interior : (N-1)(N-2)(N-3)/6 */

    for ( uint_type p = 1; p < nOrder-2; ++p )
        for ( uint_type q = 1; q < nOrder-1-p; ++q )
            for ( uint_type r = 1; r < nOrder-p-q; ++r )
            {
                ++G_i;
                ublas::row( res,G_i ) =
                    ublas::element_prod(
                        ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2[p],q ) ),
                        ublas::row( psi_3[p][q],r )
                    );
            }

    return res;
}

template<uint16_type Dim,
         uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename BoundaryAdapted<Dim, Degree,  T, StoragePolicy>::vector_matrix_type
BoundaryAdapted<Dim, Degree,  T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<3> )
{
    vector_matrix_type res( 3 );
    res[0].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[1].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[2].resize( convex_type::polyDims( nOrder ), __pts().size2() );

    FEELPP_ASSERT( __pts().size1() == 3 )( __pts().size1() ).error( "invalid space dimension" );

    details::etas<TETRAHEDRON, value_type> etas( __pts );
    vector_type eta1s = ublas::row( etas(), 0 );
    vector_type eta2s = ublas::row( etas(), 1 );
    vector_type eta3s = ublas::row( etas(), 2 );

    matrix_type psi_1( M_pfunc.evaluate_1( eta1s ) );
    matrix_type dpsi_1( M_pfunc.derivate_1( eta1s ) );

    vector_matrix_type psi_2( M_pfunc.evaluate_2( eta2s ) );
    vector_matrix_type dpsi_2( M_pfunc.derivate_2( eta2s ) );

    vector_vector_matrix_type psi_3( M_pfunc.evaluate_3( eta3s ) );
    vector_vector_matrix_type dpsi_3( M_pfunc.derivate_3( eta3s ) );

    ublas::scalar_vector<value_type> ones( __pts().size2(), value_type( 1.0 ) );

    /**
     * We must save the Jacobian matrix \f$ Jac(i,j) = \frac{\partial\eta_j}{\partial \xi_i} \f$
     * for all points and all cases.
     * We choose to build \f$Jac \in \mathcal M_{9 \times pts().size2()} \f$
     * such that ublas::row(Jac,3*(j-1)+i-1) = \f$\frac{\partial\eta_j}{\partial \xi_i}\f$
     **/

    matrix_type Jac( 9, __pts().size2() );

    /** \f$ d\eta_1/d\xi_i \f$ **/

    ublas::row( Jac,0 ) = /** -2/(\xi_2 + \xi_3) **/
        ublas::element_div( ( -2.0 )*ones, ( ublas::row( __pts(), 1 ) + ublas::row( __pts(), 2 ) ) );

    ublas::row( Jac,1 ) = /** 2(1+\xi_1)/(\xi_2 + \xi_3)^2 **/
        ublas::element_div( 2.0*( ones + ublas::row( __pts(), 0 ) ) ,  ( ublas::row( __pts(), 1 ) + ublas::row( __pts(), 2 ) ) );

    ublas::row( Jac,1 ) =
        ublas::element_div( ublas::row( Jac,1 ),  ublas::row( __pts(), 1 )+ublas::row( __pts(), 2 ) );

    ublas::row( Jac,2 ) = /** 2(1+\xi_1)/(\xi_2 + \xi_3)^2 **/
        ublas::row( Jac,1 );

    /** \f$ d\eta_2/d\xi_i \f$ **/

    ublas::row( Jac,3 ) = /** 0 **/
        ublas::scalar_vector<value_type>( __pts().size2(), value_type( 0.0 ) );

    ublas::row( Jac,4 ) = /** 2/(1-\xi_3) **/
        ublas::element_div( 2.0*ones, ( ones - ublas::row( __pts(), 2 ) ) );

    ublas::row( Jac,5 ) = /** 2(1+\xi_2)/(1-\xi_3)^2 **/
        ublas::element_div( 2.0*( ones + ublas::row( __pts(), 1 ) ) ,  ( ones - ublas::row( __pts(), 2 ) ) );

    ublas::row( Jac,5 ) =
        ublas::element_div( ublas::row( Jac,5 ),  ones  - ublas::row( __pts(), 2 ) );

    /** \f$ d\eta_3/d\xi_i \f$ **/

    ublas::row( Jac,6 ) = /** 0 **/
        ublas::row( Jac,3 );

    ublas::row( Jac,7 ) = /** 0 **/
        ublas::row( Jac,3 );

    ublas::row( Jac,8 ) = ones; /** 1 **/


    /** Usefull intermediate matrix to simplify the construction **/

    matrix_type d1 ( res[1].size1(), res[1].size2() );
    matrix_type d2 ( res[1].size1(), res[1].size2() );
    matrix_type d3 ( res[1].size1(), res[1].size2() );

    /* 4 Vertex */

    ublas::row( d1,0 ) =
        ublas::element_prod(
            ublas::element_prod( ublas::row( dpsi_1,0 ) , ublas::row( psi_2[0],0 ) ),
            ublas::row( psi_3[0][0],0 )
        );
    ublas::row( d2,0 ) =
        ublas::element_prod(
            ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( dpsi_2[0],0 ) ),
            ublas::row( psi_3[0][0],0 )
        );
    ublas::row( d3,0 ) =
        ublas::element_prod(
            ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2[0],0 ) ),
            ublas::row( dpsi_3[0][0],0 )
        );

    ublas::row( d1,1 ) =
        ublas::element_prod(
            ublas::element_prod( ublas::row( dpsi_1,nOrder ) , ublas::row( psi_2[nOrder],0 ) ),
            ublas::row( psi_3[nOrder][0],0 )
        );
    ublas::row( d2,1 ) =
        ublas::element_prod(
            ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( dpsi_2[nOrder],0 ) ),
            ublas::row( psi_3[nOrder][0],0 )
        );
    ublas::row( d3,1 ) =
        ublas::element_prod(
            ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2[nOrder],0 ) ),
            ublas::row( dpsi_3[nOrder][0],0 )
        );

    ublas::row( d1,2 ) = ublas::scalar_vector<value_type>( __pts().size2(), value_type( 0.0 ) );

    ublas::row( d2,2 ) =
        ublas::element_prod(
            ublas::row( dpsi_2[0],nOrder ), ublas::row( psi_3[0][nOrder],0 )
        );

    ublas::row( d3,2 ) =
        ublas::element_prod(
            ublas::row( psi_2[0],nOrder ), ublas::row( dpsi_3[0][nOrder],0 )
        );

    ublas::row( d1,3 ) =
        ublas::scalar_vector<value_type>( __pts().size2(), value_type( 0.0 ) );

    ublas::row( d2,3 ) =
        ublas::scalar_vector<value_type>( __pts().size2(), value_type( 0.0 ) );

    ublas::row( d3,3 ) =
        ublas::scalar_vector<value_type>( __pts().size2(), value_type( 0.5 ) );

    /* Global index */
    uint_type G_i = 3;

    /* Edges */
    //@{
    /* AB : (nOrder-1) */
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( d1,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( dpsi_1,nOrder ) , ublas::row( psi_2[nOrder],q ) ),
                ublas::row( psi_3[nOrder][q],0 )
            );
        ublas::row( d2,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( dpsi_2[nOrder],q ) ),
                ublas::row( psi_3[nOrder][q],0 )
            );
        ublas::row( d3,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2[nOrder],q ) ),
                ublas::row( dpsi_3[nOrder][q],0 )
            );
    }

    /* AC : (nOrder-1)*/
    for ( uint_type q = 1; q < nOrder; ++q )
    {
        ++G_i;
        ublas::row( d1,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( dpsi_1,0 ) , ublas::row( psi_2[0],q ) ),
                ublas::row( psi_3[0][q],0 )
            );
        ublas::row( d2,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( dpsi_2[0],q ) ),
                ublas::row( psi_3[0][q],0 )
            );
        ublas::row( d3,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2[0],q ) ),
                ublas::row( dpsi_3[0][q],0 )
            );

        /* Ensure Consistent anticlockwise orientation of the
           edges */
        if ( q%2==0 )
        {
            ublas::row( d1,G_i )*= value_type( -1.0 );
            ublas::row( d2,G_i )*= value_type( -1.0 );
            ublas::row( d3,G_i )*= value_type( -1.0 );
        }
    }

    /* BC : (nOrder-1)*/
    for ( uint_type p = 1; p < nOrder; ++p )
    {
        ++G_i;
        ublas::row( d1,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2[p],0 ) ),
                ublas::row( psi_3[p][0],0 )
            );
        ublas::row( d2,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2[p],0 ) ),
                ublas::row( psi_3[p][0],0 )
            );
        ublas::row( d3,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2[p],0 ) ),
                ublas::row( dpsi_3[p][0],0 )
            );
    }

    /* AD : (nOrder-1)*/
    for ( uint_type r = 1; r < nOrder; ++r )
    {
        ++G_i;
        ublas::row( d1,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( dpsi_1,0 ) , ublas::row( psi_2[0],0 ) ),
                ublas::row( psi_3[0][0],r )
            );
        ublas::row( d2,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( dpsi_2[0],0 ) ),
                ublas::row( psi_3[0][0],r )
            );
        ublas::row( d3,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2[0],0 ) ),
                ublas::row( dpsi_3[0][0],r )
            );
    }

    /* BD : (nOrder-1)*/
    for ( uint_type r = 1; r < nOrder; ++r )
    {
        ++G_i;
        ublas::row( d1,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( dpsi_1,nOrder ) , ublas::row( psi_2[nOrder],0 ) ),
                ublas::row( psi_3[nOrder][0],r )
            );
        ublas::row( d2,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( dpsi_2[nOrder],0 ) ),
                ublas::row( psi_3[nOrder][0],r )
            );
        ublas::row( d3,G_i ) =
            ublas::element_prod(
                ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2[nOrder],0 ) ),
                ublas::row( dpsi_3[nOrder][0],r )
            );
    }

    /* CD : (nOrder-1)*/
    for ( uint_type r = 1; r < nOrder; ++r )
    {
        ++G_i;
        ublas::row( d1,G_i ) = ublas::scalar_vector<value_type>( __pts().size2(), value_type( 0.0 ) );
        ublas::row( d2,G_i ) =
            ublas::element_prod( ublas::row( dpsi_2[0],nOrder ),ublas::row( psi_3[0][nOrder],r ) );
        ublas::row( d3,G_i ) =
            ublas::element_prod( ublas::row( psi_2[0],nOrder ),ublas::row( dpsi_3[0][nOrder],r ) );
    }

    //@}


    /* Faces  */
    //@{

    /* BCE : (N-1)(N-2)/2 */
    for ( uint_type q = 1; q < nOrder; ++q )
        for ( uint_type r = 1; r < nOrder-q; ++r )
        {
            ++G_i;
            ublas::row( d1,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( dpsi_1,nOrder ) , ublas::row( psi_2[nOrder],q ) ),
                    ublas::row( psi_3[nOrder][q],r )
                );
            ublas::row( d2,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( dpsi_2[nOrder],q ) ),
                    ublas::row( psi_3[nOrder][q],r )
                );
            ublas::row( d3,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,nOrder ) , ublas::row( psi_2[nOrder],q ) ),
                    ublas::row( dpsi_3[nOrder][q],r )
                );
        }

    /* ACE : (N-1)(N-2)/2 */
    for ( uint_type q = 1; q < nOrder; ++q )
        for ( uint_type r = 1; r < nOrder-q; ++r )
        {
            ++G_i;
            ublas::row( d1,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( dpsi_1,0 ) , ublas::row( psi_2[0],q ) ),
                    ublas::row( psi_3[0][q],r )
                );
            ublas::row( d2,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( dpsi_2[0],q ) ),
                    ublas::row( psi_3[0][q],r )
                );
            ublas::row( d3,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,0 ) , ublas::row( psi_2[0],q ) ),
                    ublas::row( dpsi_3[0][q],r )
                );
        }

    /* ABE : (N-1)(N-2)/2 */
    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type r = 1; r < nOrder-p; ++r )
        {
            ++G_i;
            ublas::row( d1,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2[p],0 ) ),
                    ublas::row( psi_3[p][0],r )
                );
            ublas::row( d2,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2[p],0 ) ),
                    ublas::row( psi_3[p][0],r )
                );
            ublas::row( d3,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2[p],0 ) ),
                    ublas::row( dpsi_3[p][0],r )
                );
        }


    /* ABC : (N-1)(N-2)/2 */
    for ( uint_type p = 1; p < nOrder; ++p )
        for ( uint_type q = 1; q < nOrder-p; ++q )
        {
            ++G_i;
            ublas::row( d1,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2[p],q ) ),
                    ublas::row( psi_3[p][q],0 )
                );
            ublas::row( d2,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2[p],q ) ),
                    ublas::row( psi_3[p][q],0 )
                );
            ublas::row( d3,G_i ) =
                ublas::element_prod(
                    ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2[p],q ) ),
                    ublas::row( dpsi_3[p][q],0 )
                );
        }

    //@}


    /* Interior : (N-1)(N-2)(N-3)/6 */

    for ( uint_type p = 1; p < nOrder-2; ++p )
        for ( uint_type q = 1; q < nOrder-1-p; ++q )
            for ( uint_type r = 1; r < nOrder-p-q; ++r )
            {
                ++G_i;
                ublas::row( d1,G_i ) =
                    ublas::element_prod( ublas::element_prod( ublas::row( dpsi_1,p ) , ublas::row( psi_2[p],q ) ), ublas::row( psi_3[p][q],r ) );
                ublas::row( d2,G_i ) =
                    ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( dpsi_2[p],q ) ), ublas::row( psi_3[p][q],r ) );
                ublas::row( d3,G_i ) =
                    ublas::element_prod( ublas::element_prod( ublas::row( psi_1,p ) , ublas::row( psi_2[p],q ) ), ublas::row( dpsi_3[p][q],r ) );
            }


    /** Adding Jacobian contribution **/

#if defined (FEELPP_HAS_QD_REAL)

    /* Warning : if qd are enable we can not use
     * the ublas::diagonal_matrix type
     */

    matrix_type j0( __pts().size2(), __pts().size2() );
    matrix_type j1( __pts().size2(), __pts().size2() );
    matrix_type j2( __pts().size2(), __pts().size2() );
    matrix_type j4( __pts().size2(), __pts().size2() );
    matrix_type j5( __pts().size2(), __pts().size2() );

    for ( uint_type i=0; i < __pts().size2(); ++i )
        for ( uint_type j=i; j < __pts().size2(); ++j )
        {
            if ( j==i )
            {
                j0( i,i ) = ublas::row( Jac,0 )( i );
                j1( i,i ) = ublas::row( Jac,1 )( i );
                j2( i,i ) = ublas::row( Jac,2 )( i );
                j4( i,i ) = ublas::row( Jac,4 )( i );
                j5( i,i ) = ublas::row( Jac,5 )( i );
            }

            else
            {
                j0( i,j ) = value_type( 0.0 );
                j1( i,j ) = value_type( 0.0 );
                j2( i,j ) = value_type( 0.0 );
                j4( i,j ) = value_type( 0.0 );
                j5( i,j ) = value_type( 0.0 );
                j0( j,i ) = value_type( 0.0 );
                j1( j,i ) = value_type( 0.0 );
                j2( j,i ) = value_type( 0.0 );
                j4( j,i ) = value_type( 0.0 );
                j5( j,i ) = value_type( 0.0 );
            }
        }

#else
    ublas::diagonal_matrix<value_type> j0( __pts().size2() );
    ublas::diagonal_matrix<value_type> j1( __pts().size2() );
    ublas::diagonal_matrix<value_type> j2( __pts().size2() );
    ublas::diagonal_matrix<value_type> j4( __pts().size2() );
    ublas::diagonal_matrix<value_type> j5( __pts().size2() );

    for ( uint_type i=0; i < __pts().size2(); ++i )
    {
        j0( i,i ) = ublas::row( Jac,0 )( i );
        j1( i,i ) = ublas::row( Jac,1 )( i );
        j2( i,i ) = ublas::row( Jac,2 )( i );
        j4( i,i ) = ublas::row( Jac,4 )( i );
        j5( i,i ) = ublas::row( Jac,5 )( i );
    }

#endif

    res[0] = ublas::prod( d1, j0 );
    res[1] = ublas::prod( d1, j1 ) + ublas::prod( d2, j4 );
    res[2] = ublas::prod( d1, j2 ) + ublas::prod( d2, j5 ) + d3;

    return res;
}


}
#endif /* __BoundAdapted_H */
