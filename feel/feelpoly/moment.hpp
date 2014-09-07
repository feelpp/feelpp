/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-12-05

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2011 Universite Joseph Fourier

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
   \file moment.hpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2005-12-05
 */
#ifndef __Moment_H
#define __Moment_H 1

#include <vector>

#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feelalg/lu.hpp>
#include <feel/feelpoly/expansions.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelpoly/polynomialset.hpp>


#include <feel/feelpoly/continuity.hpp>
#include <feel/feelpoly/discontinuous.hpp>

namespace Feel
{

/**
 * \class Moment
 * \brief Moment polynomial basis
 *
 * This class represents the Moment polynomials up to degree \c
 * Degree on a simplex in dimension \c Dim.
 *
 *
 * The moment polynomials in 1D are the canonical
 * basis to represent any polynomial set.
 *
 * \ingroup Polynomial
 * @author Gilles Steiner
 * @see
 */
template<uint16_type Dim,
         uint16_type Degree,
         //         typename NormalizationPolicy = Normalized<false>,
         typename TheConvex=Simplex<Dim>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
class Moment
{
public:

    static const uint16_type nDim = Dim;
    static const uint16_type nRealDim = Dim;
    static const uint16_type nOrder = Degree;
    static const uint16_type nConvexOrder = mpl::if_<mpl::bool_<TheConvex::is_simplex>,
                             mpl::int_<nDim+nOrder+1>,
                             mpl::int_<nOrder+2> >::type::value;
    //    static const bool is_normalized = NormalizationPolicy::is_normalized;
    static const uint16_type convex_is_simplex = TheConvex::is_simplex;
    static const uint16_type convex_is_hypercube = TheConvex::is_hypercube;
    static const bool is_product = false;

    /** @name Typedefs
     */
    //@{

    typedef Moment<Dim, Degree,TheConvex,T, StoragePolicy> self_type;

    typedef self_type basis_type;

    /*
     * numerical type
     */
    typedef T value_type;

    template<uint16_type order, typename V = value_type>
    struct Convex
    {
        typedef typename mpl::if_<mpl::bool_<TheConvex::is_simplex>,
                mpl::identity<Simplex<nDim, order, /*nRealDim*/nDim> >,
                mpl::identity<Hypercube<nDim, order, /*nRealDim*/nDim> > >::type::type type;
        typedef typename mpl::if_<mpl::bool_<TheConvex::is_simplex>,
                mpl::identity<Reference<Simplex<nDim, order, /*nRealDim*/nDim>, nDim, order, nDim/*nRealDim*/, V > >,
                mpl::identity<Reference<Hypercube<nDim, order, /*nRealDim*/nDim>, nDim, order, nDim/*nRealDim*/, V > > >::type::type reference_type;
    };

    /*
     * Geometry where the polynomials are defined and constructed
     */
    typedef TheConvex convex_type;
    typedef typename Convex<nOrder>::reference_type reference_convex_type;
    //typedef typename reference_convex_type::points_type points_type;

    /*
     * storage policy
     */
    typedef StoragePolicy<value_type> storage_policy;
    typedef typename storage_policy::matrix_type matrix_type;
    typedef typename storage_policy::vector_matrix_type vector_matrix_type;
    typedef typename storage_policy::matrix_node_type matrix_node_type;
    typedef typename storage_policy::node_type node_type;
    typedef typename storage_policy::points_type points_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Moment()
        :
        M_refconvex(),
        M_pts( M_refconvex.points() ),//M_refconvex.makePoints( nDim, 0 ) ),
        M_coord( convex_type::polyDims( nOrder ),nDim )
    {
        //std::cout << "[polydims: " << convex_type::polyDims( nOrder ) << "\n";
        if ( nDim == 1 )
            for ( uint32_type i=0; i < M_coord.size1(); ++i )
                M_coord( i,0 ) = i;

        else if ( nDim == 2 )
        {
            uint32_type G_i = 0;

            if ( convex_is_simplex )
            {
                for ( int32_type n = 0; n < nOrder+1; ++n )
                    for ( int32_type i = n; i >= 0; --i )
                    {
                        M_coord( G_i,0 ) = i;
                        M_coord( G_i,1 ) = n-i;
                        ++G_i;
                    }
            }

            else if ( convex_is_hypercube )
            {
                for ( int32_type n = 0; n < nOrder+1; ++n )
                    for ( int32_type i = 0; i < nOrder+1; ++i )
                    {
                        M_coord( G_i,0 ) = i;
                        M_coord( G_i,1 ) = n;
                        ++G_i;
                    }
            }
        }

        else if ( nDim == 3 )
        {
            uint32_type G_i = 0;

            if ( convex_is_simplex )
            {
                for ( int32_type n = 0; n < nOrder+1; ++n )
                    for ( int32_type i = n; i >= 0; --i )	   // x id
                        for ( int32_type j = n-i; j >= 0 ; --j )
                        {
                            M_coord( G_i,0 ) = i;
                            M_coord( G_i,1 ) = j;
                            M_coord( G_i,2 ) = n-i-j;
                            ++G_i;
                        }
            }

            else if ( convex_is_hypercube )
            {
                for ( int32_type n = 0; n < nOrder+1; ++n )
                    for ( int32_type i = 0; i < nOrder+1; ++i )
                        for ( int32_type j = 0; j < nOrder+1; ++j )
                        {
                            M_coord( G_i,0 ) = j;
                            M_coord( G_i,1 ) = i;
                            M_coord( G_i,1 ) = n;
                            ++G_i;
                        }
            }
        }

        //std::cout << "nDim: " << nDim << "\n";
        //std::cout << "orders: " << M_coord << "\n";
        //std::cout << "pts: " << M_pts << "\n";


        matrix_type A( ublas::trans( evaluate( M_pts ) ) );
        //std::cout << "A=" << A << "\n";
        matrix_type D = ublas::identity_matrix<value_type>( A.size1(), A.size2()  );
        //std::cout << "D=" << D << "\n";
        LU<matrix_type> lu( A );
        matrix_type C = lu.solve( D );
        //std::cout << "C=" << C << "\n";
        vector_matrix_type d ( derivate( M_pts ) );
        M_D.resize( d.size() );

        for ( size_type i = 0; i < d.size(); ++i )
        {
            M_D[i] = ublas::prod( C, ublas::trans( d[i] ) );
            //std::cout << "DM=" << M_D[i] << "\n";
        }
    }
    Moment( Moment const & d )
        :
        M_refconvex(),
        M_pts( d.M_pts ),
        M_coord( d.M_coord ),
        M_D( d.M_D )
    {}

    ~Moment()
    {}

    //@}

    /**
     * Access to the points of the reference convex associated
     **/

    points_type points()
    {
        return M_pts;
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

    //ublas::matrix_column<matrix_type const> operator()( node_type const& pt ) const
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

    ublas::matrix<int16_type> coord()
    {
        return M_coord;
    }

    /**
     * \return the maximum degree of the Moment polynomial to be
     * constructed
     */
    uint16_type degree() const
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
     * \return true if the Moment polynomials are normalized, false
     * otherwise
     */

    /**
     *  bool isNormalized() const { return is_normalized; }
     *
     **/

    /**
     * \return the \c familyName()
     */
    std::string familyName() const
    {
        return "moment";
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Moment polynomials is an orthonormal basis, the coefficients
     * of the polynomials of the basis are the canonical vectors and
     * represented by the identity matrix (lines are polynomials and
     * columns are the polynomial basis )
     *
     * This function is correct only if we use the Moment polynomials
     * as a basis
     */
    matrix_type coeff() const
    {
#if 0
        std::cout << "[Moment::coeff] coeff = "
                  << ublas::identity_matrix<value_type>( reference_convex_type::polyDims( nOrder ), M_pts.size2() )
                  << "\n";
#endif
        return ublas::identity_matrix<value_type>( reference_convex_type::polyDims( nOrder ), M_pts.size2() );
    }


    /**
     * evaluate the Moment polynomials at a set of points \p __pts
     *
     * \arg __x is a set of points
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
     * \brief derivatives of Moment polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    matrix_type const& d( uint16_type i ) const
    {
        return M_D[i];
    }

    template<template<uint16_type> class PolySetType = Scalar>
    Polynomial<self_type,PolySetType> pick( int i, int c = 0 ) const
    {
        size_type dim_p = convex_type::polyDims( nDim );
        int nComponents = PolySetType<Dim>::nComponents;
        matrix_type coeff( nComponents, dim_p );
        coeff = ublas::scalar_matrix<value_type>( coeff.size1(), coeff.size2(), 0. );
        coeff( c, i ) = 1;
        //std::cout << "pick.coeff(" << i << "," << c << ") = "  << coeff << "\n";
        return Polynomial<self_type, PolySetType>( *this, coeff, true );
        //std::cout << "p1.coeff(" << i << "," << c << ") = "  << p.coeff() << "\n";
        //return p;
    }

    /**
     * SP_integrate provide the integration of the \f$ \prod_i Comp(i,:)  \f$
     * where \f$ Comp(i,j) \f$ is the \f$j^e\f$ coefficient of the \f$i^e\f$
     * function in the moment Basis.
     * The S stands for "simplex product" and the dimension of the domain is
     * the same as the Basis space dimension.
     **/
    value_type SP_integrate( ublas::matrix<value_type> const& __Comp )
    {
        return SP_integrate( __Comp, mpl::int_<nDim>() );
    }


    /**
     * S_integrate provide the integration of the \f$ \prod_i Comp(i,:)  \f$
     * where \f$ Comp(i,j) \f$ is the \f$j^e\f$ coefficient of the \f$i^e\f$
     * function in the moment Basis.
     * The S stands for "simplex" and the argument Dim is the space dimension of the domain.
     **/
    value_type S_integrate( ublas::matrix<value_type> const& __Comp )
    {
        return S_integrate( __Comp, mpl::int_<nDim>() );
    }




    //@}

private:

    /**
     * Evaluation at a set of points of the expansion basis in 1D
     */
    matrix_type
    evaluate( points_type const& __pts, mpl::int_<1> ) const
    {
        matrix_type m ( nOrder+1 ,__pts.size2()  );
        ublas::vector<value_type> pts( ublas::row( __pts,0 ) );

        for ( uint32_type i = 0; i < m.size2(); ++i )
            m( 0,i ) = 1;

        for ( uint32_type i = 1; i < m.size1(); ++i )
            ublas::row( m,i ) = ublas::element_prod( ublas::row( m,i-1 ), pts );

        return m;
    }

    /**
     * Copy of evaluate 1D with vector type points
     **/

    matrix_type
    evaluate( ublas::vector<value_type> const& __pts ) const
    {
        matrix_type m ( nOrder+1 ,__pts.size()  );

        for ( uint32_type i = 0; i < m.size2(); ++i )
            m( 0,i ) = 1;

        for ( uint32_type i = 1; i < m.size1(); ++i )
            ublas::row( m,i ) = ublas::element_prod( ublas::row( m,i-1 ), __pts );

        return m;
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
        matrix_type E ( evaluate( __pts ) );

        for ( uint32_type i = 0; i < E.size2(); ++i )
            ublas::row( D[0], 0 )( i ) = value_type( 0.0 );

        for ( uint32_type i = 1; i < E.size1(); ++i )
            ublas::row( D[0], i ) = value_type( i )*ublas::row( E,i-1 );

        return D;
    }

    /**
     * derivation at a set of points of the expansion basis in 1D
     */
    vector_matrix_type
    derivate( ublas::vector<value_type> const& __pts ) const
    {
        vector_matrix_type D( 1 );

        D[0].resize( nOrder+1, __pts.size() );
        matrix_type E ( evaluate( __pts ) );

        for ( uint32_type i = 0; i < E.size2(); ++i )
            ublas::row( D[0], 0 )( i ) = value_type( 0.0 );

        for ( int32_type i = 1; i < E.size1(); ++i )
            ublas::row( D[0], i ) = value_type( i )*ublas::row( E,i-1 );

        return D;
    }

    /**
     * Evaluation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    matrix_type evaluate( points_type const& __pts, mpl::int_<2> ) const;

    /**
     * derivation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    template<typename AE>
    vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<2> ) const;

    /**
     * Evaluation at a set of points of the expansion basis in 3D on
     * the tetrahedron
     */
    matrix_type evaluate( points_type const& __pts, mpl::int_<3> ) const;

    /**
     * derivation at a set of points of the expansion basis in 3D on
     * the tetrahedron
     */
    template<typename AE>
    vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<3> ) const;

    /**
     * Exact integration on a simplex product in the moment basis in 1D.
     **/

    value_type SP_integrate( ublas::matrix<value_type> const& __Comp, mpl::int_<1> )
    {
        value_type res( 0.0 );

        if ( __Comp.size1() == 1 ) // Only 1 function to integrate : \f$ \int u \f$
        {
            for ( int ki = 0; ki <__Comp.size2(); ki += 2 )
            {
                res+= 2.0*__Comp( 0,ki ) / value_type( ki+1 );
            }
        }

        else if ( __Comp.size1() == 2 ) // 2 functions to integrate : \f$ \int uv \f$
        {
            for ( int k1 = 0; k1 <__Comp.size2(); ++k1 )
                for ( int k2 = 0; k2 <__Comp.size2(); ++k2 )
                {
                    if ( ( k1+k2 )%2 == 0 )
                        res+= 2.0 * __Comp( 0,k1 ) * __Comp( 1,k2 )  / value_type( k1+k2+1 );
                }
        }

        else if ( __Comp.size1() == 3 ) // 3 functions to integrate : \f$ \int uvw \f$
        {
            for ( int k1 = 0; k1 <__Comp.size2(); ++k1 )
                for ( int k2 = 0; k2 <__Comp.size2(); ++k2 )
                    for ( int k3 = 0; k3 <__Comp.size2(); ++k3 )
                    {
                        if ( ( k1+k2+k3 )%2 == 0 )
                            res+= 2.0 * __Comp( 0,k1 ) * __Comp( 1,k2 )* __Comp( 2,k3 )  / value_type( k1+k2+k3+1 );
                    }
        }

        return res;
    }


    /**
     * Exact integration on a simplex product in the moment basis in 2D.
     **/

    value_type SP_integrate( ublas::matrix<value_type> const& __Comp, mpl::int_<2> )
    {
        value_type res( 0.0 );

        if ( __Comp.size1() == 1 ) // Only 1 function to integrate : \f$ \int u \f$
        {
            for ( int ki = 0; ki <__Comp.size2(); ++ki )
            {
                int a( this->coord()( ki,0 ) );
                int b( this->coord()( ki,1 ) );

                if ( ( a%2 == 0 ) && ( b%2 == 0 ) )
                    res+= 4.0*__Comp( 0,ki ) / value_type( a+1 ) / value_type( b+1 );
            }
        }

        else if ( __Comp.size1() == 2 ) // 2 functions to integrate : \f$ \int uv \f$
        {
            for ( int k1 = 0; k1 <__Comp.size2(); ++k1 )
                for ( int k2 = 0; k2 <__Comp.size2(); ++k2 )
                {
                    int a( this->coord()( k1,0 ) + this->coord()( k2,0 ) );
                    int b( this->coord()( k2,1 ) + this->coord()( k2,1 ) );

                    if ( ( a%2 == 0 ) && ( b%2 == 0 ) )
                        res+= 4.0 * __Comp( 0,k1 ) * __Comp( 1,k2 )  / value_type( a+1 ) / value_type( b+1 );
                }
        }

        else if ( __Comp.size1() == 3 ) // 3 functions to integrate : \f$ \int uvw \f$
        {
            for ( int k1 = 0; k1 <__Comp.size2(); ++k1 )
                for ( int k2 = 0; k2 <__Comp.size2(); ++k2 )
                    for ( int k3 = 0; k3 <__Comp.size2(); ++k3 )
                    {
                        int a( this->coord()( k1,0 ) + this->coord()( k2,0 ) + this->coord()( k3,0 ) );
                        int b( this->coord()( k1,1 ) + this->coord()( k2,1 ) + this->coord()( k3,1 ) );

                        if ( ( a%2 == 0 ) && ( b%2 == 0 ) )
                            res+= 4.0 * __Comp( 0,k1 ) * __Comp( 1,k2 ) * __Comp( 2,k3 ) / value_type( a+1 ) / value_type( b+1 );
                    }
        }

        return res;
    }
    /**
     * Exact integration on a simplex product in the moment basis in 3D.
     **/

    value_type SP_integrate( ublas::matrix<value_type> const& __Comp, mpl::int_<3> )
    {
        value_type res( 0.0 );

        if ( __Comp.size1() == 1 ) // Only 1 function to integrate : \f$ \int u \f$
        {
            for ( int ki = 0; ki <__Comp.size2(); ++ki )
            {
                int a( this->coord()( ki,0 ) );
                int b( this->coord()( ki,1 ) );
                int c( this->coord()( ki,2 ) );

                if ( ( a%2 == 0 ) && ( b%2 == 0 ) && ( c%2 == 0 ) )
                    res+= 8.0*__Comp( 0,ki ) / value_type( a+1 ) / value_type( b+1 ) / value_type( c+1 );
            }
        }

        else if ( __Comp.size1() == 2 ) // 2 functions to integrate : \f$ \int uv \f$
        {
            for ( int k1 = 0; k1 <__Comp.size2(); ++k1 )
                for ( int k2 = 0; k2 <__Comp.size2(); ++k2 )
                {
                    int a( this->coord()( k1,0 ) + this->coord()( k2,0 ) );
                    int b( this->coord()( k1,1 ) + this->coord()( k2,1 ) );
                    int c( this->coord()( k1,2 ) + this->coord()( k2,2 ) );

                    if ( ( a%2 == 0 ) && ( b%2 == 0 ) && ( c%2 == 0 ) )
                        res+= 8.0*__Comp( 0,k1 )* __Comp( 1,k2 ) / value_type( a+1 ) / value_type( b+1 ) / value_type( c+1 );
                }
        }

        else if ( __Comp.size1() == 3 ) // 3 functions to integrate : \f$ \int uvw \f$
        {
            for ( int k1 = 0; k1 <__Comp.size2(); ++k1 )
                for ( int k2 = 0; k2 <__Comp.size2(); ++k2 )
                    for ( int k3 = 0; k3 <__Comp.size2(); ++k3 )
                    {
                        int a( this->coord()( k1,0 ) + this->coord()( k2,0 ) + this->coord()( k3,0 ) );
                        int b( this->coord()( k1,1 ) + this->coord()( k2,1 ) + this->coord()( k3,1 ) );
                        int c( this->coord()( k1,2 ) + this->coord()( k2,2 ) + this->coord()( k3,2 ) );

                        if ( ( a%2 == 0 ) && ( b%2 == 0 ) && ( c%2 == 0 ) )
                            res+= 8.0*__Comp( 0,k1 )* __Comp( 1,k2 )* __Comp( 2,k3 )/value_type( a+1 )/value_type( b+1 )/value_type( c+1 );
                    }
        }

        return res;
    }

    /**
     * Exact integration on a 1D-simplex.
     **/

    value_type S_integrate( ublas::matrix<value_type> const& __Comp, mpl::int_<1> )
    {
        return SP_integrate(  __Comp, mpl::int_<nDim>() );
    }

    /**
     * Exact integration on a 2D-simplex.
     **/

    value_type S_integrate( ublas::matrix<value_type> const& __Comp, mpl::int_<2> )
    {
        value_type res( 0.0 );

        if ( __Comp.size1() == 1 ) // Only 1 function to integrate : \f$ \int u \f$
        {
            for ( uint32_type ki = 0; ki <__Comp.size2(); ++ki )
            {
                int i( this->coord()( ki,0 ) );
                int j( this->coord()( ki,1 ) );

                // i and j parity

                int pi = i%2;
                int pj = j%2;

                if ( ( i%2 == 0 ) || ( ( i+j )%2 != 0 ) )
                {
                    if ( pi==0 && pj==0 ) // (-1)/(j+1)*(-2/(i+1))
                        res+= __Comp( 0,ki )/value_type( j+1 ) * ( 2.0/value_type( i+1 ) ) ;

                    else if ( pi==0 && pj==1 ) // 1/(j+1)*(2/(i+j+1)-2/(i+1))
                        res+= __Comp( 0,ki )/value_type( j+1 ) * ( 2.0/value_type( i+j+2 )-2.0/value_type( i+1 ) ) ;

                    else if ( pi==1 && pj==0 ) // (-1)/(j+1)*(2/(i+j+1))
                        res-= __Comp( 0,ki )/value_type( j+1 ) * ( 2.0/value_type( i+j+2 ) ) ;
                }
            }
        }

        else if ( __Comp.size1() == 2 ) // 2 function2 to integrate : \f$ \int uv \f$
        {
            for ( uint32_type k1 = 0; k1 <__Comp.size2(); ++k1 )
                for ( uint32_type k2 = 0; k2 <__Comp.size2(); ++k2 )
                {
                    int i( this->coord()( k1,0 ) + this->coord()( k2,0 ) );
                    int j( this->coord()( k1,1 ) + this->coord()( k2,1 ) );

                    // i and j parity

                    int pi = i%2;
                    int pj = j%2;

                    if ( ( i%2 == 0 ) || ( ( i+j )%2 != 0 ) )
                    {
                        if ( pi==0 && pj==0 ) // (-1)/(j+1)*(-2/(i+1))
                            res+= __Comp( 0,k1 ) * __Comp( 1,k2 )/value_type( j+1 ) * ( 2.0/value_type( i+1 ) ) ;

                        else if ( pi==0 && pj==1 ) // 1/(j+1)*(2/(i+j+1)-2/(i+1))
                            res+= __Comp( 0,k1 ) * __Comp( 1,k2 )/value_type( j+1 ) * ( 2.0/value_type( i+j+2 )-2.0/value_type( i+1 ) ) ;

                        else if ( pi==1 && pj==0 ) // (-1)/(j+1)*(2/(i+j+1))
                            res-= __Comp( 0,k1 ) * __Comp( 1,k2 )/value_type( j+1 ) * ( 2.0/value_type( i+j+2 ) ) ;
                    }
                }
        }

        else if ( __Comp.size1() == 3 ) // 3 functions to integrate : \f$ \int uvw \f$
        {
            for ( uint32_type k1 = 0; k1 <__Comp.size2(); ++k1 )
                for ( uint32_type k2 = 0; k2 <__Comp.size2(); ++k2 )
                    for ( uint32_type k3 = 0; k3 <__Comp.size2(); ++k3 )
                    {
                        int i( this->coord()( k1,0 ) + this->coord()( k2,0 )+ this->coord()( k3,0 ) );
                        int j( this->coord()( k1,1 ) + this->coord()( k2,1 )+ this->coord()( k3,1 ) );

                        // i and j parity

                        int pi = i%2;
                        int pj = j%2;

                        if ( ( i%2 == 0 ) || ( ( i+j )%2 != 0 ) )
                        {
                            if ( pi==0 && pj==0 ) // (-1)/(j+1)*(-2/(i+1))
                                res+= __Comp( 0,k1 ) * __Comp( 1,k2 ) * __Comp( 2,k3 )/value_type( j+1 ) * ( 2.0/value_type( i+1 ) ) ;

                            else if ( pi==0 && pj==1 ) // 1/(j+1)*(2/(i+j+1)-2/(i+1))
                                res+= __Comp( 0,k1 ) * __Comp( 1,k2 ) * __Comp( 2,k3 )/value_type( j+1 ) * ( 2.0/value_type( i+j+2 )-2.0/value_type( i+1 ) ) ;

                            else if ( pi==1 && pj==0 ) // (-1)/(j+1)*(2/(i+j+1))
                                res-= __Comp( 0,k1 ) * __Comp( 1,k2 ) * __Comp( 2,k3 )/value_type( j+1 ) * ( 2.0/value_type( i+j+2 ) ) ;
                        }
                    }
        }

        return res;
    }

    /**
     * Exact integration on a 3D-simplex.
     **/

    value_type S_integrate( ublas::matrix<value_type> const& __Comp, mpl::int_<3> )
    {
        return 0.0;
    }



private:
    reference_convex_type M_refconvex;
    points_type M_pts;

    /**
     * In this matrix I will save the degree \f$ i,j,k \f$ in each direction
     * for the moment basis polynomials of type \f$ x^i y^j z^k \f$.
     */

    ublas::matrix<int16_type> M_coord;

    /**
     * Derivation matrix
     */
    std::vector<matrix_type> M_D;
};

template<uint16_type Dim,
         uint16_type Degree,
         typename TheConvex,
         typename T,
         template<class> class StoragePolicy>
typename Moment<Dim, Degree,TheConvex, T, StoragePolicy>::matrix_type
Moment<Dim, Degree,TheConvex, T, StoragePolicy>::evaluate( points_type const& __pts, mpl::int_<2> ) const
{
    matrix_type res( convex_type::polyDims( nOrder ), __pts.size2() );

    ublas::vector<value_type> pts_x = ublas::row( __pts, 0 );
    ublas::vector<value_type> pts_y = ublas::row( __pts, 1 );

    matrix_type E_x( evaluate( pts_x ) );
    matrix_type E_y( evaluate( pts_y ) );

    //  std::cout << "\nDimensions E_x = " << E_x.size1() << " x " << E_x.size2() << std::endl;
    //  std::cout << "\nDimensions E_y = " << E_y.size1() << " x " << E_y.size2() << std::endl;
    //  std::cout << "Order = " << nOrder+1 << std::endl;

    uint32_type G_i=0; // Global indice

    if ( convex_is_hypercube )
    {
        for ( int32_type n = 0; n < nOrder+1; ++n )
            for ( int32_type i = 0; i < nOrder+1; ++i )

            {
                ublas::row( res, G_i )=ublas::element_prod( ublas::row( E_x, i ) , ublas::row( E_y, n ) );
                ++G_i;
            }
    }
    else if ( convex_is_simplex )
    {
        for ( int32_type n = 0; n < nOrder+1; ++n )
            for ( int32_type i = n; i >= 0; --i )
            {
                ublas::row( res, G_i )=ublas::element_prod( ublas::row( E_x, i ) , ublas::row( E_y, n-i ) );
                ++G_i;
            }
    }

    return res;
}

template<uint16_type Dim,
         uint16_type Degree,
         typename TheConvex,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename Moment<Dim, Degree, TheConvex,  T, StoragePolicy>::vector_matrix_type
Moment<Dim, Degree, TheConvex,  T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<2> ) const
{
    vector_matrix_type res( 2 );
    res[0].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[1].resize( convex_type::polyDims( nOrder ), __pts().size2() );

    ublas::vector<value_type> pts_x = ublas::row( __pts(), 0 );
    ublas::vector<value_type> pts_y = ublas::row( __pts(), 1 );

    matrix_type E_x( evaluate( pts_x ) );
    matrix_type E_y( evaluate( pts_y ) );

    vector_matrix_type D_x( derivate( pts_x ) );
    vector_matrix_type D_y( derivate( pts_y ) );

    uint32_type G_i=0;

    if ( convex_is_simplex )
    {
        for ( int32_type n = 0; n < nOrder+1; ++n )
            for ( int32_type i = n; i >= 0; --i )
            {
                ublas::row( res[0], G_i )=ublas::element_prod( ublas::row( D_x[0], i ) , ublas::row( E_y, n-i ) );
                ublas::row( res[1], G_i )=ublas::element_prod( ublas::row( E_x, i ) , ublas::row( D_y[0], n-i ) );
                ++G_i;
            }
    }
    else if ( convex_is_hypercube )
    {
        for ( int32_type n = 0; n < nOrder+1; ++n )
            for ( int32_type i = 0; i < nOrder+1; ++i )
            {
                ublas::row( res[0], G_i )=ublas::element_prod( ublas::row( D_x[0], i ) , ublas::row( E_y, n ) );
                ublas::row( res[1], G_i )=ublas::element_prod( ublas::row( E_x, i ) , ublas::row( D_y[0], n ) );
                ++G_i;
            }
    }

    return res;
}

template<uint16_type Dim,
         uint16_type Degree,
         typename TheConvex,
         typename T,
         template<class> class StoragePolicy>
typename Moment<Dim, Degree, TheConvex,  T, StoragePolicy>::matrix_type
Moment<Dim, Degree, TheConvex, T, StoragePolicy>::evaluate( points_type const& __pts, mpl::int_<3> ) const
{
    matrix_type res( convex_type::polyDims( nOrder ), __pts.size2() );

    FEELPP_ASSERT( __pts.size1() == 3 )( __pts.size1() ).error( "invalid space dimension" );

    ublas::vector<value_type> pts_x = ublas::row( __pts, 0 );
    ublas::vector<value_type> pts_y = ublas::row( __pts, 1 );
    ublas::vector<value_type> pts_z = ublas::row( __pts, 2 );

    matrix_type E_x( evaluate( pts_x ) );
    matrix_type E_y( evaluate( pts_y ) );
    matrix_type E_z( evaluate( pts_z ) );

    uint32_type G_i=0; // Global indice

    for ( int32_type n = 0; n < nOrder+1; ++n )
        for ( int32_type i = n; i >= 0; --i )	   // x id
            for ( int32_type j = n-i; j >= 0 ; --j )
            {
                ublas::vector<value_type> tmp( ublas::element_prod( ublas::row( E_x, i ) , ublas::row( E_y, j ) ) );
                ublas::row( res, G_i )=ublas::element_prod( tmp , ublas::row( E_z, n-i-j ) );
                ++G_i;
            }

    return res;
}

template<uint16_type Dim,
         uint16_type Degree,
         typename TheConvex,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename Moment<Dim, Degree, TheConvex,  T, StoragePolicy>::vector_matrix_type
Moment<Dim, Degree, TheConvex,  T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<3> ) const
{
    vector_matrix_type res( 3 );
    res[0].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[1].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[2].resize( convex_type::polyDims( nOrder ), __pts().size2() );

    FEELPP_ASSERT( __pts().size1() == 3 )( __pts().size1() ).error( "invalid space dimension" );

    ublas::vector<value_type> pts_x = ublas::row( __pts(), 0 );
    ublas::vector<value_type> pts_y = ublas::row( __pts(), 1 );
    ublas::vector<value_type> pts_z = ublas::row( __pts(), 2 );

    matrix_type E_x( evaluate( pts_x ) );
    matrix_type E_y( evaluate( pts_y ) );
    matrix_type E_z( evaluate( pts_z ) );

    vector_matrix_type D_x( derivate( pts_x ) );
    vector_matrix_type D_y( derivate( pts_y ) );
    vector_matrix_type D_z( derivate( pts_z ) );

    uint32_type G_i=0;

    for ( int32_type n = 0; n < nOrder+1; ++n )
        for ( int32_type i = n; i >= 0; --i )	   // x id
            for ( int32_type j = n-i; j >=0 ; --j )
            {
                ublas::vector<value_type> tmp( ublas::element_prod( ublas::row( D_x[0], i ) , ublas::row( E_y, j ) ) );
                ublas::row( res[0], G_i )=ublas::element_prod( tmp , ublas::row( E_z, n-i-j ) );

                tmp = ublas::element_prod( ublas::row( E_x, i ) , ublas::row( D_y[0], j ) );
                ublas::row( res[1], G_i )=ublas::element_prod( tmp , ublas::row( E_z, n-i-j ) );

                tmp = ublas::element_prod( ublas::row( E_x, i ) , ublas::row( E_y, j ) );
                ublas::row( res[2], G_i )=ublas::element_prod( tmp , ublas::row( D_z[0], n-i-j ) );

                ++G_i;
            }

    return res;
}

namespace fem
{
/// \cond DETAIL
namespace detail
{
/**
 * \internal
 * \class MomentPolynomialSet
 * \brief a set of orthonormal polynomials over a convex
 *
 * On the simplicies we use the Dubiner basis
 *
 */
template<uint16_type Dim,
         uint16_type Order,
         uint16_type RealDim,
         template<uint16_type> class PolySetType = Scalar,
         typename T = double,
         template<uint16_type,uint16_type,uint16_type> class Convex = Simplex>
class MomentPolynomialSet
    :
public PolynomialSet<Moment<Dim, Order, Convex<Dim,1,Dim>, T, StorageUBlas>, PolySetType >
{
    typedef PolynomialSet<Moment<Dim, Order, Convex<Dim,1,Dim>, T, StorageUBlas>, PolySetType > super;
public:

    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;
    static const uint16_type nRealDim = RealDim;
    static const bool isTransformationEquivalent = true;
    typedef MomentPolynomialSet<Dim, Order,RealDim, PolySetType, T, Simplex> self_type;
    typedef self_type component_basis_type;

    typedef typename super::polyset_type polyset_type;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const bool is_continuous = false;
    static const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const bool is_product = true;
    static const bool isContinuous = false;
    typedef Discontinuous continuity_type;

    typedef typename super::component_type component_type;

    typedef T value_type;
    typedef Moment<Dim, Order, Convex<Dim,1,Dim>, T, StorageUBlas> basis_type;
    typedef Convex<Dim, Order, /*RealDim*/Dim> convex_type;
    template<int O>
    struct convex
    {
        typedef Convex<Dim, O, /*RealDim*/Dim> type;
    };
    typedef Reference<convex_type, nDim, nOrder, nDim/*nRealDim*/, value_type> reference_convex_type;

    typedef typename super::polynomial_type polynomial_type;

    //!< Number of degrees of freedom per vertex
    static const uint16_type nDofPerVertex = convex_type::nbPtsPerVertex;
    //!< Number of degrees  of freedom per edge
    static const uint16_type nDofPerEdge = convex_type::nbPtsPerEdge;
    //!< Number of degrees  of freedom per face
    static const uint16_type nDofPerFace = convex_type::nbPtsPerFace;

    //!< Number of degrees  of freedom per volume
    static const uint16_type nDofPerVolume = convex_type::nbPtsPerVolume;

    static const uint16_type nLocalDof = convex_type::numPoints;

    static const uint16_type nDof = nLocalDof;
    static const uint16_type nNodes = nDof;
    static const uint16_type nDofGrad = super::nDim*nDof;
    static const uint16_type nDofHess = super::nDim*super::nDim*nDof;
    //typedef typename matrix_node<value_type>::type points_type;
    typedef StorageUBlas<value_type> storage_policy;
    typedef typename storage_policy::matrix_type matrix_type;
    typedef typename storage_policy::vector_matrix_type vector_matrix_type;
    typedef typename storage_policy::matrix_node_type matrix_node_type;
    typedef typename storage_policy::points_type points_type;
    typedef typename storage_policy::node_type node_type;
    MomentPolynomialSet()
        :
        super( basis_type() )

    {
        ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );

        if ( is_tensor2 )
            std::cout << "[orthonormalpolynomialset] m = " << m << "\n";

        if ( !( ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
                                       m ) < 1e-10 ) )
            std::cout << "m1=" << m << "\n"
                      << "m2=" << polyset_type::toMatrix( polyset_type::toType( m ) ) << "\n"
                      << ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) - m ) << "\n";

        FEELPP_ASSERT( ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
                                              m ) < 1e-10 )( m ).warn ( "invalid transformation" );
        this->setCoefficient( polyset_type::toType( m ), true );
    }

    MomentPolynomialSet<Dim, Order, RealDim, Scalar,T, Simplex > toScalar() const
    {
        return MomentPolynomialSet<Dim, Order, RealDim, Scalar,T, Simplex >();
    }

    /**
     * \return the family name of the polynomial set
     */
    std::string familyName() const
    {
        return "dubiner";
    }


    points_type points() const
    {
        return points_type();
    }
    points_type points( int f ) const
    {
        return points_type();
    }
};

template<uint16_type Dim,
         uint16_type Order,
         uint16_type RealDim,
         template<uint16_type> class PolySetType,
         typename T,
         template<uint16_type,uint16_type,uint16_type> class Convex>
const uint16_type MomentPolynomialSet<Dim, Order, RealDim, PolySetType,T, Convex>::nLocalDof;

} // detail
/// \encond

template<uint16_type Order,
         template<uint16_type Dim> class PolySetType = Scalar>
class MomentPolynomialSet
{
public:
    template<uint16_type N,
             typename T = double,
             typename Convex = Simplex<N> >
    struct apply
    {
        typedef typename mpl::if_<mpl::bool_<Convex::is_simplex>,
                mpl::identity<detail::MomentPolynomialSet<N,Order,N,PolySetType,T,Simplex> >,
                mpl::identity<detail::MomentPolynomialSet<N,Order,N,PolySetType,T,Hypercube> > >::type::type result_type;
        typedef result_type type;
    };

    typedef MomentPolynomialSet<Order,Scalar> component_basis_type;
};
}


}
#endif /* __Moment_H */
