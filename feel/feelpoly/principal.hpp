/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
       Date: 2005-12-14

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
   \file principal.hpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2005-12-14
 */
#ifndef __principal_H
#define __principal_H 1

#include <vector>


#include <feel/feelalg/lu.hpp>
#include <feel/feelpoly/expansions.hpp>
#include <feel/feelpoly/policy.hpp>

namespace Feel
{

/**
 * \class Principal
 * \brief Principal modified functions
 *
 * This class is a useful class to construct multidimensionnal
 * boundary adapted expansions basis.
 * It is constructed following the section 3.2.3.1 of
 * the book from Sherwin and Karniadakis "Spectral/hp
 * element methods for computational fluid dynamics".
 *
 * \ingroup Polynomial
 * @author Gilles Steiner
 * @see dubiner.hpp
 */

template< uint16_type Degree, typename T = double,
          template<class> class StoragePolicy = StorageUBlas>
class Principal
{
public:

    static const uint16_type nOrder = Degree;

    /** @name Typedefs
     */
    //@{

    typedef Principal<Degree, T, StoragePolicy> self_type;

    /*
     * numerical type
     */
    typedef T value_type;

    /*
     * storage policy
     */
    typedef StoragePolicy<value_type> storage_policy;
    typedef typename storage_policy::matrix_type matrix_type;
    typedef typename storage_policy::vector_matrix_type vector_matrix_type;
    typedef typename storage_policy::matrix_node_type matrix_node_type;
    typedef typename storage_policy::node_type node_type;
    typedef typename storage_policy::vector_vector_matrix_type  vector_vector_matrix_type;
    typedef typename storage_policy::vector_type vector_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Principal()
        : M_a( value_type( 1.0 ) ), M_b( value_type( 1.0 ) )
    {}

    Principal( value_type a,value_type b )
        : M_a( a ), M_b( b )
    {}


    ~Principal()
    {}

    //@}


    self_type const& operator=( self_type const& d )
    {
        if ( this != &d )
        {
            M_a = d.M_a;
            M_b = d.M_b;
        }

        return *this;
    }

    /** @name Accessors
     */
    //@{

    uint16_type degree() const
    {
        return nOrder;
    }

    //@}

    /** @name  Methods
     */
    //@{

    value_type a() const
    {
        return M_a;
    }
    value_type b() const
    {
        return M_b;
    }

    /**
     * evaluate the Principal functions at a set of points \p __pts
     *
     * \arg __pts is a set of points in one dimension
     */

    /** One order principal function \f$ \psi_i(z) \f$ **/

    matrix_type evaluate_1( vector_type const& __pts ) const;

    /** Second order principal function \f$ \psi_{ij}(z) \f$ **/

    vector_matrix_type evaluate_2( vector_type const& __pts ) const;


    /** Third order principal function \f$ \psi_{ijk}(z) \f$ **/

    vector_vector_matrix_type evaluate_3( vector_type const& __pts ) const;
    //@}


    matrix_type derivate_1(  vector_type const& __pts ) const;

    vector_matrix_type derivate_2( vector_type  const& __pts ) const;

    vector_vector_matrix_type derivate_3( vector_type  const& __pts ) const;

private:

    value_type M_a;
    value_type M_b;


};

template<uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename Principal<Degree, T, StoragePolicy>::matrix_type
Principal<Degree,  T, StoragePolicy>::evaluate_1( vector_type const& __pts ) const
{
    //  std::cout <<"[principal]1D evaluation ..."<< std::endl;
    matrix_type J ( JacobiBatchEvaluation<nOrder-1,value_type>( M_a, M_b, __pts ) );
    matrix_type D( nOrder+1,__pts.size() );

    vector_type ones( ublas::scalar_vector<value_type>( __pts.size(), value_type( 1.0 ) ) );

    ublas::row( D,0 ) = value_type( 0.5 )*( ones - __pts );
    ublas::row( D,nOrder ) = value_type( 0.5 )*( ones + __pts );

    vector_type tmp( ublas::element_prod( ublas::row( D,0 ),ublas::row( D,D.size1()-1 ) ) );

    for ( uint16_type i = 1; i < D.size1()-1; ++i )
        ublas::row( D, i ) = ublas::element_prod( tmp, ublas::row( J,i-1 ) );

    return D;
}


template<uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename Principal<Degree,  T, StoragePolicy>::vector_matrix_type
Principal<Degree,  T, StoragePolicy>::evaluate_2( vector_type const& __pts ) const
{
    //  std::cout <<"[principal]2D evaluation ..."<< std::endl;
    matrix_type psi_1( evaluate_1( __pts ) );

    vector_matrix_type m( nOrder+1 );
    vector_matrix_type J( nOrder-1 );

    m[0].resize( nOrder+1 ,__pts.size() );
    m[0] = psi_1;
    m[nOrder] = psi_1;

    vector_type ones( ublas::scalar_vector<value_type>( __pts.size(), value_type( 1.0 ) ) );

    vector_type tmp1 = ( ones - __pts )/value_type( 2.0 ); // (1-x)/2
    vector_type tmp2 = ( ones + __pts )/value_type( 2.0 ); // (1+x)/2
    vector_type tmp_i = tmp1;

    for ( uint16_type i= 1; i <= J.size() ; ++i )
    {
        m[i].resize( nOrder ,__pts.size() );
        J[i-1] = JacobiBatchEvaluation<nOrder-1,value_type>( value_type( 2.0*i+1.0 ), M_b, __pts );

        tmp_i=ublas::element_prod( tmp_i,tmp1 );

        ublas::row( m[i],0 ) = tmp_i;
        vector_type tmp3( ublas::element_prod( tmp2,tmp_i ) );

        for ( int16_type j=1; j < nOrder; ++j )
        {
            ublas::row( m[i],j ) = ublas::element_prod( tmp3, ublas::row( J[i-1],j-1 ) );
        }
    }

    return m;
}


template<uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename Principal<Degree,  T, StoragePolicy>::vector_vector_matrix_type
Principal<Degree,  T, StoragePolicy>::evaluate_3( vector_type const& __pts ) const
{
    //  std::cout <<"[principal]3D evaluation ..."<< std::endl;
    vector_vector_matrix_type m( nOrder+1 );

    vector_matrix_type psi_2( evaluate_2( __pts ) );

    // i=0, 0 <= j <= nOrder , 0 <= k <= nOrder
    m[0].resize( nOrder+1 );
    m[0] = psi_2;

    // i=nOrder, 0 <= j <= nOrder , 0 <= k <= nOrder
    m[nOrder].resize( nOrder+1 );
    m[nOrder] = psi_2;

    vector_type ones( ublas::scalar_vector<value_type>( __pts.size(), value_type( 1.0 ) ) );
    vector_type tmp1 = ( ones - __pts )/value_type( 2.0 ); // (1-x)/2
    vector_type tmp2 = ( ones + __pts )/value_type( 2.0 ); // (1+x)/2
    vector_type tmp_i = tmp1;

    vector_vector_matrix_type J( nOrder+1 );

    for ( int16_type i=1; i < nOrder; ++i )
    {
        m[i].resize( nOrder+1 );
        J[i].resize( nOrder+1 );

        // 1 <= i <= nOrder-1,  j = 0 , 0 <= k <= nOrder
        m[i][0].resize( psi_2[i].size1(), psi_2[i].size2() );
        m[i][0] = psi_2[i];

        // 1 <= i <= nOrder-1,  j = nOrder , 0 <= k <= nOrder
        m[i][nOrder].resize( psi_2[i].size1(),psi_2[i].size2() );
        m[i][nOrder] = psi_2[i];

        tmp_i=ublas::element_prod( tmp_i,tmp1 );

        vector_type tmp_i_j( tmp_i );

        for ( int16_type j=1; j < nOrder; ++j )
        {
            // 1 <= i <= nOrder-1,  1 <= j <= nOrder-1 , 0 <= k <= nOrder-1
            m[i][j].resize( nOrder,__pts.size() );
            J[i][j] = JacobiBatchEvaluation<nOrder-1,value_type>( ( 2.0*i+2.0*j+1.0 ), 1.0, __pts );
            tmp_i_j = ublas::element_prod( tmp_i_j,tmp1 );

            ublas::row( m[i][j],0 ) = tmp_i_j; // k=0

            for ( int16_type k=1; k < nOrder; ++k ) // 1 <= k <= nOrder-1
            {
                vector_type tmp3( ublas::element_prod( tmp2,tmp_i_j ) );
                ublas::row( m[i][j],k ) = ublas::element_prod( tmp3, ublas::row( J[i][j],k-1 ) );
            }
        }
    }

    return m;
}


template<uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename Principal<Degree,  T, StoragePolicy>::matrix_type
Principal<Degree,  T, StoragePolicy>::derivate_1( vector_type const& __pts ) const
{
    //  std::cout <<"[principal]1D derivation ..."<< std::endl;
    matrix_type D(  nOrder+1, __pts.size() );
    matrix_type J ( JacobiBatchDerivation<nOrder-1,value_type>( M_a, M_b, __pts ) );
    vector_type demi( ublas::scalar_vector<value_type>( __pts.size(), value_type( 0.5 ) ) );

    ublas::row( D,0 ) = -demi;
    ublas::row( D,nOrder ) = demi;

    vector_type ones( ublas::scalar_vector<value_type>( __pts.size(), value_type( 1.0 ) ) );


    matrix_type E ( JacobiBatchEvaluation<nOrder-1,value_type>( M_a, M_b, __pts ) );

    vector_type tmp( ublas::element_prod( 0.5*( ones - __pts ),0.5*( ones + __pts ) ) );

    for ( uint16_type i = 1; i < D.size1()-1; ++i )
        ublas::row( D, i ) = ublas::element_prod( tmp, ublas::row( J,i-1 ) )  - 0.5*ublas::element_prod( __pts, ublas::row( E,i-1 ) ); // \f[ = \frac{1-x}{2}\frac{1+x}{2} \frac{d}{dx}P_{i-1}^{(1,1)}(x) - \frac{x}{2}P_{i-1}^{(1,1)}(x) \f]

    return D;
}

template<uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename Principal<Degree,  T, StoragePolicy>::vector_matrix_type
Principal<Degree,  T, StoragePolicy>::derivate_2( vector_type const& __pts ) const
{
    //  std::cout <<"[principal]2D derivation ..."<< std::endl;
    matrix_type dpsi_1( derivate_1( __pts ) );

    vector_matrix_type m( nOrder+1 );
    vector_matrix_type J( nOrder-1 );
    vector_matrix_type dJ( nOrder-1 );

    m[0].resize( nOrder+1 ,__pts.size() );
    m[0] = dpsi_1;
    m[nOrder] = dpsi_1;

    vector_type ones( ublas::scalar_vector<value_type>( __pts.size(), value_type( 1.0 ) ) );

    vector_type tmp1 = ( ones - __pts )/value_type( 2.0 ); // (1-x)/2
    vector_type tmp2 = ( ones + __pts )/value_type( 2.0 ); // (1+x)/2

    vector_type tmp_prod = ublas::element_prod( tmp1,tmp2 );

    vector_type tmp_i = tmp1;

    for ( uint16_type i= 1; i <= J.size() ; ++i )
    {
        m[i].resize( nOrder ,__pts.size() );
        J[i-1]  = JacobiBatchEvaluation<nOrder-1,value_type>( ( 2.0*i+1.0 ), M_b, __pts );
        dJ[i-1] = JacobiBatchDerivation<nOrder-1,value_type>( ( 2.0*i+1.0 ), M_b, __pts );

        ublas::row( m[i],0 ) = ( - value_type( i ) - 1.0 ) *  tmp_i / 2.0;

        for ( int16_type j=1; j < nOrder; ++j )
        {
#if 0
            vector_type tmp3(
                ublas::element_prod( tmp1, ublas::row( J[i-1],j-1 ) )
                + 2.0*element_prod( tmp_prod, ublas::row( dJ[i-1],j-1 ) )
                - ( value_type( i ) + 1.0 ) * ublas::element_prod( tmp2,ublas::row( J[i-1],j-1 ) )
            );
            ublas::row( m[i],j ) = ( ublas::element_prod( tmp_i, tmp3 ) ) / value_type( 2.0 ) ;
#else
            vector_type A( ublas::element_prod( ublas::row( m[i],0 ), tmp2 ) );
            A = ublas::element_prod( A, ublas::row( J[i-1],j-1 ) );

            vector_type B( ublas::element_prod( 0.5*tmp_i , ublas::row( J[i-1],j-1 ) ) );

            vector_type C( ublas::element_prod( tmp_i , tmp_prod ) );
            C = ublas::element_prod( C ,ublas::row( dJ[i-1],j-1 ) );

            ublas::row( m[i],j ) = A + B + C;
#endif
        }

        tmp_i=ublas::element_prod( tmp_i,tmp1 );
    }

    return m;
}


template<uint16_type Degree,
         typename T,
         template<class> class StoragePolicy>
typename Principal<Degree,  T, StoragePolicy>::vector_vector_matrix_type
Principal<Degree,  T, StoragePolicy>::derivate_3( vector_type const& __pts ) const
{
    //  std::cout <<"[principal]3D derivation ..."<< std::endl;
    vector_vector_matrix_type m( nOrder+1 );

    vector_matrix_type dpsi_2( derivate_2( __pts ) );

    // i=0, 0 <= j <= nOrder , 0 <= k <= nOrder
    m[0].resize( nOrder+1 );
    m[0] = dpsi_2;

    // i=nOrder, 0 <= j <= nOrder , 0 <= k <= nOrder
    m[nOrder].resize( nOrder+1 );
    m[nOrder] = dpsi_2;

    vector_type ones( ublas::scalar_vector<value_type>( __pts.size(), value_type( 1.0 ) ) );
    vector_type tmp1 = ( ones - __pts )/value_type( 2.0 ); // (1-x)/2
    vector_type tmp2 = ( ones + __pts )/value_type( 2.0 ); // (1+x)/2
    vector_type tmp_i = ones;
    vector_type tmp_prod = ublas::element_prod( tmp1,tmp2 );

    vector_vector_matrix_type J( nOrder+1 );
    vector_vector_matrix_type dJ( nOrder+1 );

    for ( int16_type i=1; i < nOrder; ++i )
    {
        m[i].resize( nOrder+1 );
        J[i].resize( nOrder+1 );
        dJ[i].resize( nOrder+1 );

        // 1 <= i <= nOrder-1,  j = 0 , 0 <= k <= nOrder
        m[i][0].resize( dpsi_2[i].size1(), dpsi_2[i].size2() );
        m[i][0] = dpsi_2[i];

        // 1 <= i <= nOrder-1,  j = nOrder , 0 <= k <= nOrder
        m[i][nOrder].resize( dpsi_2[i].size1(),dpsi_2[i].size2() );
        m[i][nOrder] = dpsi_2[i];

        tmp_i=ublas::element_prod( tmp_i,tmp1 ); // [(1-x) / 2]^i

        vector_type tmp_i_j( tmp_i );

        for ( int16_type j=1; j < nOrder; ++j )
        {
            // 1 <= i <= nOrder-1,  1 <= j <= nOrder-1 , 0 <= k <= nOrder-1
            m[i][j].resize( nOrder,__pts.size() );
            J[i][j] = JacobiBatchEvaluation<nOrder-1,value_type>( ( 2.0*i+2.0*j+1.0 ), 1.0, __pts );
            dJ[i][j] = JacobiBatchDerivation<nOrder-1,value_type>( ( 2.0*i+2.0*j+1.0 ), 1.0, __pts );

            tmp_i_j=ublas::element_prod( tmp_i_j,tmp1 );

            ublas::row( m[i][j],0 ) = -value_type( i+j+1.0 )*tmp_i_j / value_type( 2.0 ); // -(i+j+1)/2*((1-x)/2)^(i+j)

            for ( int16_type k=1; k < nOrder; ++k ) // 1 <= k <= nOrder-1
            {
                vector_type tmp3( ublas::element_prod( tmp1, ublas::row( J[i][j],k-1 ) )
                                  + 2.0*element_prod( tmp_prod, ublas::row( dJ[i][j],k-1 ) )
                                  - ( value_type( i+j+1.0 ) ) * ublas::element_prod( tmp2,ublas::row( J[i][j],k-1 ) )  );

                ublas::row( m[i][j],k ) = ( ublas::element_prod( tmp_i_j, tmp3 ) ) / 2.0 ;
            }

        }
    }

    return m;
}


} /* Feel */

#endif /* __Principal_H */
