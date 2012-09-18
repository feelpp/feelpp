/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-03-09

  Copyright (C) 2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier Grenoble 1

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
   \file basis.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-03-09
 */
#ifndef __Basis_H
#define __Basis_H 1

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace Feel
{
/**
 * \class Basis
 * \brief Base class for basis
 *
 *  \ingroup Polynomial
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename tag, typename T>
class Basis
{
public:


    /** @name Typedefs
     */
    //@{
    typedef T value_type;
    typedef ublas::matrix<value_type,ublas::row_major> matrix_type;
#if 0
    static const uint16_type nDim = PTraits::nDim;
    static const uint16_type nOrder = PTraits::nOrder;



    typedef P sub_type;
    typedef PTraits traits_type;

    typedef typename traits_type::value_type value_type;

    typedef typename traits_type::convex_type convex_type;
    typedef typename traits_type::reference_convex_type reference_convex_type;

    typedef typename traits_type::diff_pointset_type diff_pointset_type;

    typedef typename traits_type::storage_policy storage_policy;
    typedef typename traits_type::matrix_type matrix_type;
    typedef typename traits_type::vector_matrix_type vector_matrix_type;
    typedef typename traits_type::matrix_node_type matrix_node_type;
    typedef typename traits_type::points_type points_type;
    typedef typename traits_type::node_type node_type;


#if defined( FEELPP_HAS_QD_REAL )
    typedef typename traits_type::template ChangeValueType<qd_real>::type qd_basis_type;
    typedef typename traits_type::template ChangeValueType<qd_real>::traits_type::diff_pointset_type qd_diff_pointset_type;
#endif // FEELPP_HAS_QD_REAL
#endif
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     * call differentiation matrix static construction
     */
    template<typename PrimalBasis>
    Basis( PrimalBasis const& p )
    {
        initDerivation( p );
    }

    /**
     * copy constructor
     * no need to do something, everything is static
     */
    Basis( Basis const & b )
    {}

    /**
     * destructor, nothing to do
     */
    virtual ~Basis()
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
     * \brief derivatives of Dubiner polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    static matrix_type const& d( uint16_type i )
    {
        return _S_D[i];
    }

    /**
     * \brief derivatives of Dubiner polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    static matrix_type const& derivate( uint16_type i )
    {
        return _S_D[i];
    }



    //@}



protected:

    template<typename PrimalBasis>
    static void initDerivation( PrimalBasis const& basis );

private:

    /**
     * \c true if differentation matrix initialized, \c false
     * otherwise
     */
    static  bool _S_has_derivation;

    /**
     * Derivation matrix
     * \note construct it only once per dubiner polynomials
     */
    //static std::vector<matrix_type> _S_D;
    static std::vector<matrix_type> _S_D;
};

template<typename Tag, typename T>
bool Basis<Tag, T>::_S_has_derivation = false;

template<typename Tag, typename T>
std::vector<typename Basis<Tag, T>::matrix_type>  Basis<Tag, T>::_S_D;



#if defined( FEELPP_HAS_QD_REAL)
template<typename PrimalBasis>
static void initDerivation( PrimalBasis const& basis )
{
    initDerivation( basis, mpl::bool_<boost::is_same<value_type,double>::value>() );
}

template<typename PrimalBasis>
static void initDerivation( PrimalBasis const& basis, mpl::bool_<false> )
{

    if ( _S_has_derivation == false )
    {
        _S_has_derivation = true;

        reference_convex_type refconvex;

        // constructor pointset for differentiation only in
        // the interior(1)
        diff_pointset_type diff_pts( 1 );

        matrix_type A( sub_type::evaluate( diff_pts.points() ) );
#if 1
        matrix_type D = ublas::identity_matrix<value_type>( A.size1(), A.size2()  );

        LU<matrix_type> lu( A );

        matrix_type C = lu.solve( D );

        vector_matrix_type d ( sub_type::derivate( diff_pts.points() ) );
        _S_D.resize( d.size() );

        for ( size_type i = 0; i < d.size(); ++i )
        {
            _S_D[i] = ublas::prod( d[i], C );
            glas::clean( _S_D[i] );
        }

#else
        ublas::permutation_matrix<value_type> perm( A.size1() );
        ublas::lu_factorize( A );
        _S_D = sub_type::derivate( diff_pts.points() );

        for ( size_type i = 0; i < _S_D.size(); ++i )
        {
            ublas::lu_substitute( ublas::matrix_expression<matrix_type>( _S_D[i] ), A );
            glas::clean( _S_D[i] );
        }

#endif

    }
}

template<typename PrimalBasis>
static void initDerivation( PrimalBasis const& basis, mpl::bool_<true> )
{
    if ( _S_has_derivation == false )
    {
        _S_has_derivation = true;

        typename qd_basis_type::reference_convex_type refconvex;

        qd_diff_pointset_type diff_pts( 1 );

        typename qd_basis_type::matrix_type A( qd_basis_type::evaluate( diff_pts.points() ) );
#if 1
        typename qd_basis_type::matrix_type D = ublas::identity_matrix<value_type>( A.size1(), A.size2()  );
        LU<typename qd_basis_type::matrix_type> lu( A );
        typename qd_basis_type::matrix_type C = lu.solve( D );

        typename qd_basis_type::vector_matrix_type d ( qd_basis_type::derivate( diff_pts.points() ) );
        std::vector<typename qd_basis_type::matrix_type> _qd_derivatives_matrix;

        _qd_derivatives_matrix.resize( d.size() );
        _S_D.resize( _qd_derivatives_matrix.size() );

        for ( size_type i = 0; i < d.size(); ++i )
        {
            _qd_derivatives_matrix[i] = ublas::prod( d[i], C );
            glas::clean( _qd_derivatives_matrix[i] );

            _S_D[i].resize( _qd_derivatives_matrix[i].size1(),_qd_derivatives_matrix[i].size2() );

            for ( size_type j = 0; j < _S_D[i].size1(); ++j )
                for ( size_type k = 0; k < _S_D[i].size2(); ++k )
                    _S_D[i]( j,k ) = value_type( _qd_derivatives_matrix[i]( j,k ) );

            glas::clean( _S_D[i] );
        }

#else
        ublas::permutation_matrix<typename qd_basis_type::value_type> perm( A.size1() );
        ublas::lu_factorize( A );
        typename qd_basis_type::vector_matrix_type d ( qd_basis_type::derivate( diff_pts.points() ) );

        _S_D.resize( d.size() );

        for ( size_type i = 0; i < d.size(); ++i )
        {
            ublas::lu_substitute( ublas::matrix_expression<matrix_type>( d[i] ), A );

            _S_D[i].resize( d[i].size1(),d[i].size2() );

            for ( size_type j = 0; j < _S_D[i].size1(); ++j )
                for ( size_type k = 0; k < _S_D[i].size2(); ++k )
                    _S_D[i]( j,k ) = value_type( d[i]( j,k ) );

            glas::clean( _S_D[i] );
        }

#endif
    }
}

#else
template<typename Tag, typename T>
template<typename PrimalBasis>
void
Basis<Tag, T>::initDerivation( PrimalBasis const& basis )
{
    typedef typename PrimalBasis::traits_type traits_type;
    typedef typename traits_type::value_type value_type;

    typedef typename traits_type::convex_type convex_type;
    typedef typename traits_type::reference_convex_type reference_convex_type;

    typedef typename traits_type::diff_pointset_type diff_pointset_type;

    typedef typename traits_type::storage_policy storage_policy;
    typedef typename traits_type::matrix_type matrix_type;
    typedef typename traits_type::vector_matrix_type vector_matrix_type;
    typedef typename traits_type::matrix_node_type matrix_node_type;
    typedef typename traits_type::points_type points_type;
    typedef typename traits_type::node_type node_type;

    if ( _S_has_derivation == false )
    {
        _S_has_derivation = true;

        reference_convex_type refconvex;
        // constructor pointset for differentiation only in
        // the interior(1)
        diff_pointset_type diff_pts( 1 );
        matrix_type A( basis.evaluate( diff_pts.points() ) );

#if 1
        matrix_type D = ublas::identity_matrix<value_type>( A.size1(), A.size2()  );
        LU<matrix_type> lu( A );
        matrix_type C = lu.solve( D );

        vector_matrix_type d ( basis.derivate( diff_pts.points() ) );
        _S_D.resize( d.size() );

        for ( size_type i = 0; i < d.size(); ++i )
        {
            _S_D[i] = ublas::prod( d[i], C );
            glas::clean( _S_D[i] );
        }

#else
        ublas::permutation_matrix<double> perm( A.size1() );
        ublas::lu_factorize( A );
        _S_D = basis.derivate( diff_pts.points() );

        for ( size_type i = 0; i < d.size(); ++i )
        {
            ublas::lu_substitute( _S_D[i], A );
            glas::clean( _S_D[i] );
        }

#endif
    }
}
#endif // FEELPP_HAS_QD_REAL


} // Feel
#endif /* __Basis_H */
