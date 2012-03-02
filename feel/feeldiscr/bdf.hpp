/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2006-12-30

   Copyright (C) 2006 Universit� Joseph Fourier (Grenoble)

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
   \file bdf.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-12-30
*/
#ifndef _BDF_H
#define _BDF_H

#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <boost/shared_array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/utility.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

enum BDFTimeScheme { BDF_ORDER_ONE=1, BDF_ORDER_TWO, BDF_ORDER_THREE, BDF_ORDER_FOUR, BDF_MAX_ORDER = 4 };

/**
 * \class Bdf
 * \ingroup SpaceTime
 * \brief Backward differencing formula time discretization
 *
 * A differential equation of the form
 *
 * \f$ M u' = A u + f \f$
 *
 * is discretized in time as
 *
 * \f$ M p'(t_{k+1}) = A u_{k+1} + f_{k+1} \f$
 *
 * where p denotes the polynomial of order n in t that interpolates
 * (t_i,u_i) for i = k-n+1,...,k+1.
 *
 * The approximative time derivative \f$ p'(t_{k+1}) \f$ is a linear
 * combination of state vectors u_i:
 *
 * \f$ p'(t_{k+1}) = \frac{1}{\Delta t} (\alpha_0 u_{k+1} - \sum_{i=0}^n \alpha_i u_{k+1-i} )\f$
 *
 * Thus we have
 *
 * \f$ \frac{\alpha_0}{\Delta t} M u_{k+1} = A u_{k+1} + f + M \bar{p} \f$
 *
 * with
 *
 * \f$ \bar{p} = \frac{1}{\Delta t} \sum_{i=1}^n \alpha_i u_{k+1-i} \f$
 *
 * This class stores the n last state vectors in order to be able to
 * calculate \f$ \bar{p} \f$. It also provides alpha_i
 * and can extrapolate the the new state from the n last states with a
 * polynomial of order n-1:
 *
 * \f$ u_{k+1} \approx \sum_{i=0}^{n-1} \beta_i u_{k-i} \f$
 */
template<typename SpaceType>
class Bdf
{
public:

    typedef SpaceType space_type;
    typedef boost::shared_ptr<space_type>  space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::return_type return_type;
    typedef typename element_type::value_type value_type;
    typedef boost::shared_ptr<element_type> unknown_type;
    typedef std::vector< unknown_type > unknowns_type;
    typedef typename node<value_type>::type node_type;

    /**
     * Constructor
     *
     * @param space approximation space
     * @param n order of the BDF
     */
    Bdf( space_ptrtype const& space );

    ~Bdf();


    /**
     * get the order of the BDF
     *
     *
     * @return order of the BDF
     */
    BDFTimeScheme order() const { return _M_order; }

    /**
       Initialize all the entries of the unknown vector to be derived with the
       vector u0 (duplicated)
    */
    void initialize( element_type const& u0 );

    /**
       Initialize all the entries of the unknown vector to be derived with a
       set of vectors uv0
    */
    void initialize( unknowns_type const& uv0 );

    /**
       Update the vectors of the previous time steps by shifting on the right
       the old values.
       @param u_curr current (new) value of the state vector
    */
    template<typename container_type>
    void shiftRight( typename space_type::template Element<value_type, container_type> const& u_curr );

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
    //! formula
    element_type derivate( int, scalar_type dt ) const;

    //! Returns the right hand side \f$ \bar{p} \Delta t \f$ of the time
    //! derivative formula. The timestep is taken into account elsewhere,
    //! e. g. in the mass matrix.
    element_type derivate( int n ) const;

    //! Compute the polynomial extrapolation approximation of order n-1 of
    //! u^{n+1} defined by the n stored state vectors
    element_type extrapolate( int n ) const;

    //! Return the i-th coefficient of the time derivative alpha_i
    double derivateCoefficient( int n, size_type i, double dt = 1.0 ) const;

    //! Return the i-th coefficient of the time extrapolation beta_i
    double extrapolateCoefficient( int n, size_type i, double dt = 1.0 ) const;

    //! Return a vector with the last n state vectors
    unknowns_type const& unknowns() const;

    //! Return a vector with the last n state vectors
    element_type& unknown( int i );

    template<typename container_type>
    void setUnknown( int i,  typename space_type::template Element<value_type, container_type> const& e )
    {
        _M_unknowns[i]->assign( e );
    }

    void showMe( std::ostream& __out = std::cout ) const;

private:

    //! space
    space_ptrtype _M_space;

    //! time order
    BDFTimeScheme _M_order;

    //! Size of the unknown vector
    size_type _M_size;

    //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
    std::vector<ublas::vector<double> > _M_alpha;

    //! Coefficients \f$ \beta_i \f$ of the extrapolation
    std::vector<ublas::vector<double> > _M_beta;

    //! Last n state vectors
    unknowns_type _M_unknowns;
};

template <typename SpaceType>
Bdf<SpaceType>::Bdf( space_ptrtype const& __space )
    :
    _M_space( __space ),
    _M_size( 0 ),
    _M_alpha( BDF_MAX_ORDER ),
    _M_beta( BDF_MAX_ORDER )
{
    for( int i = 0; i < BDF_MAX_ORDER; ++i )
        {
            _M_alpha[ i ].resize( i+2 );
            _M_beta[ i ].resize( i+1 );
        }

    for ( int i = 0; i < BDF_MAX_ORDER; ++i )
        {
            if (  i == 0 ) // BDF_ORDER_ONE:
                {
                    _M_alpha[i][ 0 ] = 1.; // Backward Euler
                    _M_alpha[i][ 1 ] = 1.;
                    _M_beta[i][ 0 ] = 1.; // u^{n+1} \approx u^n
                }
            else if ( i == 1 ) // BDF_ORDER_TWO:
                {
                    _M_alpha[i][ 0 ] = 3. / 2.;
                    _M_alpha[i][ 1 ] = 2.;
                    _M_alpha[i][ 2 ] = -1. / 2.;
                    _M_beta[i][ 0 ] = 2.;
                    _M_beta[i][ 1 ] = -1.;
                }
            else if ( i == 2 ) // BDF_ORDER_THREE:
                {
                    _M_alpha[i][ 0 ] = 11. / 6.;
                    _M_alpha[i][ 1 ] = 3.;
                    _M_alpha[i][ 2 ] = -3. / 2.;
                    _M_alpha[i][ 3 ] = 1. / 3.;
                    _M_beta[i][ 0 ] = 3.;
                    _M_beta[i][ 1 ] = -3.;
                    _M_beta[i][ 2 ] = 1.;
                }
            else if ( i == 3 ) /// BDF_ORDER_FOUR:
                {
                    _M_alpha[i][ 0 ] = 25. / 12.;
                    _M_alpha[i][ 1 ] = 4.;
                    _M_alpha[i][ 2 ] = -3.;
                    _M_alpha[i][ 3 ] = 4. / 3.;
                    _M_alpha[i][ 4 ] = -1. / 4.;
                    _M_beta[i][ 0 ] = 4.;
                    _M_beta[i][ 1 ] = -6.;
                    _M_beta[i][ 2 ] = 4.;
                    _M_beta[i][ 3 ] = -1.;
                }
        }

    _M_unknowns.resize( BDF_MAX_ORDER );
    for ( uint8_type __i = 0; __i < ( uint8_type )BDF_MAX_ORDER; ++__i )
        {
            _M_unknowns[__i] = unknown_type( new element_type( _M_space ) );
            _M_unknowns[__i]->zero();
        }
}


template <typename SpaceType>
Bdf<SpaceType>::~Bdf()
{}


template <typename SpaceType>
void
Bdf<SpaceType>::initialize( element_type const& u0 )
{
    std::for_each( _M_unknowns.begin(), _M_unknowns.end(), *boost::lambda::_1 = u0 );
}

template <typename SpaceType>
void
Bdf<SpaceType>::initialize( unknowns_type const& uv0 )
{
    // Check if uv0 has the right dimensions
    //FEELPP_ASSERT( uv0.size() == uint16_type(_M_order) ).error( "Initial data set are not enough for the selected BDF" );

    std::copy( uv0.begin(), uv0.end(), _M_unknowns.begin() );
}

template <typename SpaceType>
const
typename Bdf<SpaceType>::unknowns_type&
Bdf<SpaceType>::unknowns() const
{
    return _M_unknowns;
}

template <typename SpaceType>
typename Bdf<SpaceType>::element_type&
Bdf<SpaceType>::unknown( int i )
{
    //Debug() << "[Bdf::unknown] id: " << i << " l2norm = " << _M_unknowns[i]->l2Norm() << "\n";
    return *_M_unknowns[i];
}




template <typename SpaceType>
double
Bdf<SpaceType>::derivateCoefficient( int n, size_type i, double dt ) const
{
    FEELPP_ASSERT( i < size_type(n + 1) ).error( "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_alpha[n-1][ i ]/dt;
}

template <typename SpaceType>
double
Bdf<SpaceType>::extrapolateCoefficient( int n, size_type i, double dt ) const
{
    FEELPP_ASSERT( i < n ).error( "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_beta[n-1][ i ];
}

template <typename SpaceType>
void
Bdf<SpaceType>::showMe( std::ostream& __out ) const
{
#if 0
    __out << "*** BDF Time discretization of order " << _M_order << " ***"
          << "  size : " << _M_unknowns[0]->size() << "\n"
          << " alpha : " << _M_alpha << "\n"
          << "  beta : " << _M_beta << "\n";
#endif
}

template <typename SpaceType>
template<typename container_type>
void
Bdf<SpaceType>::shiftRight( typename space_type::template Element<value_type, container_type> const& __new_unk )
{
    using namespace boost::lambda;
    typename unknowns_type::reverse_iterator __it = boost::next( _M_unknowns.rbegin() );
    std::for_each( _M_unknowns.rbegin(), boost::prior( _M_unknowns.rend() ),
                   (*lambda::_1 = *(*lambda::var( __it )), ++lambda::var( __it ) ) );
    // u(t^{n}) coefficient is in _M_unknowns[0]
    *_M_unknowns[0] = __new_unk;

    /*    int i = 0;
    BOOST_FOREACH( boost::shared_ptr<element_type>& t, _M_unknowns  )
        {
            //Debug() << "[Bdf::shiftright] id: " << i << " l2norm = " << t->l2Norm() << "\n";
            ++i;
        }
    */
}

template <typename SpaceType>
typename Bdf<SpaceType>::element_type
Bdf<SpaceType>::derivate( int n, scalar_type dt ) const
{
    element_type __t( _M_space );
    __t.zero();

    __t.add( (1./dt), derivate( n ) );

    return __t;
}

template <typename SpaceType>
typename Bdf<SpaceType>::element_type
Bdf<SpaceType>::derivate( int n ) const
{
    element_type __t( _M_space );
    __t.zero();

    FEELPP_ASSERT( __t.size() == _M_space->nDof() )( __t.size() )( _M_space->nDof() ).error( "invalid space element size" );
    FEELPP_ASSERT( __t.size() == _M_unknowns[0]->size() )( __t.size() )( _M_unknowns[0]->size() ).error( "invalid space element size" );
    for ( uint8_type i = 0;i < n;++i )
        __t.add( _M_alpha[n-1][ i+1 ], *_M_unknowns[i] );

    return __t;
}

template <typename SpaceType>
typename Bdf<SpaceType>::element_type
Bdf<SpaceType>::extrapolate( int n ) const
{
    element_type __t( _M_space );
    __t.zero();

    FEELPP_ASSERT( __t.size() == _M_space->nDof() )( __t.size() )( _M_space->nDof() ).error( "invalid space element size" );
    FEELPP_ASSERT( __t.size() == _M_unknowns[0]->size() )( __t.size() )( _M_unknowns[0]->size() ).error( "invalid space element size" );

    for ( uint8_type i = 0;i < n;++i )
        __t.add(  _M_beta[n-1][ i ],  *_M_unknowns[ i ] );

    return __t;
}

}
#endif
