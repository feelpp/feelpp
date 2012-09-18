/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file dirscalingmatrix.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef __DirScalingMatrix_H
#define __DirScalingMatrix_H 1

#include <algorithm>


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace Feel
{
using namespace boost::numeric::ublas;

/*!
  \class DirScalingMatrix
  \brief implements the directional Scaling Matrix for directionally-scaled trust region

  @author Christophe Prud'homme
  @see
*/
template<typename NumType>
class DirScalingMatrix
{
public:


    /** @name Typedefs
     */
    //@{

    typedef enum mode_type { NO_JACOBIAN, WITH_JACOBIAN };

    typedef NumType value_type;

    typedef vector<value_type> vector_type;
    typedef banded_matrix<value_type> matrix_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    DirScalingMatrix()
        :
        M_lb(),
        M_ub(),
        M_lb_ub(),
        M_value(),
        M_jacobian(),
        M_trust_region_active( false )
    {}

    DirScalingMatrix( vector_type const& __lb, vector_type const& __ub )
        :
        M_lb( __lb ),
        M_ub( __ub ),
        M_lb_ub( M_ub - M_lb ),
        M_value( __lb.size(), __lb.size(), 0, 0 ),
        M_jacobian( __lb.size(), __lb.size(), 0, 0 ),
        M_trust_region_active( false )
    {
        GST_SMART_ASSERT( __lb.size() == __ub.size() )( __lb )( __ub )( "inconsistent bounds definition" );
        GST_SMART_ASSERT( *std::min_element( M_lb_ub.begin(), M_lb_ub.end() ) >= 0 )
        ( M_lb )( M_ub )( "lower and upper bounds are not properly defined" );
    }
    DirScalingMatrix( DirScalingMatrix const & )
    {}
    ~DirScalingMatrix()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    value_type zeta( vector_type const& __x ) const;
    value_type zeta() const
    {
        return M_zeta;
    }

    matrix_type const& operator()() const
    {
        return M_value;
    }

    matrix_type const& jacobian() const
    {
        return M_jacobian;
    }

    bool isTrustRegionActive() const
    {
        return M_trust_region_active;
    }
    //@}

    /** @name  Mutators
     */
    //@{

    void update( value_type const&,  vector_type const&, vector_type const&, mode_type = WITH_JACOBIAN );

    void setBounds( vector_type const& __lb, vector_type const& __up )
    {
        GST_SMART_ASSERT( __lb.size() == __up.size() )( __lb )( __up )( "inconsistent bounds definition" );
        M_lb = __lb;
        M_ub = __up;
        M_lb_ub = __up - __lb;

        GST_SMART_ASSERT( *std::min_element( M_lb_ub.begin(), M_lb_ub.end() ) >= 0 )
        ( M_lb )( M_ub )( "lower and upper bounds are not properly defined" );
    }
    //@}

    /** @name  Methods
     */
    //@{

    //@}



protected:

    vector_type distanceToLB( vector_type const& __x ) const
    {
        GST_SMART_ASSERT( __x.size() == M_lb.size() )( __x )( M_lb )( "inconsistent bounds definition" );
        return element_div( __x - M_lb, M_lb_ub );
    }

    vector_type distanceToUB( vector_type const& __x ) const
    {
        GST_SMART_ASSERT( __x.size() == M_ub.size() )( __x )( M_ub )( "inconsistent bounds definition" );
        return element_div( __x - M_ub, M_lb_ub );
    }

private:

    vector_type M_lb;
    vector_type M_ub;
    vector_type M_lb_ub;

    matrix_type M_value;
    matrix_type M_jacobian;

    mutable value_type M_zeta;

    bool M_trust_region_active;
};

template<typename NumType>
typename DirScalingMatrix<NumType>::value_type
DirScalingMatrix<NumType>::zeta( vector_type const& __x ) const
{
    M_zeta = 0.9;//M_options.zeta_min;

    M_zeta = std::max( M_zeta, std::max( ublas::norm_inf( distanceToLB( __x ) ),
                                         ublas::norm_inf( distanceToUB( __x ) ) ) );
    return M_zeta;
}
template<typename NumType>
void
DirScalingMatrix<NumType>::update( value_type const& __Delta,
                                   vector_type const& __x,
                                   vector_type const& __s,
                                   mode_type __mode )
{
    GST_SMART_ASSERT( __x.size() == M_lb.size() )( __x )( M_lb )( "inconsistent bounds definition" );
    GST_SMART_ASSERT( __x.size() == M_ub.size() )( __x )( M_ub )( "inconsistent bounds definition" );

    M_value.resize( __x.size(), __x.size(), 0, 0 );
    M_jacobian.resize( __x.size(), __x.size(), 0, 0 );


    M_zeta = zeta( __x );

    vector_type __dl = distanceToLB( __x );
    vector_type __du = distanceToUB( __x );

    if ( M_zeta * std::min( norm_inf( __dl ), norm_inf( __du ) ) > __Delta )
    {
        // we are in the trust region
        M_value = identity_matrix<value_type>( M_value.size1(), M_value.size2() );
    }

    else
    {
        M_trust_region_active = true;

        for ( size_t __i = 0; __i < __x.size(); ++__i )
        {
            if ( __s ( __i ) < 0 )
                M_value ( __i, __i ) = M_zeta * std::min( 1. , __dl( __i ) ) / __Delta;

            else
                M_value ( __i, __i ) = M_zeta * std::min( 1. , __du( __i ) ) / __Delta;
        }
    }

    if ( __mode == WITH_JACOBIAN )
    {
        for ( size_t i = 0; i < __x.size(); i++ )
        {
            if ( ( __s( i ) < 0 ) && ( __x( i ) < M_lb( i ) + __Delta ) )
            {
                M_jacobian ( i, i ) = M_zeta / __Delta;
            }

            else if	( ( __s ( i ) > 0 ) && ( __x ( i ) > M_ub( i ) - __Delta ) )
            {
                M_jacobian ( i, i ) = -M_zeta / __Delta;
            }

            else
            {
                M_jacobian ( i, i ) = 0;
            }
        }
    }
}

}
#endif /* __DirScalingMatrix_H */

