/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-29

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
   \file expansions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-29
 */
#ifndef __expansions_H
#define __expansions_H 1

#include <feel/feelalg/glas.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feelpoly/jacobi.hpp>

namespace Feel
{
enum TheShape { LINE = 1, TRIANGLE = 2, TETRAHEDRON = 3 };

template<TheShape sh>
struct DimFromShape
{
    static const uint16_type value = mpl::if_<mpl::equal_to<mpl::int_<sh>, mpl::int_<LINE> >,
                             mpl::int_<1>,
                             typename mpl::if_<mpl::equal_to<mpl::int_<sh>, mpl::int_<TRIANGLE> >,
                             mpl::int_<2>,
                             mpl::int_<3> >::type>::type::value;
};

/// \cond DETAIL
/**
 * expansions function as described in Karniadakis/Sherwin book
 */
namespace details
{

/**
 * collapsed coordinates
 */
template<TheShape sh,  typename T = double>
struct xi
{
};

template<typename T>
struct xi<TRIANGLE, T>
{
    typedef T value_type;
    typedef typename node<value_type>::type node_type;

    xi()
        :
        M_xi( 2 )
    {}

    xi( node_type const& eta )
        :
        M_xi( 2 )
    {
        M_xi[0] = 0.5*( 1.0 + eta[0] )*( 1.0 - eta[1] ) - 1.0;
        M_xi[1] = eta[1];
    }

    node_type const& operator()( node_type const& eta )
    {
        M_xi[0] = 0.5*( 1.0 + eta[0] )*( 1.0 - eta[1] ) - 1.0;
        M_xi[1] = eta[1];
        return M_xi;
    }
    node_type const& operator()() const
    {
        return M_xi;
    }
    node_type M_xi;
};

template<typename T>
struct xi<TETRAHEDRON, T>
{
    typedef T value_type;
    typedef typename node<value_type>::type node_type;

    xi()
        :
        M_xi( 3 )
    {}

    xi( node_type const& eta )
        :
        M_xi( 3 )
    {
        M_xi[0] = 0.25*( 1.0 + eta[0] )*( 1.0 - eta[1] )*( 1.0 - eta[2] ) - 1.0;
        M_xi[1] = 0.5*( 1.0 + eta[1] )*( 1.0 - eta[2] ) - 1.0;
        M_xi[2] = eta[2];
    }

    node_type const& operator()( node_type const& eta )
    {
        M_xi[0] = 0.25*( 1.0 + eta[0] )*( 1.0 - eta[1] )*( 1.0 - eta[2] ) - 1.0;
        M_xi[1] = 0.5*( 1.0 + eta[1] )*( 1.0 - eta[2] ) - 1.0;
        M_xi[2] = eta[2];
        return M_xi;
    }
    node_type const& operator()() const
    {
        return M_xi;
    }
    node_type M_xi;
};


/**
 * collapsed coordinates
 */
template<TheShape sh,  typename T = double>
struct eta
{
};
template<TheShape sh,  typename T = double>
struct etas
{
};

template<typename T>
struct eta<TRIANGLE, T>
{
    typedef T value_type;
    typedef typename node<value_type>::type node_type;

    eta()
        :
        M_eta( 2 )
    {}

    eta( node_type const& xi )
        :
        M_eta( 2 )
    {
        if ( xi[1] == 1.0 )
            M_eta[0] = -1.0;

        else
            M_eta[0] = 2.0 * ( 1.0 + xi[0] ) / ( 1.0 - xi[1] ) - 1.0;

        M_eta[1] = xi[1];
    }

    node_type const& operator()( node_type const& xi )
    {
        if ( xi[1] == 1.0 )
            M_eta[0] = -1.0;

        else
            M_eta[0] = 2.0 * ( 1.0 + xi[0] ) / ( 1.0 - xi[1] ) - 1.0;

        M_eta[1] = xi[1];
        return M_eta;
    }
    node_type const& operator()() const
    {
        return M_eta;
    }
    node_type M_eta;
};
template<typename T>
struct etas<TRIANGLE, T>
{
    typedef T value_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;

    etas()
        :
        M_eta()
    {}

    etas( matrix_node_type const& xi )
        :
        M_eta( xi.size1(), xi.size2() )
    {
        for ( size_type i = 0; i < xi.size2(); ++i )
        {
            if ( xi( 1, i ) == 1.0 )
                M_eta( 0, i ) = -1.0;

            else
                M_eta( 0, i ) = 2.0 * ( 1.0 + xi( 0, i ) ) / ( 1.0 - xi( 1, i ) ) - 1.0;

            M_eta( 1, i ) = xi( 1, i );
        }

    }

    matrix_node_type const& operator()( matrix_node_type const& xi )
    {
        M_eta.resize( xi.size1(), xi.size2() );

        for ( size_type i = 0; i < xi.size2(); ++i )
        {
            if ( xi( 1, i ) == 1.0 )
                M_eta( 0, i ) = -1.0;

            else
                M_eta( 0, i ) = 2.0 * ( 1.0 + xi( 0, i ) ) / ( 1.0 - xi( 1, i ) ) - 1.0;

            M_eta( 1, i ) = xi( 1, i );
        }

        return M_eta;
    }
    matrix_node_type const& operator()() const
    {
        return M_eta;
    }
    matrix_node_type M_eta;
};

template<typename T>
struct etas<TETRAHEDRON, T>
{
    typedef T value_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;

    etas()
        :
        M_eta()
    {}

    etas( matrix_node_type const& xi )
        :
        M_eta( xi.size1(), xi.size2() )
    {
        for ( size_type i = 0; i < xi.size2(); ++i )
        {
            if ( xi( 1, i ) + xi( 2, i ) == 0. )
                M_eta( 0, i ) = 1.;

            else
                M_eta( 0, i ) = -2. * ( 1. + xi( 0, i ) ) / ( xi( 1, i ) + xi( 2, i ) ) - 1.;

            if ( xi( 2, i ) == 1. )
                M_eta( 1, i ) = -1.;

            else
                M_eta( 1, i ) = 2. * ( 1. + xi( 1, i ) ) / ( 1. - xi( 2, i ) ) - 1.;

            M_eta( 2, i ) = xi( 2, i );
        }
    }

    matrix_node_type const& operator()( matrix_node_type const& xi )
    {
        M_eta.resize( xi.size1(), xi.size2() );

        for ( size_type i = 0; i < xi.size2(); ++i )
        {
            if ( xi( 1, i ) + xi( 2, i ) == 0. )
                M_eta( 0, i ) = 1.;

            else
                M_eta( 0, i ) = -2. * ( 1. + xi( 0, i ) ) / ( xi( 1, i ) + xi( 2, i ) ) - 1.;

            if ( xi( 2, i ) == 1. )
                M_eta( 1, i ) = -1.;

            else
                M_eta( 1, i ) = 2. * ( 1. + xi( 1, i ) ) / ( 1. - xi( 2, i ) ) - 1.;

            M_eta( 2, i ) = xi( 2, i );
        }

        return M_eta;
    }
    matrix_node_type const& operator()() const
    {
        return M_eta;
    }
    matrix_node_type M_eta;
};

/**
 * \class psitilde
 * \brief implements parts of the Dubiner polynomials
 *
 * psitilde allows to compute parts of the Dubnier polynomials namely
 * \f$ \left(\frac{1-\eta_2}{2}\right)^i P^{1 i + 1,0}_j(\eta_2) \f$
 * where \f$ P^{2 i + 1,0}_j(\eta_2) \f$ is the Jacobi polynomial of
 * degree \f$ j\f$ with weights \f$ 2 i + 1 \f$ and \f$ 0 \f$.
 *
 * The cordinates must be cartesian : in a simplex one must then
 * change from wrapped coordinates to cartesian coordinates
 */
template<typename T>
class psitilde
{
public:
    typedef T value_type;
    typedef Feel::dyna::Jacobi<T> P;
    psitilde()
        :
        p( 0, 0.0, 0.0 ),
        M_b( 0.0 )
    {}
    psitilde( int i )
        :
        p( i, 0.0,  0.0 ),
        M_b( i )
    {
    }
    psitilde( int i, int j )
        :
        p( j, 2*i+1, 0.0 ),
        M_b( i )
    {
    }
    psitilde( int i, int j, int k )
        :
        p( k, 2*( i+j+1 ), 0.0 ),
        M_b( i+j )
    {
    }
    /**
     * \brief 1D case
     * \arg x coordinate in [-1;1]
     * \return \f$ P^{0,0}_i(x) \f$
     */
    value_type a( value_type const& x ) const
    {
        return p( x );
    }
    /**
     * \brief 2D case
     * \arg x coordinate in [-1;1]
     * \return \f$ \left(\frac{1-x}{2}\right)^{i} P^{2 i + 1),0}_j(x) \f$
     */
    value_type b( value_type const& x ) const
    {
        return math::pow( 0.5 *( 1.0-x ), M_b )*p( x );
    }
    /**
     * \brief 3D case
     * \arg x coordinate in [-1;1]
     * \return \f$ \left(\frac{1-x}{2}\right)^{i+j} P^{2( i + j + 1),0}_k(x) \f$
     */
    value_type c( value_type const& x ) const
    {
        return math::pow( 0.5 *( 1.0-x ), M_b )*p( x );
    }
    P p;
    value_type M_b;
};

/**
 * \class scalings
 * \brief scaling factors for Dubiner basis
 *
 * \author Christophe Prud'homme
 */
template<uint16_type N, typename T = double>
struct scalings
{
    typedef T value_type;
    typedef ublas::matrix<value_type> matrix_type;
    typedef typename matrix_node<value_type>::type matrix_node_type;

    /**
     * \brief evaluates scaling factors for Dubiner polynomials
     *
     * pts are coordinates in [-1;1]
     */
    scalings( ublas::vector<value_type> const& pts )
        :
        M_s( N+1,  pts.size() )
    {
#if 0
        ublas::vector<value_type> one( ublas::scalar_vector<value_type>( pts.size(), 1.0 ) );
        ublas::row( M_s, 0 ) = one;
#else
        ublas::row( M_s, 0 ) = ublas::scalar_vector<value_type>( pts.size(), 1.0 );
#endif

        if ( N > 0 )
        {
            ublas::row( M_s, 1 ) = 0.5 * ( ublas::row( M_s, 0 ) - pts );

            for ( uint16_type k = 2; k < N+1; ++k )
            {
                ublas::row( M_s, k ) = ublas::element_prod( ublas::row( M_s, k-1 ),
                                        ublas::row( M_s, 1 ) );
            }
        }
    }
    matrix_type const& operator()() const
    {
        return M_s;
    }

    matrix_type M_s;
};
} // details
/// \endcond


}

#endif /* __expansions_H */
