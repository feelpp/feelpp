/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-06-23

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
   \file imsimplex.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-06-23
 */
#ifndef __IMSimplex_H
#define __IMSimplex_H 1

#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelpoly/gauss.hpp>


namespace Feel
{
namespace ublas = boost::numeric::ublas;

/// \cond detail
/**
 * \internal
 * \namespace Feel::detail
 */
namespace detail
{
using namespace boost::assign;

/**
 * \internal
 * \class IMTetrahedra
 */
template<int Order, typename T>
struct IMTriangle
{
};

template<typename T>
struct IMTriangle<1,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 1;
    static const uint16_type nPoints = 1;

    IMTriangle()
        :
        q ( 1 )
    {
        m += 1;
        w += 1;
        push_back( q[0] ).repeat( 3, value_type( 1 )/value_type( 3 ) );
    }


    std::vector<uint16_type> m;
    std::vector<value_type> w;
    std::vector<std::vector<value_type> > q;
};

template<typename T> struct IMTriangle<0,T>: public IMTriangle<1,T> {};

template<typename T>
struct IMTriangle<2,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 2;
    static const uint16_type nPoints = 3;

    IMTriangle()
        :
        q( 1 )
    {
        m += 3;
        w += value_type( 1 )/value_type( 3 );
        q[0] +=  value_type( 2 )/value_type( 3 ), repeat( 2, value_type( 1 )/value_type( 6 ) );
    }


    std::vector<uint16_type> m;
    std::vector<value_type> w;
    std::vector<std::vector<value_type> > q;
};

template<>
struct IMTriangle<3,double>
{
    typedef double value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 3;
    static const uint16_type nPoints = 6;
    IMTriangle()
        :
        q( 1 )
    {
        m += 6;
        w += value_type( 1 )/value_type( 6 );
        q[0] +=  value_type( 0.659027622374092 ), value_type( 0.231933368553031 ), value_type( 0.109039009072877 );
    }

    std::vector<uint16_type> m;
    std::vector<value_type> w;
    std::vector<std::vector<value_type> > q;
};

template<>
struct IMTriangle<4,double>
{
    typedef double value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 4;
    static const uint16_type nPoints = 6;
    IMTriangle()
        :
        q( 2 )
    {
        m += 3,3;
        w += 0.109951743655322, 0.223381589678011;
        q[0] += 0.816847572980459, 0.091576213509771, 0.091576213509771;
        q[1] += 0.108103018168070, 0.445948490915965, 0.445948490915965;
    }

    std::vector<uint16_type> m;
    std::vector<value_type> w;
    std::vector<std::vector<value_type> > q;
};

template<>
struct IMTriangle<5,double>
{
    typedef double value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 5;
    static const uint16_type nPoints = 7;
    IMTriangle()
        :
        q( 3 )
    {
        m += 1, 3, 3;
        w += 0.22500000000000, 0.125939180544827, 0.132394152788506;
        q[0] += value_type( 1.0 )/value_type( 3.0 ), value_type( 1.0 )/value_type( 3.0 ), value_type( 1.0 )/value_type( 3.0 );
        q[1] += 0.797426985353087, 0.101286507323456, 0.101286507323456;
        q[2] += 0.059715871789770, 0.470142064105115, 0.470142064105115;
    }

    std::vector<uint16_type> m;
    std::vector<value_type> w;
    std::vector<std::vector<value_type> > q;
};

template<>
struct IMTriangle<6,double>
{
    typedef double value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 6;
    static const uint16_type nPoints = 12;
    IMTriangle()
        :
        q( 3 )
    {
        m += 3, 3, 6;
        w += 0.050844906370207, 0.116786275726379, 0.082851075618374;
        q[0] += 0.873821971016996, 0.063089014491502, 0.063089014491502;
        q[1] += 0.501426509658179, 0.249286745170910, 0.249286745170910;
        q[2] += 0.636502499121399, 0.310352451033784, 0.053145049844817;
    }

    std::vector<uint16_type> m;
    std::vector<value_type> w;
    std::vector<std::vector<value_type> > q;
};

/// \endcond

/**
 * \internal
 * \class IMTetrahedra
 */
template<int Order, typename T>
struct IMTetrahedra
{
};

template<typename T>
struct IMTetrahedra<1,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 1;
    static const uint16_type nPoints = 1;

    IMTetrahedra()
        :
        q ( 1 )
    {
        m += 1;
        w += 1;
        q[0] += value_type( 1 )/value_type( 4 );
    }

    std::vector<uint16_type> m;
    std::vector<value_type> w;
    std::vector<std::vector<value_type> > q;
};
template<typename T> struct IMTetrahedra<0,T>: public IMTetrahedra<1,T> {};
template<typename T>
struct IMTetrahedra<2,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 2;
    static const uint16_type nPoints = 4;
    IMTetrahedra()
        :
        q( 1 )
    {
        m += 4;
        w += value_type( 1 )/value_type( 4 );
        q[0] += 0.5854101966249685, 0.1381966011250105;
    }

    std::vector<uint16_type> m;
    std::vector<value_type> w;
    std::vector<std::vector<value_type> > q;
};

template<typename T>
struct IMTetrahedra<4,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 4;
    static const uint16_type nPoints = 16;
    IMTetrahedra()
        :
        q( 2 )
    {
        m += 4, 12;
        w += 0.05037379410012282, 0.06654206863329239;
        q[0] += 0.7716429020672371, 0.7611903264425430e-01;
        q[1] += 0.1197005277978019, 0.7183164526766925e-01, 0.4042339134672644;
    }

    std::vector<uint16_type> m;
    std::vector<value_type> w;
    std::vector<std::vector<value_type> > q;
};
template<typename T> struct IMTetrahedra<3,T>: public IMTetrahedra<4,T> {};
template<typename T>
struct IMTetrahedra<6,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 6;
    static const uint16_type nPoints = 29;
    IMTetrahedra()
        :
        q( 4 )
    {
        m += 1, 4, 12, 12;
        w += 0.9040129046014750e-01, 0.1911983427899124e-01,
             0.4361493840666568e-01, 0.2581167596199161e-01;

        q[0] += 0.25;
        q[1] += 0.8277192480479295, 0.5742691731735683e-01;
        q[2] += 0.5135188412556341e-01, 0.4860510285706072, 0.2312985436519147;
        q[3] += 0.2967538129690260, 0.6081079894015281, 0.4756909881472290e-01;
    }

    std::vector<uint16_type> m;
    std::vector<value_type> w;
    std::vector<std::vector<value_type> > q;
};
template<typename T> struct IMTetrahedra<5,T>: public IMTetrahedra<6,T> {};
} // detail


/**
 * \class IMSimplex
 * \brief brief description
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
 */
template<int Dim,int Order, typename T>
class IMSimplex
    :
public PointSetQuadrature<Simplex<Dim,1> , Order, T>
{
    typedef PointSetQuadrature<Simplex<Dim,1> , Order, T> super;
public:

    /** @name Typedefs
     */
    //@{
    typedef T value_type;
    typedef ublas::matrix<value_type,ublas::column_major> matrix_type;
    typedef ublas::vector<value_type> vector_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<Dim>,mpl::int_<2> >,
            mpl::identity<detail::IMTriangle<Order,T> >,
            mpl::identity<detail::IMTetrahedra<Order,T> > >::type::type quad_type;

#if 1
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<Dim>,mpl::int_<2> >,
            mpl::identity<Gauss<Simplex<Dim-1,1>,Order,T> >,
            mpl::identity<IMSimplex<Dim-1,Order,T> > >::type::type face_quad_type;
#else
    typedef IMSimplex<Dim-1,Order,T> face_quad_type;
#endif
    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;

    //@}

    /** @name Constructors, destructor
     */
    //@{
    IMSimplex( bool transform = true )
        :
        super( quad_type::nPoints ),
        M_quad()
    {
        //std::cout << "npoints=" << quad_type::nPoints << "\n";
        for ( size_type i=0, wi = 0; i< M_quad.q.size(); i++ )
        {
            permute( M_quad.m[i], M_quad.q[i], wi );

            for ( int j=0; j<M_quad.m[i]; j++, wi++ )
            {
                this->M_w( wi ) = factor()*M_quad.w[i];
            }
        }


        if ( transform )
        {
            //std::cout << "order = " << Order << " | " << this->M_points << "\n";
            typedef GeoND<nDim,GeoEntity<Simplex<nDim,1> >, value_type > element_type;
            RealToReference<element_type, Simplex,value_type> real_to_ref( element() );
            this->setPoints(  real_to_ref( this->M_points ) );
            this->M_w /= real_to_ref.J();
            //std::cout << "jac = " << real_to_ref.J() << "\n";
        }

        //std::cout << "this->M_points=" << this->M_points << "\n"
        //<< "  w=" << this->M_w << "\n";

        boost::shared_ptr<GT_Lagrange<Dim,1,Simplex,T> > gm( new GT_Lagrange<Dim, 1, Simplex, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Simplex<Dim, 1>, Dim, 1>(), gm, face_qr );
    }

    ~IMSimplex()
    {}

    //@}

    T factor() const
    {
        return ( Dim==2 )?T( 1 )/T( 2 ):T( 1 )/T( 6 );
    }


    /** @name  Methods
     */
    //@{
    GeoND<nDim,GeoEntity<Simplex<nDim,1> >, value_type >  element() const
    {
        return element( mpl::int_<nDim>() );
    }

    GeoND<nDim,GeoEntity<Simplex<nDim,1> >, value_type >  element( mpl::int_<2> ) const
    {
        // setup the convex
        typename node<value_type>::type pt( nDim );
        typedef GeoND<nDim,GeoEntity<Simplex<nDim,1> >, value_type > element_type;
        typedef typename element_type::point_type point_type;
        element_type elem;
        // lower left
        pt[0] = value_type( 0 );
        pt[1] = value_type( 0 );
        point_type p1( pt );
        elem.setPoint( 0, p1 );

        // lower right
        pt[0] = value_type( 1 );
        pt[1] = value_type( 0 );
        point_type p2( pt );
        elem.setPoint( 1, p2 );

        // upper
        pt[0] = value_type( 0 );
        pt[1] = value_type( 1 );
        point_type p3( pt );
        elem.setPoint( 2, p3 );
        return elem;

    }

    GeoND<nDim,GeoEntity<Simplex<nDim,1> >, value_type >  element( mpl::int_<3> ) const
    {
        // setup the convex
        typename node<value_type>::type pt( nDim );
        typedef GeoND<nDim,GeoEntity<Simplex<nDim,1> >, value_type > element_type;
        typedef typename element_type::point_type point_type;
        element_type elem;
        // lower left
        pt[0] = value_type( 0 );
        pt[1] = value_type( 0 );
        pt[2] =  value_type( 0 );
        point_type p1( pt );
        elem.setPoint( 0, p1 );

        // lower right
        pt[0] = value_type( 1 );
        pt[1] = value_type( 0 );
        pt[2] = value_type( 0 );
        point_type p2( pt );
        elem.setPoint( 1, p2 );

        //
        pt[0] = value_type( 0 );
        pt[1] = value_type( 1 );
        pt[2] = value_type( 0 );
        point_type p3( pt );
        elem.setPoint( 2, p3 );

        //
        pt[0] = value_type( 0 );
        pt[1] = value_type( 0 );
        pt[2] = value_type( 1 );
        point_type p4( pt );
        elem.setPoint( 3, p4 );
        return elem;

    }

    template<typename QVec>
    void
    permute( int m, QVec const& q, int wi )
    {
        permute( m, q, wi, mpl::int_<Dim>() );
    }
    template<typename QVec>
    void
    permute( int m, QVec const& q, int wi, mpl::int_<2> )
    {
        switch ( m )
        {
        case 1:
        {
            this->M_points( 0, wi ) = q[ 0 ];
            this->M_points( 1, wi ) = q[ 1 ];
            wi++;
        }
        break;

        case 3:
        {
            this->M_points( 0, wi ) = q[ 0 ];
            this->M_points( 1, wi ) = q[ 1 ];
            wi++;

            this->M_points( 0, wi ) = q[ 1 ];
            this->M_points( 1, wi ) = q[ 0 ];
            wi++;

            this->M_points( 0, wi ) = q[ 2 ];
            this->M_points( 1, wi ) = q[ 1 ];
            wi++;

        }
        break;

        case 6:
        {
            //qPerm[0] = tuple(q[0], q[1], q[2]);
            this->M_points( 0, wi ) = q[ 0 ];
            this->M_points( 1, wi ) = q[ 1 ];
            wi++;

            //qPerm[1] = tuple(q[0], q[2], q[1]);
            this->M_points( 0, wi ) = q[ 0 ];
            this->M_points( 1, wi ) = q[ 2 ];
            wi++;

            //qPerm[2] = tuple(q[1], q[0], q[2]);
            this->M_points( 0, wi ) = q[ 1 ];
            this->M_points( 1, wi ) = q[ 0 ];
            wi++;

            //qPerm[3] = tuple(q[1], q[2], q[0]);
            this->M_points( 0, wi ) = q[ 1 ];
            this->M_points( 1, wi ) = q[ 2 ];
            wi++;

            //qPerm[4] = tuple(q[2], q[1], q[0]);
            this->M_points( 0, wi ) = q[ 2 ];
            this->M_points( 1, wi ) = q[ 1 ];
            wi++;

            //qPerm[5] = tuple(q[2], q[0], q[1]);
            this->M_points( 0, wi ) = q[ 2 ];
            this->M_points( 1, wi ) = q[ 0 ];
            wi++;
        }
        break;
        }
    }
    template<typename QVec>
    void
    permute( int m, QVec const& q, int wi, mpl::int_<3> )
    {
        if ( m==1 )
        {
            //qPerm[0] = boost::make_tuple(q[0], q[0], q[0], q[0]);
            this->M_points( 0, wi ) = q[ 0 ];
            this->M_points( 1, wi ) = q[ 0 ];
            this->M_points( 2, wi ) = q[ 0 ];
        }

        else if ( m==4 )
        {
            //qPerm[0] = boost::make_tuple(q[0], q[1], q[1], q[1]);
            this->M_points( 0, wi ) = q[ 0 ];
            this->M_points( 1, wi ) = q[ 1 ];
            this->M_points( 2, wi ) = q[ 1 ];
            wi++;

            //qPerm[1] = boost::make_tuple(q[1], q[0], q[1], q[1]);
            this->M_points( 0, wi ) = q[ 1 ];
            this->M_points( 1, wi ) = q[ 0 ];
            this->M_points( 2, wi ) = q[ 1 ];
            wi++;

            //qPerm[2] = boost::make_tuple(q[1], q[1], q[0], q[1]);
            this->M_points( 0, wi ) = q[ 1 ];
            this->M_points( 1, wi ) = q[ 1 ];
            this->M_points( 2, wi ) = q[ 0 ];
            wi++;

            //qPerm[3] = boost::make_tuple(q[1], q[1], q[1], q[0]);
            this->M_points( 0, wi ) = q[ 1 ];
            this->M_points( 1, wi ) = q[ 1 ];
            this->M_points( 2, wi ) = q[ 1 ];
            wi++;
        }

        else if ( m==12 )
        {
            //qPerm[0] = boost::make_tuple(q[0], q[1], q[2], q[2]);
            this->M_points( 0, wi ) = q[ 0 ];
            this->M_points( 1, wi ) = q[ 1 ];
            this->M_points( 2, wi ) = q[ 2 ];
            wi++;

            //qPerm[1] = boost::make_tuple(q[0], q[2], q[1], q[2]);
            this->M_points( 0, wi ) = q[ 0 ];
            this->M_points( 1, wi ) = q[ 2 ];
            this->M_points( 2, wi ) = q[ 1 ];
            wi++;

            //qPerm[2] = boost::make_tuple(q[0], q[2], q[2], q[1]);
            this->M_points( 0, wi ) = q[ 0 ];
            this->M_points( 1, wi ) = q[ 2 ];
            this->M_points( 2, wi ) = q[ 2 ];
            wi++;


            //qPerm[3] = boost::make_tuple(q[1], q[0], q[2], q[2]);
            this->M_points( 0, wi ) = q[ 1 ];
            this->M_points( 1, wi ) = q[ 0 ];
            this->M_points( 2, wi ) = q[ 2 ];
            wi++;

            //qPerm[4] = boost::make_tuple(q[2], q[0], q[1], q[2]);
            this->M_points( 0, wi ) = q[ 2 ];
            this->M_points( 1, wi ) = q[ 0 ];
            this->M_points( 2, wi ) = q[ 1 ];
            wi++;

            //qPerm[5] = boost::make_tuple(q[2], q[0], q[2], q[1]);
            this->M_points( 0, wi ) = q[ 2 ];
            this->M_points( 1, wi ) = q[ 0 ];
            this->M_points( 2, wi ) = q[ 2 ];
            wi++;

            //qPerm[6] = boost::make_tuple(q[1], q[2], q[0], q[2]);
            this->M_points( 0, wi ) = q[ 1 ];
            this->M_points( 1, wi ) = q[ 2 ];
            this->M_points( 2, wi ) = q[ 0 ];
            wi++;

            //qPerm[7] = boost::make_tuple(q[2], q[1], q[0], q[2]);
            this->M_points( 0, wi ) = q[ 2 ];
            this->M_points( 1, wi ) = q[ 1 ];
            this->M_points( 2, wi ) = q[ 0 ];
            wi++;

            //qPerm[8] = boost::make_tuple(q[2], q[2], q[0], q[1]);
            this->M_points( 0, wi ) = q[ 2 ];
            this->M_points( 1, wi ) = q[ 2 ];
            this->M_points( 2, wi ) = q[ 0 ];
            wi++;

            //qPerm[9] = boost::make_tuple(q[1], q[2], q[2], q[0]);
            this->M_points( 0, wi ) = q[ 1 ];
            this->M_points( 1, wi ) = q[ 2 ];
            this->M_points( 2, wi ) = q[ 2 ];
            wi++;

            //qPerm[10] = boost::make_tuple(q[2], q[1], q[2], q[0]);
            this->M_points( 0, wi ) = q[ 2 ];
            this->M_points( 1, wi ) = q[ 1 ];
            this->M_points( 2, wi ) = q[ 2 ];
            wi++;

            //qPerm[11] = boost::make_tuple(q[2], q[2], q[1], q[0]);
            this->M_points( 0, wi ) = q[ 2 ];
            this->M_points( 1, wi ) = q[ 2 ];
            this->M_points( 2, wi ) = q[ 1 ];
            wi++;
        }

        else
        {
            FEELPP_ASSERT( 0 )( m ).error( "invalid multiplicity in IMSimplex<3>" );
        }

    }

    bool test()
    {
        return test( mpl::int_<nDim>() );
    }
    bool test( mpl::int_<2> )
    {
        bool pass = true;
        ublas::vector<value_type> x( ublas::row( this->M_points, 0 ) );
        ublas::vector<value_type> y( ublas::row( this->M_points, 1 ) );
        ublas::vector<value_type> w( this->M_w );

        for ( int a=0; a<=nOrder; a++ )
        {
            for ( int b=0; b<nOrder-a; b++ )
            {
                int cMax = nOrder - a - b;

                for ( int c=0; c<=cMax; c++ )
                {
                    double sum = 0.0;

                    for ( int q=0; q<w.size(); q++ )
                    {
                        sum += 0.5*w[q] * pow( x[q], ( double ) a ) * pow( y[q], ( double ) b )
                               * pow( 1.0 - x[q] - y[q], ( double ) c );
                    }

                    double err = fabs( sum - exact( a,b,c ) );
                    bool localPass = err < 1.0e-14;
                    pass = pass && localPass;

                    if ( !localPass )
                    {
                        fprintf( stderr, "order=%d m (%d, %d, %d) q=%22.15g exact=%22.15g\n", nOrder, a, b, c, sum, exact( a, b, c ) );
                        std::cerr << "error = " << err << std::endl;
                    }

                    else
                    {
                        fprintf( stderr, "order=%d m (%d, %d, %d) q=%22.15g exact=%22.15g\n", nOrder, a, b, c, sum, exact( a, b, c ) );
                    }
                }
            }
        }

        return pass;
    }
    bool test( mpl::int_<3> )
    {
        ublas::vector<value_type> x( ublas::row( this->M_points, 0 ) );
        ublas::vector<value_type> y( ublas::row( this->M_points, 1 ) );
        ublas::vector<value_type> z( ublas::row( this->M_points, 2 ) );
        ublas::vector<value_type> w( this->M_w );

        bool pass = true;
        int p = nOrder;

        for ( int a=0; a<=p; a++ )
        {
            for ( int b=0; b<p-a; b++ )
            {
                int cMax = p - a - b;

                for ( int c=0; c<=cMax; c++ )
                {
                    int dMax = p - a - b - c;

                    for ( int d=0; d<=dMax; d++ )
                    {
                        double sum = 0.0;

                        for ( int q=0; q<w.size(); q++ )
                        {
                            sum += ( 1.0/6.0 )*w[q] * pow( x[q], ( double ) a ) * pow( y[q], ( double ) b )
                                   * pow( z[q], ( double ) c )
                                   * pow( 1.0 - x[q] - y[q] - z[q], ( double ) d );
                        }

                        double err = fabs( sum - exact( a,b,c,d ) );
                        bool localPass = err < 1.0e-14;
                        pass = pass && localPass;

                        if ( !localPass )
                        {
                            fprintf( stderr, "order=%d m (%d, %d, %d %d) q=%22.15g exact=%22.15g\n", p, a, b, c, d, sum, exact( a, b, c, d ) );
                            std::cerr << "error = " << err << std::endl;
                        }

                        else
                        {
                            fprintf( stderr, "order=%d m (%d, %d, %d %d) q=%22.15g exact=%22.15g\n", p, a, b, c, d, sum, exact( a, b, c, d ) );
                        }
                    }
                }
            }
        }

        return pass;

    }

    double exact( int a, int b, int c ) const
    {
        return fact( a )*fact( b )*fact( c )/fact( a+b+c+2 );
    }
    double exact( int a, int b, int c, int d ) const
    {
        return fact( a )*fact( b )*fact( c )*fact( d )/fact( a+b+c+d+3 );
    }

    double fact( int x ) const
    {
        if ( x==0 ) return 1.0;

        return x*fact( x-1 );
    }
    //@}


private:
    quad_type M_quad;
};

#if 0

/*******      Particular Quadrature Rule from old Feel       *******/

/** 7 points Integration rule for triangle (Ref. Stroud) D of Ex = 5 **/

template< typename T >
class TriangleQuadRule
{
public:

    typedef T value_type;
    typedef typename node<value_type>::type node_type;
    typedef ublas::matrix<value_type, ublas::column_major> nodes_type;
    typedef ublas::vector<value_type> weights_type;

    nodes_type const& points() const
    {
        return _pts;
    }

    ublas::matrix_column<nodes_type const> point( uint32_type __i ) const
    {
        return ublas::column( _pts, __i );
    }

    weights_type const& weights() const
    {
        return M_w;
    }
    value_type const& weight( int q ) const
    {
        return M_w[q];
    }

    value_type integrate( boost::function<value_type( node_type const& )> const& f ) const
    {
        value_type res = 0.0;

        for ( uint16_type k = 0; k < 7 ; ++k )
        {
            res += this->weights()[k]*f( this->point( k ) );
        }

        return res;
    }

    TriangleQuadRule()
        : _pts( 2, 7 ), M_w( 7 )
    {
        value_type t7pt_x[3] = {value_type( 1.0 )/3.0, 0.10128650732345633, 0.47014206410511508};
        value_type t7pt_w[3] = {0.1125, 0.062969590272413576, 0.066197076394253090};

        for ( int i = 0; i < 3 ; ++i )
        {
            t7pt_x[i]= 2.0*t7pt_x[i]-1.0;
            t7pt_w[i]= 4.0*t7pt_w[i];
        }

        M_w( 0 ) = t7pt_w[0];
        M_w( 1 ) = t7pt_w[1];
        M_w( 2 ) = t7pt_w[1];
        M_w( 3 ) = t7pt_w[1];
        M_w( 4 ) = t7pt_w[2];
        M_w( 5 ) = t7pt_w[2];
        M_w( 6 ) = t7pt_w[2];

        _pts( 0, 0 ) = t7pt_x[0];
        _pts( 1, 0 ) = t7pt_x[0];

        _pts( 0, 1 ) = t7pt_x[1];
        _pts( 1, 1 ) = t7pt_x[1];

        _pts( 0, 2 ) = t7pt_x[1];
        _pts( 1, 2 ) = 1.0-2.0*t7pt_x[1];

        _pts( 0, 3 ) = 1.0-2.0*t7pt_x[1];
        _pts( 1, 3 ) = t7pt_x[1];

        _pts( 0, 4 ) = t7pt_x[2];
        _pts( 1, 4 ) = t7pt_x[2];

        _pts( 0, 5 ) = t7pt_x[2];
        _pts( 1, 5 ) = 1.0-2.0*t7pt_x[2];

        _pts( 0, 6 ) = 1.0-2.0*t7pt_x[2];
        _pts( 1, 6 ) = t7pt_x[2];
    }

    ~TriangleQuadRule() {}

protected:
    nodes_type _pts;
    weights_type M_w;
};


#endif // 0
}

#endif /* __IMSimplex_H */
