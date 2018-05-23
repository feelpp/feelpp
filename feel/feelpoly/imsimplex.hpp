/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=1
*/

template<typename T>
struct IMTriangle<1,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 1;
    static const uint16_type nPoints = 1;

    IMTriangle();

    std::vector<value_type> q;
};

template<typename T> struct IMTriangle<0,T>: public IMTriangle<1,T> {};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=2
*/
template<typename T>
struct IMTriangle<2,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 2;
    static const uint16_type nPoints = 3;

    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=3
*/

template<typename T>
struct IMTriangle<3,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 3;
    static const uint16_type nPoints = 4;
    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=4
*/

template<typename T>
struct IMTriangle<4,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 4;
    static const uint16_type nPoints = 6;
    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=5
*/
template<typename T>
struct IMTriangle<5,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 5;
    static const uint16_type nPoints = 7;
    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=6
*/
template<typename T>
struct IMTriangle<6,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 6;
    static const uint16_type nPoints = 12;
    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=7
*/
template<typename T>
struct IMTriangle<7,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 7;
    static const uint16_type nPoints = 13;
    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=8
*/
template<typename T>
struct IMTriangle<8,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 8;
    static const uint16_type nPoints = 16;
    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=9
*/
template<typename T>
struct IMTriangle<9,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 9;
    static const uint16_type nPoints = 19;
    IMTriangle();

    std::vector<value_type> q;
};
/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=10
*/
template<typename T>
struct IMTriangle<10,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 10;
    static const uint16_type nPoints = 25;
    IMTriangle();

    std::vector<value_type> q;
};
/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=11
*/

template<typename T>
struct IMTriangle<11,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 11;
    static const uint16_type nPoints = 27;
    IMTriangle();

    std::vector<value_type> q;
};
/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=12
*/
template<typename T>
struct IMTriangle<12,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 12;
    static const uint16_type nPoints = 33;
    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=13
*/
template<typename T>
struct IMTriangle<13,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 13;
    static const uint16_type nPoints = 37;
    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=14
*/
template<typename T>
struct IMTriangle<14,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 14;
    static const uint16_type nPoints = 42;
    IMTriangle();

    std::vector<value_type> q;
};
/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=15
*/
template<typename T>
struct IMTriangle<15,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 15;
    static const uint16_type nPoints = 48;
    IMTriangle();

    std::vector<value_type> q;
};
/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=16
*/
template<typename T>
struct IMTriangle<16,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 16;
    static const uint16_type nPoints = 52;
    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=17
*/
template<typename T>
struct IMTriangle<17,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 17;
    static const uint16_type nPoints = 61;
    IMTriangle();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=18
*/
template<typename T>
struct IMTriangle<18,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 18;
    static const uint16_type nPoints = 70;
    IMTriangle();

    std::vector<value_type> q;
};
/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=19
*/
template<typename T>
struct IMTriangle<19,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 19;
    static const uint16_type nPoints = 73;
    IMTriangle();

    std::vector<value_type> q;
};
/*
  Gauss  quadrature  points  and  weights  on  the  reference  triangle  order  p=20
 */
template<typename T>
struct IMTriangle<20,T>
{
    typedef T value_type;
    static const uint16_type nDim = 2;
    static const uint16_type nOrder = 20;
    static const uint16_type nPoints = 79;
    IMTriangle();

    std::vector<value_type> q;
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

/**
 Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=1
*/
template<typename T>
struct IMTetrahedra<1,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 1;
    static const uint16_type nPoints = 1;

    IMTetrahedra();

    std::vector<value_type> q;
};
template<typename T> struct IMTetrahedra<0,T> : public IMTetrahedra<1,T> {};
/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=2
 */
template<typename T>
struct IMTetrahedra<2,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 2;
    static const uint16_type nPoints = 4;
    IMTetrahedra();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=3
 */
template<typename T>
struct IMTetrahedra<3,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 3;
    static const uint16_type nPoints = 5;
    IMTetrahedra();

    std::vector<value_type> q;
};
/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=4
*/
template<typename T> struct IMTetrahedra<4,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 4;
    static const uint16_type nPoints = 11;
    IMTetrahedra();

    std::vector<value_type> q;
};
/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=5
*/
template<typename T> struct IMTetrahedra<5,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 5;
    static const uint16_type nPoints = 14;
    IMTetrahedra();

    std::vector<value_type> q;
};
/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=6
 */
template<typename T>
struct IMTetrahedra<6,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 6;
    static const uint16_type nPoints = 24;
    IMTetrahedra();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=7
*/
template<typename T>
struct IMTetrahedra<7,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 7;
    static const uint16_type nPoints = 31;
    IMTetrahedra();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=8
*/
template<typename T>
struct IMTetrahedra<8,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 8;
    static const uint16_type nPoints = 43;
    IMTetrahedra();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=9
*/
template<typename T>
struct IMTetrahedra<9,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 9;
    static const uint16_type nPoints = 53;
    IMTetrahedra();

    std::vector<value_type> q;
};

/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=11
*/
template<typename T>
struct IMTetrahedra<11,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 11;
    static const uint16_type nPoints = 126;
    IMTetrahedra();

    std::vector<value_type> q;
};
template<typename T> struct IMTetrahedra<10,T>: public IMTetrahedra<11,T> {};
/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=13
*/
template<typename T>
struct IMTetrahedra<13,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 13;
    static const uint16_type nPoints = 210;
    IMTetrahedra();

    std::vector<value_type> q;
};
template<typename T> struct IMTetrahedra<12,T>: public IMTetrahedra<13,T> {};
/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=15
*/
template<typename T>
struct IMTetrahedra<15,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 15;
    static const uint16_type nPoints = 330;
    IMTetrahedra();

    std::vector<value_type> q;
};
template<typename T> struct IMTetrahedra<14,T>: public IMTetrahedra<15,T> {};
/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=17
*/
template<typename T>
struct IMTetrahedra<17,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 17;
    static const uint16_type nPoints = 495;
    IMTetrahedra();

    std::vector<value_type> q;
};
template<typename T> struct IMTetrahedra<16,T>: public IMTetrahedra<17,T> {};
/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=19
*/
template<typename T>
struct IMTetrahedra<19,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 19;
    static const uint16_type nPoints = 715;
    IMTetrahedra();

    std::vector<value_type> q;
};
template<typename T> struct IMTetrahedra<18,T>: public IMTetrahedra<19,T> {};
/*
  Gauss  quadrature  constants  for  the  reference  tetrahedron  order  p=21
*/
template<typename T>
struct IMTetrahedra<21,T>
{
    typedef T value_type;
    static const uint16_type nDim = 3;
    static const uint16_type nOrder = 21;
    static const uint16_type nPoints = 1001;
    IMTetrahedra();

    std::vector<value_type> q;
};
template<typename T> struct IMTetrahedra<20,T>: public IMTetrahedra<21,T> {};
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
    typedef Simplex<Dim,1> convex_type;
    typedef T value_type;
    typedef ublas::matrix<value_type,ublas::column_major> matrix_type;
    typedef ublas::vector<value_type> vector_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<Dim>,mpl::int_<2> >,
            mpl::identity<Feel::detail::IMTriangle<Order,T> >,
            mpl::identity<Feel::detail::IMTetrahedra<Order,T> > >::type::type quad_type;

#if 0
#if 0
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<Dim>,mpl::int_<2> >,
            mpl::identity<Gauss<Simplex<Dim-1,1>,Order,T> >,
            mpl::identity<Gauss<Simplex<Dim-1,1>,Order,T> > >::type::type face_quad_type;
    //mpl::identity<IMSimplex<Dim-1,Order,T> > >::type::type face_quad_type;
#else
    typedef Gauss<Simplex<Dim-1,1>,Order,T> face_quad_type;
#endif
#else
    using face_quad_type =  typename mpl::if_<mpl::equal_to<mpl::int_<Dim>,mpl::int_<3>>, mpl::identity<IMSimplex<2,Order,T> >, mpl::identity<Gauss<Simplex<Dim-1,1>,Order,T> > >::type::type;
    //typedef IMSimplex<((Dim==1) || (Dim ==0))?0:Dim-1,Order,T> face_quad_type;
#endif
    typedef IMSimplex<Dim,Order,T> parent_quadrature_type;
    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;
    static const uint16_type nQuadPoints = quad_type::nPoints;

    //@}

    /** @name Constructors, destructor
     */
    //@{
    IMSimplex(  )
        :
        super( quad_type::nPoints ),
        M_quad()
    {
        for ( size_type i=0; i< quad_type::nPoints; i++ )
        {

            for ( int j = 0; j < Dim; ++j )
            {
                this->M_points( j, i ) = M_quad.q[( Dim+1 )*i+j];
            }

            this->M_w( i ) = M_quad.q[( Dim+1 )*i+Dim];
        }

        boost::shared_ptr<GT_Lagrange<Dim,1,Dim,Simplex,T> > gm( new GT_Lagrange<Dim, 1, Dim, Simplex, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Simplex<Dim, 1>, Dim, 1>(), gm, face_qr );
    }

    ~IMSimplex()
    {}

    //@}

    /**
     **
     */
    T factor() const
    {
        return ( Dim==2 )?T( 1 )/T( 2 ):T( 1 )/T( 6 );
    }


    /** @name  Methods
     */
    //@{
    IMSimplex& operator=( IMSimplex const & i )
    {
        if ( this != &i )
        {
            M_quad = i.M_quad;
        }

        return *this;
    }
    //@}


private:
    quad_type M_quad;
};

}

#endif /* __IMSimplex_H */
