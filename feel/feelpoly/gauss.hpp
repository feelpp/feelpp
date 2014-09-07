/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-12-30

  Copyright (C) 2006 Universite Joseph Fourier (Grenoble)

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
   \file gauss.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-12-30
 */
#ifndef __Gauss_H
#define __Gauss_H 1

#include <feel/feelpoly/quadpoint.hpp>

namespace Feel
{
/*!
 * \class Gauss
 * \brief Gauss quadrature points
 *
 * \code
 * // generate a Gauss point set that would integrate exactly linear
 * // functions using double precision numerical type
 * Gauss<Simplex,1,double> gauss;
 * \endcode
 *
 * \ingroup Polynomial
 * @author Gilles Steiner
 * @author Christophe Prud'homme
 */
template<class Convex, uint16_type Integration_Degree, typename T>
class Gauss : public PointSetQuadrature<Convex, Integration_Degree, T>  {};

template< uint16_type Integration_Degree, typename T>
class Gauss<Simplex<0,1> , Integration_Degree ,T >  : public PointSetQuadrature<Simplex<0,1> , Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Simplex<0,1> , Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    static const uint16_type Degree = invalid_uint16_type_value;
    static const uint32_type Npoints = 1;

    Gauss()
        :
        super( Npoints )
    {

    }

    ~Gauss() {}

    FEELPP_DEFINE_VISITABLE();
};

/// \cond detail
template< uint16_type Integration_Degree, typename T>
class Gauss<Simplex<1,1> , Integration_Degree ,T >  : public PointSetQuadrature<Simplex<1,1> , Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Simplex<1,1> , Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    static const uint16_type Degree = ( Integration_Degree+1 )/2+1;
    static const uint32_type Npoints = Degree;

    typedef Gauss<Simplex<0,1>,Integration_Degree, T> face_quad_type;

    Gauss()
        :
        super( Npoints )
    {
        ublas::vector<T> px( Npoints );

        details::gaussjacobi<Npoints, T, ublas::vector<T>, ublas::vector<T> >( this->M_w, px );
        ublas::row( this->M_points, 0 ) = px;

        boost::shared_ptr<GT_Lagrange<1,1,1,Simplex,T> > gm( new GT_Lagrange<1, 1, 1, Simplex, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Simplex<1, 1>, 1, 1>(), gm, face_qr );


    }

    ~Gauss() {}

    FEELPP_DEFINE_VISITABLE();
};

/** Gauss Quadrature on a triangle **/

template< uint16_type Integration_Degree, typename T>
class Gauss<Simplex<2,1> , Integration_Degree ,T >  : public PointSetQuadrature<Simplex<2,1> , Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Simplex<2,1> , Integration_Degree, T> super;
    typedef typename super::return_type return_type;

    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    typedef Gauss<Simplex<1,1>,Integration_Degree, T> face_quad_type;


    static const uint16_type Degree = ( Integration_Degree+1 )/2+1;
    static const uint32_type Npoints = Degree*Degree;

    Gauss()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Degree );
        weights_type px( Degree );
        details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );

        weights_type wy( Degree );
        weights_type py( Degree );
        details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wy, py, 1.0, 0.0 );

        // coordinate in cartesian space

#if 0
        std::cout<<"[Debug quadpoint] Npoints = " << Npoints << std::endl ;
        std::cout<<"[Debug quadpoint] _pts.size2() = " << this->M_points.size2() << std::endl ;
        std::cout<<"[Debug quadpoint] Degree = " << Degree << std::endl ;
#endif

        node_type eta( 2 );
        details::xi<TRIANGLE, value_type> to_xi;

        for ( int i = 0,  k = 0; i < Degree; ++i )
        {
            for ( int j = 0; j < Degree; ++j, ++k )
            {
#if 0

                if ( j%100==0 )
                    std::cout<<"[Debug quadpoint] i = " << i << " ; j = " << j << " ; k = " << k << std::endl ;

#endif

                // computes the weight of the k-th node
                this->M_w( k ) = 0.5 * wx( i ) * wy( j );
                // use expansion for the collapsed triangle to compute the points
                // coordinates (from cartesian to collapsed coordinates)
                eta( 0 ) = px( i );
                eta( 1 ) = py( j );

                ublas::column( this->M_points, k ) = to_xi( eta );
            }
        }

        boost::shared_ptr<GT_Lagrange<2,1,2,Simplex,T> > gm( new GT_Lagrange<2, 1, 2, Simplex, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Simplex<2, 1>, 2, 1>(), gm, face_qr );

    }

    ~Gauss() {}

    FEELPP_DEFINE_VISITABLE();
};


/** Gauss Quadrature on a tetrahedra **/

template< uint16_type Integration_Degree, typename T>
class Gauss<Simplex<3,1> , Integration_Degree ,T >  : public PointSetQuadrature<Simplex<3,1> , Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Simplex<3,1> , Integration_Degree, T> super;

    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    typedef Gauss<Simplex<2,1>,Integration_Degree, T> face_quad_type;


    static const uint16_type Degree = ( Integration_Degree+1 )/2+1;
    static const uint32_type Npoints = Degree*Degree*Degree;

    Gauss()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Degree );
        weights_type px( Degree );
        details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );

        weights_type wy( Degree );
        weights_type py( Degree );
        details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wy, py, 1.0, 0.0 );

        weights_type wz( Degree );
        weights_type pz( Degree );
        details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wz, pz, 2.0, 0.0 );

        // coordinate in cartesian space
        node_type eta( 3 );
        details::xi<TETRAHEDRON, value_type> to_xi;

        for ( int i = 0,  k = 0; i < Degree; ++i )
        {
            for ( int j = 0; j < Degree; ++j )
            {
                for ( int l = 0; l < Degree; ++l, ++k )
                {
                    // computes the weight of the k-th node
                    this->M_w( k ) = 0.125 * wx( i ) * wy( j ) * wz( l );
                    // use expansion for the collapsed triangle to compute the points
                    // coordinates (from cartesian to collapsed coordinates)
                    eta( 0 ) = px( i );
                    eta( 1 ) = py( j );
                    eta( 2 ) = pz( l );
                    ublas::column( this->M_points, k ) = to_xi( eta );
                }
            }
        }

        boost::shared_ptr<GT_Lagrange<3, 1, 3, Simplex, T> > gm( new GT_Lagrange<3, 1, 3, Simplex, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Simplex<3, 1>,3,1>(), gm, face_qr );
    }

    ~Gauss() {}

    FEELPP_DEFINE_VISITABLE();
};

/** Gauss Quadrature on Simplex Product **/

template< uint16_type Integration_Degree, typename T>
class Gauss<Hypercube<1,1>, Integration_Degree ,T >
    :
public PointSetQuadrature<Hypercube<1,1>, Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Hypercube<1,1>, Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    static const uint16_type Degree = ( Integration_Degree+1 )/2+1;
    static const uint32_type Npoints = Degree;

    typedef Gauss<Simplex<0,1>,Integration_Degree, T> face_quad_type;

    Gauss(  )
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Degree );
        weights_type px( Degree );
        details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );
#if 0
        VLOG(1) << "[gauss<SP<2,1>] jacobi p = " << px << "\n";
        VLOG(1) << "[gauss<SP<2,1>] jacobi w = " << wx << "\n";
#endif

        for ( int i = 0; i < Degree; ++i )
        {
            // computes the weight of the k-th node
            this->M_w( i ) = wx( i );
            this->M_points( 0, i ) = px( i );
        }


#if 0
        VLOG(1) << "[gauss<SP<2,1>] p = " << this->M_points << "\n";
        VLOG(1) << "[gauss<SP<2,1>] w = " << this->M_w << "\n";
#endif


        boost::shared_ptr<GT_Lagrange<1,1,1, Hypercube,T> > gm( new GT_Lagrange<1, 1, 1, Hypercube, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Hypercube<1, 1>, 1, 1>(), gm, face_qr );

    }

    ~Gauss() {}

    FEELPP_DEFINE_VISITABLE();
};
/** Gauss Quadrature on the quadrangle [-1,1]x[-1,1] **/

template< uint16_type Integration_Degree, typename T>
class Gauss<Hypercube<2,1>, Integration_Degree ,T >
    :
public PointSetQuadrature<Hypercube<2,1>, Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Hypercube<2,1>, Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    typedef Gauss<Hypercube<1,1>,Integration_Degree, T> face_quad_type;

    static const uint16_type Degree = ( Integration_Degree+1 )/2+1;
    static const uint32_type Npoints = Degree*Degree;

    Gauss()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Degree );
        weights_type px( Degree );
        details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );
#if 0
        VLOG(1) << "[gauss<SP<2,1>] jacobi p = " << px << "\n";
        VLOG(1) << "[gauss<SP<2,1>] jacobi w = " << wx << "\n";
#endif

        for ( int i = 0,  k = 0; i < Degree; ++i )
        {
            for ( int j = 0; j < Degree; ++j, ++k )
            {
                // computes the weight of the k-th node
                this->M_w( k ) = wx( i ) * wx( j );
                this->M_points( 0, k ) = px( i );
                this->M_points( 1, k ) = px( j );
            }
        }

#if 0
        VLOG(1) << "[gauss<SP<2,1>] p = " << this->M_points << "\n";
        VLOG(1) << "[gauss<SP<2,1>] w = " << this->M_w << "\n";
#endif
        boost::shared_ptr<GT_Lagrange<2, 1, 2, Hypercube, T> > gm( new GT_Lagrange<2, 1, 2, Hypercube, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Hypercube<2, 1>,2,1>(), gm, face_qr );

    }

    ~Gauss() {}

    FEELPP_DEFINE_VISITABLE();
};

/** Gauss Quadrature on the hexahedra [-1,1]x[-1,1]x[-1,1] **/

template< uint16_type Integration_Degree, typename T>
class Gauss<Hypercube<3,1>, Integration_Degree ,T >
    :
public PointSetQuadrature<Hypercube<3,1>, Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Hypercube<3,1>, Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;
    typedef Gauss<Hypercube<2,1>,Integration_Degree, T> face_quad_type;
    static const uint16_type Degree = ( Integration_Degree+1 )/2+1;
    static const uint32_type Npoints = Degree*Degree*Degree;

    Gauss()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Degree );
        weights_type px( Degree );
        details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );

        for ( int i = 0,  k = 0; i < Degree; ++i )
        {
            for ( int j = 0; j < Degree; ++j )
            {
                for ( int l = 0; l < Degree ; ++l, ++k )
                {
                    // computes the weight of the k-th node
                    this->M_w( k ) = wx( i ) * wx( j ) * wx( l );
                    this->M_points( 0, k ) = px( i );
                    this->M_points( 1, k ) = px( j );
                    this->M_points( 2, k ) = px( l );
                }
            }
        }

        boost::shared_ptr<GT_Lagrange<3, 1, 3, Hypercube, T> > gm( new GT_Lagrange<3, 1, 3, Hypercube, T> );
        boost::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Hypercube<3, 1>,3,1>(), gm, face_qr );

    }

    ~Gauss() {}

    FEELPP_DEFINE_VISITABLE();
};


template< uint16_type Integration_Degree, typename T>
class Gauss<Hypercube<4,1>, Integration_Degree ,T >
    :
public PointSetQuadrature<Hypercube<4,1>, Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Hypercube<4,1>, Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    static const uint16_type Degree = ( Integration_Degree+1 )/2+1;
    static const uint32_type Npoints = Degree*Degree*Degree*Degree;

    Gauss()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Degree );
        weights_type px( Degree );
        details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );

        for ( int i = 0,  k = 0; i < Degree; ++i )
        {
            for ( int j = 0; j < Degree; ++j )
            {
                for ( int l = 0; l < Degree ; ++l )
                {
                    for ( int r = 0; r < Degree ; ++r, ++k )
                    {
                        // computes the weight of the k-th node
                        this->M_w( k ) = wx( i ) * wx( j ) * wx( l ) * wx( r );
                        this->M_points( 0, k ) = px( i );
                        this->M_points( 1, k ) = px( j );
                        this->M_points( 2, k ) = px( l );
                        this->M_points( 3, k ) = px( r );
                    }
                }
            }
        }
    }

    ~Gauss() {}

    FEELPP_DEFINE_VISITABLE();
};


template< uint16_type Integration_Degree, typename T>
class Gauss<Hypercube<5,1>, Integration_Degree ,T >
    :
public PointSetQuadrature<Hypercube<5,1>, Integration_Degree, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Hypercube<5,1>, Integration_Degree, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    static const uint16_type Degree = ( Integration_Degree+1 )/2+1;
    static const uint32_type Npoints = Degree*Degree*Degree*Degree*Degree;

    Gauss()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Degree );
        weights_type px( Degree );
        details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );

        for ( int i = 0,  k = 0; i < Degree; ++i )
        {
            for ( int j = 0; j < Degree; ++j )
            {
                for ( int l = 0; l < Degree ; ++l )
                {
                    for ( int r = 0; r < Degree ; ++r )
                    {
                        for ( int s = 0; s < Degree ; ++s, ++k )
                        {
                            // computes the weight of the k-th node
                            this->M_w( k ) = wx( i ) * wx( j ) * wx( l ) * wx( r ) * wx ( s );
                            this->M_points( 0, k ) = px( i );
                            this->M_points( 1, k ) = px( j );
                            this->M_points( 2, k ) = px( l );
                            this->M_points( 3, k ) = px( r );
                            this->M_points( 4, k ) = px( s );
                        }
                    }
                }
            }
        }
    }

    ~Gauss() {}

    FEELPP_DEFINE_VISITABLE();
};
/// \endcond
} // Feel

#endif /* __Gauss_H */
