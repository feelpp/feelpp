/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Thomas Lantz <thomas.lantz@etu.unistra.fr>
       Date: 2015-03-09
@
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
   \file multiscalequadrature.hpp
   \author Thomas Lantz 
   \date 2015-04-27
 */
#ifndef __TestMultiScale_H
#define __TestMultiScale_H 1

#include <feel/feelpoly/pointsetquadrature.hpp>



namespace Feel
{
/*!
 * \class MultiScaleQuadrature
 * \brief quadrature points by rectangle approximation for multi scale approach
 *
 * \code
 * // generate a  point set that would integrate linear
 * // functions ( and images ) using double precision numerical type
 * MultiScaleQuadrature<Simplex,1,double> msq;
 * \endcode
 *
 * @author Thomas Lantz
 */

template<class Convex, int Integration_Degree, typename T>
class MultiScaleQuadrature : public PointSetQuadrature<Convex, T>  {};

template< int Integration_Degree, typename T>
class MultiScaleQuadrature<Simplex<0,1> , Integration_Degree ,T >  : public PointSetQuadrature<Simplex<0,1>,  T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Simplex<0,1> , T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    static const int Degree = invalid_uint16_type_value;
    static const uint32_type Npoints = 1;

    MultiScaleQuadrature()
        :
        super( Npoints )
    {

    }

    ~MultiScaleQuadrature() {}

    FEELPP_DEFINE_VISITABLE();
};

/** Multi Scale Quadrature on the segment [-1,1] **/

template< int Integration_Degree, typename T>
class MultiScaleQuadrature<Hypercube<1,1>, Integration_Degree ,T >
    :
public PointSetQuadrature<Hypercube<1,1>, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Hypercube<1,1>, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    static const int Degree = 1;
    static const uint32_type Npoints = Integration_Degree;

    typedef MultiScaleQuadrature<Simplex<0,1>,Integration_Degree, T> face_quad_type;


    MultiScaleQuadrature()
        :
        super( std::pow(2,ioption("msi.level"))+1  )
    {

        int  gridsize = std::pow(2,ioption("msi.level"));

        // build rules in x and y direction
        weights_type wx( gridsize+1  );
        weights_type px( gridsize+1  );

        double tmp=-1;
        for ( int i = 0; i <=gridsize ; i++ )
        {
            // computes the weight of the k-th node
            this->M_w( i ) = 2./(gridsize+1) ;// wx( i );
            this->M_points( 0, i ) = tmp ;
            tmp+=2./gridsize ;
        }

        std::shared_ptr<GT_Lagrange<1,1,1, Hypercube,T> > gm( new GT_Lagrange<1, 1, 1, Hypercube, T> );
        std::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Hypercube<1, 1>, 1, 1>(), gm, face_qr );

    }
/*
    MultiScaleQuadrature()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Npoints );
        weights_type px( Npoints );
        //details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );
#if 0
        VLOG(1) << "[gauss<SP<2,1>] jacobi p = " << px << "\n";
        VLOG(1) << "[gauss<SP<2,1>] jacobi w = " << wx << "\n";
#endif
        double tmp=-1;
        for ( int i = 0; i < Npoints; i++ )
        {
            // computes the weight of the k-th node
            this->M_w( i ) = 2./Npoints ;// wx( i );
            this->M_points( 0, i ) = tmp  ;
            tmp+=2./Npoints;
        }


#if 0
        VLOG(1) << "[gauss<SP<2,1>] p = " << this->M_points << "\n";
        VLOG(1) << "[gauss<SP<2,1>] w = " << this->M_w << "\n";
#endif


        std::shared_ptr<GT_Lagrange<1,1,1, Hypercube,T> > gm( new GT_Lagrange<1, 1, 1, Hypercube, T> );
        std::shared_ptr<face_quad_type> face_qr( new face_quad_type );
        // construct face quadratures
        this->constructQROnFace( Reference<Hypercube<1, 1>, 1, 1>(), gm, face_qr );

    }
*/

    ~MultiScaleQuadrature() {}

    FEELPP_DEFINE_VISITABLE();
};


/** Multi Scale Quadrature on the quadrangle [-1,1]x[-1,1] **/

template< int Integration_Degree, typename T>
class MultiScaleQuadrature<Hypercube<2,1>, Integration_Degree ,T >
    :
public PointSetQuadrature<Hypercube<2,1>, T>
{
public :
    typedef T value_type;

    typedef PointSetQuadrature<Hypercube<2,1>, T> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef typename super::weights_type weights_type;

    typedef MultiScaleQuadrature<Hypercube<1,1>,Integration_Degree, T> face_quad_type;

    static const int Degree = Integration_Degree;
    static const uint32_type Npoints =Integration_Degree*Integration_Degree;


    MultiScaleQuadrature( )
        :
        super( (1+std::pow(2,ioption("msi.level")))*(std::pow(2,ioption("msi.level"))+1) )
    {

        int gridsize=std::pow(2,ioption("msi.level"));

        // build rules in x and y direction
        weights_type wx( (gridsize+1)*(gridsize+1) );
        //weights_type px( (gridsize*gridsize );
        
        double tmpx=-1.;
        double tmpy=1.;
        for ( int i = 0,  k = 0; i <= gridsize; i++ )
        {
            for ( int j = 0; j <=  gridsize; j++, ++k )
            {
                // computes the weight of the k-th node
                this->M_w( k ) = 4./( (gridsize+1)*(gridsize+1)) ;//wx( i ) * wx( j );
                this->M_points( 0, k ) = tmpx ;
                this->M_points( 1, k ) = tmpy ;
                tmpx+=2./gridsize;
            }
            tmpx=-1.;
            tmpy-=2./gridsize;
        }
        std::cout << "quadrature points:" << this->M_points << std::endl; 
        std::shared_ptr<GT_Lagrange<2, 1, 2, Hypercube, T> > gm( new GT_Lagrange<2, 1, 2, Hypercube, T> );
        std::shared_ptr<face_quad_type> face_qr( new face_quad_type());
        // construct face quadratures
        this->constructQROnFace( Reference<Hypercube<2, 1>,2,1>(), gm, face_qr );

    }

/*

    MultiScaleQuadrature()
        :
        super( Npoints )
    {
        // build rules in x and y direction
        weights_type wx( Npoints );
        weights_type px( Npoints );
        //details::gaussjacobi<Degree,T, ublas::vector<T>, ublas::vector<T> >( wx, px, 0.0, 0.0 );
#if 0
        VLOG(1) << "[gauss<SP<2,1>] jacobi p = " << px << "\n";
        VLOG(1) << "[gauss<SP<2,1>] jacobi w = " << wx << "\n";
#endif
        double tmpx=-1;
        double tmpy=-1;
        for ( int i = 0,  k = 0; i < Degree; i++ )
        {
            for ( int j = 0; j < Degree; j++, ++k )
            {
                // computes the weight of the k-th node
                this->M_w( k ) = 4.*(1./Npoints) ;//wx( i ) * wx( j );
                this->M_points( 0, k ) = tmpx ;
                this->M_points( 1, k ) = tmpy ;
                tmpy+=2./Degree;

            }
            tmpx+=2./Degree;
        }

#if 0
        VLOG(1) << "[gauss<SP<2,1>] p = " << this->M_points << "\n";
        VLOG(1) << "[gauss<SP<2,1>] w = " << this->M_w << "\n";
#endif
        std::shared_ptr<GT_Lagrange<2, 1, 2, Hypercube, T> > gm( new GT_Lagrange<2, 1, 2, Hypercube, T> );
        std::shared_ptr<face_quad_type> face_qr( new face_quad_type);
        // construct face quadratures
        this->constructQROnFace( Reference<Hypercube<2, 1>,2,1>(), gm, face_qr );

    }
    */

    ~MultiScaleQuadrature() {}

    FEELPP_DEFINE_VISITABLE();
};
}
#endif /* __TestMultiScale_H */
