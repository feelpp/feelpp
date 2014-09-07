/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-12-30

  Copyright (C) 2006-2011 Université Joseph Fourier

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
   \file im.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-12-30
 */
#ifndef __im_H
#define __im_H 1
#include <iostream>

#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/glas.hpp>

#include <feel/feelpoly/gauss.hpp>
#include <feel/feelpoly/gausslobatto.hpp>
#include <feel/feelpoly/imsimplex.hpp>
#include <feel/feelpoly/imexact.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

/**
 * degree policy that defines the integration method that used jacobi
 * polynomials of degree N. It computes also the exactness (2*N-1)
 */
template<uint16_type N>
struct JacobiDegree
{
    static const uint16_type jacobi_degree = N;
    static const uint16_type integration_degree = 2*jacobi_degree-1;
};
/**
 * order policy that computes integrals of polynomials up to degree
 * N. It selects the proper jacobi polynomials degree.
 */
template<uint16_type N>
struct IntegrationDegree
{
    static const uint16_type integration_degree = N;
    static const uint16_type jacobi_degree = ( integration_degree+1 )/2+1;
};


/**
 * 
 * \brief integration method interface class
 *
 * \code
 * // integrate exactly linear functions in 2D over triangles using
 * // Gauss points and double precision
 * IM<2,1,double,Simplex,Gauss> im;
 *
 * // integrate exactly polynomials of order 5in 3D over hexahedrons using
 * // GaussLobatto points and double precision
 * IM<3,5,double,Hypercube,GaussLobatto> im;
 *
 *
 * \endcode
 **/
template<int Dim,
         uint16_type Order,
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Entity = Simplex,
         template<class Convex, uint16_type O, typename T2> class QPS = Gauss,
         template<uint16_type N> class DegreePolicy = IntegrationDegree>
class IMGeneral
    :
public QPS<Entity<Dim,1,Dim>, DegreePolicy<Order>::integration_degree, T>
{
    typedef QPS<Entity<Dim,1,Dim>, DegreePolicy<Order>::integration_degree, T> super;
public:
    static const bool is_exact = false;
    typedef DegreePolicy<Order> degree_policy_type;
    static const uint16_type nDim = Dim;
    static const uint16_type nNodes = degree_policy_type::jacobi_degree;
    static const uint16_type nOrder = degree_policy_type::integration_degree;
    static const uint16_type nQuadPoints = super::Npoints;

    //typedef typename super::convex_type convex_type;
    typedef Entity<Dim,1,Dim> convex_type;
    typedef typename super::value_type value_type;
    typedef typename super::node_type node_type;
    typedef typename super::weights_type weights_type;
    typedef typename super::nodes_type nodes_type;

    typedef boost::tuple<nodes_type, weights_type> quadrature_data_type;
    typedef typename super::face_quadrature_type face_quadrature_type;
    typedef IMGeneral<Dim,Order,T,Entity,QPS,DegreePolicy> parent_quadrature_type;

    IMGeneral()
        :
        super()
    {

    }

    ~IMGeneral() {}
    quadrature_data_type data() const
    {
        return boost::make_tuple( this->points(), this->weights() );
    }

};

template<int DIM,
         int IMORDER,
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Entity = Simplex>
struct IM
        :
    public mpl::if_<mpl::and_<mpl::or_<mpl::and_<mpl::less_equal<mpl::int_<IMORDER>,mpl::int_<20> >,
        mpl::equal_to<mpl::int_<DIM>,mpl::int_<2> > >,
        mpl::and_<mpl::less_equal<mpl::int_<IMORDER>,mpl::int_<20> >,
        mpl::equal_to<mpl::int_<DIM>,mpl::int_<3> > > >,
        mpl::bool_<Entity<DIM,1,DIM>::is_simplex> >,
        mpl::identity<IMSimplex<DIM, IMORDER, T> >,
        mpl::identity<IMGeneral<DIM, IMORDER, T, Entity> > >::type::type
{
    template<int DIM1,
             typename T1,
             template<uint16_type, uint16_type, uint16_type> class Entity1>
    struct apply
    {
        typedef typename mpl::if_<mpl::and_<mpl::or_<mpl::and_<mpl::less_equal<mpl::int_<IMORDER>,mpl::int_<20> >,
                mpl::equal_to<mpl::int_<DIM1>,mpl::int_<2> > >,
                mpl::and_<mpl::less_equal<mpl::int_<IMORDER>,mpl::int_<20> >,
                mpl::equal_to<mpl::int_<DIM1>,mpl::int_<3> > > >,
                mpl::bool_<Entity1<DIM1,1,DIM1>::is_simplex> >,
                mpl::identity<IMSimplex<DIM1, IMORDER, T1> >,
                mpl::identity<IMGeneral<DIM1, IMORDER, T1, Entity1> > >::type::type type;
    };
};

template<int IMORDER>
struct _Q
{
    static const int order = IMORDER;

    template<int DIM,
             typename T,
             template<uint16_type, uint16_type, uint16_type> class Entity>
    struct apply
    {
        typedef typename mpl::if_<mpl::and_<mpl::or_<mpl::and_<mpl::less_equal<mpl::int_<IMORDER>,mpl::int_<20> >,
                mpl::equal_to<mpl::int_<DIM>,mpl::int_<2> > >,
                mpl::and_<mpl::less_equal<mpl::int_<IMORDER>,mpl::int_<20> >,
                mpl::equal_to<mpl::int_<DIM>,mpl::int_<3> > > >,
                mpl::bool_<Entity<DIM,1,DIM>::is_simplex> >,
                mpl::identity<IMSimplex<DIM, IMORDER, T> >,
                mpl::identity<IMGeneral<DIM, IMORDER, T, Entity> > >::type::type type;
    };

    template<int DIM,
             typename T,
             template<uint16_type, uint16_type, uint16_type> class Entity>
    struct applyIMGeneral
    {
        typedef IMGeneral<DIM, IMORDER, T, Entity> type;
    };

    template<typename ContextType>
    struct applyContext
    {
        static const int DIM = ContextType::PDim;
        typedef typename ContextType::value_type T;
        typedef typename mpl::if_<mpl::and_<mpl::or_<mpl::and_<mpl::less_equal<mpl::int_<IMORDER>,mpl::int_<20> >,
                mpl::equal_to<mpl::int_<DIM>,mpl::int_<2> > >,
                mpl::and_<mpl::less_equal<mpl::int_<IMORDER>,mpl::int_<20> >,
                mpl::equal_to<mpl::int_<DIM>,mpl::int_<3> > > >,
                mpl::bool_<ContextType::element_type::is_simplex> >,
                mpl::identity<IMSimplex<DIM, IMORDER, T> >,
                typename mpl::if_<mpl::bool_<ContextType::element_type::is_simplex>,
                mpl::identity<IMGeneral<DIM, IMORDER, T, Simplex> >,
                mpl::identity<IMGeneral<DIM, IMORDER, T, Hypercube> > >::type>::type::type type;
    };
};


#if 0
template<int Dim,
         uint16_type Order,
         typename T,
         template<uint16_type,uint16_type,uint16_type> class Entity,
         template<class Convex, uint16_type O, typename T2> class QPS,
         template<uint16_type N> class DegreePolicy>
std::ostream&
operator<<( std::ostream& os,
            Feel::IM<Dim, Order, T, Entity,QPS, DegreePolicy> const& qps )
{
    os << "quadrature point set:\n"
       << "number of points: " << qps.nPoints() << "\n"
       << "points :  " << qps.points() << "\n"
       << "weights :  " << qps.weights() << "\n";

    for ( size_type f = 0; f < qps.nFaces(); ++f )
    {
        os << " o Face " << f << "\n";
        os << "   number of points: " << qps.nPointsOnFace( f ) << "\n"
           << "   points :  " << qps.fpoints( f ) << "\n"
           << "   weights :  " << qps.weights( f ) << "\n";
    }

    return os;
}

template<int Dim, int IMORDER, typename T> struct ImBestSimplex
        :
    public mpl::if_<mpl::less_equal<mpl::int_<IMORDER>,mpl::int_<6> >,
        mpl::identity<IMSimplex<Dim, IMORDER, T> >,
        mpl::identity<IM<Dim, IMORDER, T, Simplex> > >::type::type
{};
#endif
} // Feel


#endif /* __im_H */
