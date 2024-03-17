/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
struct _QBase {};
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
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Entity = Simplex>
class IMGeneral
    :
        public PointSetQuadrature<Entity<Dim,1,Dim> , T, index_type>

{
    using super = PointSetQuadrature<Entity<Dim,1,Dim> , T, index_type>;
    
public:
    static const bool is_exact = false;
    static const uint16_type nDim = Dim;
    static const uint16_type nRealDim = Dim;

    typedef Entity<Dim,1,Dim> convex_type;
    typedef typename super::value_type value_type;
    typedef typename super::node_type node_type;
    typedef typename super::weights_type weights_type;
    typedef typename super::nodes_type nodes_type;

    typedef boost::tuple<nodes_type, weights_type> quadrature_data_type;
    typedef typename super::face_quadrature_type face_quadrature_type;
    typedef IMGeneral<Dim,T,Entity> parent_quadrature_type;
    static constexpr uint16_type Dim_m_1 = (Dim > 1 )?Dim-1:0;
    using face_quad_type = IMGeneral<Dim_m_1,T, Entity>;

    IMGeneral() = default;
    
    explicit IMGeneral( uint16_type o ): super(o) 
        {
            if ( nDim > 0 )
            {
                auto gm = makeGeometricTransformation<convex_type,T>();
                auto face_qr = std::make_shared<face_quad_type>(o);
                this->constructQROnFace( makeReferenceConvex<convex_type,nDim,1,nRealDim>(), gm, face_qr );
            }
        }
    ~IMGeneral() override = default;

    IMGeneral( IMGeneral const& i ) = default;
    IMGeneral( IMGeneral && i ) = default;
    IMGeneral& operator=( IMGeneral const& ) = default;
    IMGeneral& operator=( IMGeneral && ) = default;
    
    quadrature_data_type data() const
    {
        return boost::make_tuple( this->points(), this->weights() );
    }

    void create( uint16_type order ) override
        {
            super::create(order);
            auto gm = makeGeometricTransformation<convex_type,T>();
            auto face_qr = std::make_shared<face_quad_type>(order);
         
            this->constructQROnFace( Reference<convex_type,nDim,1,nRealDim>(), gm, face_qr );
        }
};

//!
//! type of quadrature
//!
template<typename ConvexT, typename T = double>
using im_t = typename mpl::if_<mpl::bool_<ConvexT::is_simplex>,
                               mpl::identity<IMGeneral<ConvexT::nDim, T, Simplex> >,
                               mpl::identity<IMGeneral<ConvexT::nDim, T, Hypercube> > >::type::type;

//!
//! instantiate a quadrature on convex @p ConvexT to integrate up to order @p O
//!
template<typename ConvexT,typename T = double>
im_t<ConvexT,T> im( uint16_type O )
{
    return im_t<ConvexT,T>{ O };
}

template<typename IMT>
IMT im( uint16_type O, std::enable_if_t<std::is_base_of<PointSetBase,IMT>::value>* = nullptr )
{
    CHECK( O != invalid_uint16_type_value ) << "Invalid integration order";
    return IMT{ O };
}

template<typename IMT, typename QT>
IMT im( QT the_q, std::enable_if_t<std::is_base_of<_QBase,QT>::value>* = nullptr )
{
    return IMT{ the_q.order() };
}

template<typename IMT>
IMT&& im( IMT && the_im, std::enable_if_t<std::is_base_of<PointSetBase,IMT>::value>* = nullptr )
{
    return std::forward<IMT>( the_im );
}

template<typename IMT>
IMT const& im( IMT const& the_im, std::enable_if_t<std::is_base_of<PointSetBase,IMT>::value>* = nullptr )
{
    return the_im;
}


/**
 * @brief build a quadrature formula from a @p mesh and a polynomial order @p O to integrate exactly
 * @ingroup Quadrature
 * 
 * @code
 * 
 * @endcode
 * 
 * @tparam MeshT type of the mesh
 * @tparam T numerical type of the quadrature
 * @param mesh mesh to integrate on
 * @param O polynomial order integrated exactly
 * @return im_t<typename MeshT::element_type,T> the quadrature formula 
 */
template<typename MeshT,typename T = double>
im_t<typename MeshT::element_type,T> im( std::shared_ptr<MeshT> mesh, uint16_type O )
{
    return im_t<typename MeshT::element_type,T>{ O };
}

template<int DIM,
         int IMORDER,
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Entity = Simplex>
class IM
        :
        public IMGeneral<DIM, T, Entity>
{
    using super = IMGeneral<DIM, T, Entity>;
public:
    template<int DIM1,
             typename T1,
             template<uint16_type, uint16_type, uint16_type> class Entity1>
    struct apply
    {
        typedef IMGeneral<DIM1, T1, Entity1> type;
        
        
    };
    IM() : super( IMORDER ) {}
    explicit IM( uint16_type o ) : super( o ) {}
    IM( IM const& ) = default;
    IM( IM && ) = default;
    IM& operator=( IM const& ) = default;
    IM& operator=( IM && ) = default;
};


template<uint16_type IMORDER = invalid_uint16_type_value,
         template<class Convex, uint16_type O, typename T2> class QPS = Gauss>
struct _Q : public _QBase
{
    //static const int order = IMORDER;
    static const uint16_type CompileTimeOrder = IMORDER;

    template<int DIM,
             typename T,
             template<uint16_type, uint16_type, uint16_type> class Entity>
    struct Apply
    {
        typedef IMGeneral<DIM, T, Entity> type;
    };
    template<typename T,
             typename GeoEntityType>
    struct ApplyGeoEntity
    {
        static const uint16_type DIM = GeoEntityType::nRealDim;
        typedef typename mpl::if_<mpl::bool_<GeoEntityType::is_simplex>,
                                  mpl::identity<IMGeneral<DIM, T, Simplex> >,
                                  mpl::identity<IMGeneral<DIM, T, Hypercube> > >::type::type type;
    };

    template<int DIM,
             typename T,
             template<uint16_type, uint16_type, uint16_type> class Entity>
    typename Apply<DIM,T,Entity>::type apply( uint16_type O ) const
        {
            return typename Apply<DIM,T,Entity>::type( O );
        }

    template<int DIM,
             typename T,
             template<uint16_type, uint16_type, uint16_type> class Entity>
    typename Apply<DIM,T,Entity>::type apply() const
        {
            return typename Apply<DIM,T,Entity>::type( this->order() );
        }


    template<int DIM,
             typename T,
             template<uint16_type, uint16_type, uint16_type> class Entity>
    typename Apply<DIM,T,Entity>::type get() const
        {
            return typename Apply<DIM,T,Entity>::type( this->order() );
        }

    template< typename T,
              int DIM,
             template<uint16_type, uint16_type, uint16_type> class Entity>
    typename Apply<DIM,T,Entity>::type get( Entity<DIM,1,DIM> const& e ) const
        {
            return typename Apply<DIM,T,Entity>::type( this->order() );
        }
    template< typename T,
              int DIM,
              template<uint16_type, uint16_type, uint16_type> class Entity>
    typename Apply<DIM,T,Entity>::type get( Entity<DIM,1,DIM> && e ) const
        {
            return typename Apply<DIM,T,Entity>::type( this->order() );
        }

    template<typename T,
             typename GeoEntityType>
    typename ApplyGeoEntity<T,GeoEntityType>::type getGeoEntity() const
        {
            return typename ApplyGeoEntity<T,GeoEntityType>::type( this->order() );
        }

    template<int DIM,
             typename T,
             template<uint16_type, uint16_type, uint16_type> class Entity>
    struct ApplyIMGeneral
    {
        //typedef IMGeneral<DIM, IMORDER, T, Entity,QPS> type;
        typedef IMGeneral<DIM, T, Entity> type;
    };

    template<int DIM,
             typename T,
             template<uint16_type, uint16_type, uint16_type> class Entity>
    typename ApplyIMGeneral<DIM,T,Entity>::type applyIMGeneral( uint16_type O ) const 
        {
            return typename ApplyIMGeneral<DIM,T,Entity>::type( O );
        }

    template<int DIM,
             typename T,
             template<uint16_type, uint16_type, uint16_type> class Entity>
    typename ApplyIMGeneral<DIM,T,Entity>::type applyIMGeneral() const 
        {
            return typename ApplyIMGeneral<DIM,T,Entity>::type( this->order() );
        }
    
    template<typename ContextType>
    struct ApplyContext
    {
        static const int DIM = ContextType::PDim;
        typedef typename ContextType::value_type T;
#if 0
        typedef typename mpl::if_<mpl::and_<mpl::less_equal<mpl::int_<IMORDER>,mpl::int_<20> >,
                                            mpl::bool_<ContextType::element_type::is_simplex> >,
                                  mpl::identity<IMSimplex<DIM, T> >,
                                  typename mpl::if_<mpl::bool_<ContextType::element_type::is_simplex>,
                                                    mpl::identity<IMGeneral<DIM, IMORDER, T, Simplex,QPS> >,
                                                    mpl::identity<IMGeneral<DIM, IMORDER, T, Hypercube,QPS> > >::type>::type::type type;
#else
        typedef typename mpl::if_<mpl::bool_<ContextType::element_type::is_simplex>,
                                  mpl::identity<IMGeneral<DIM, T, Simplex> >,
                                  mpl::identity<IMGeneral<DIM, T, Hypercube> > >::type::type type;
#endif
    };
    template<typename ContextType>
    typename ApplyContext<ContextType>::type applyContext( uint16_type O ) const
        {
            return typename ApplyContext<ContextType>::type( O );
        }
    template<typename ContextType>
    typename ApplyContext<ContextType>::type applyContext() const
        {
            return typename ApplyContext<ContextType>::type( this->order() );
        }
    
    _Q()
        :
        M_order( (CompileTimeOrder!=invalid_uint16_type_value)?CompileTimeOrder:1 )
        {}
    explicit _Q( uint16_type M_order )
        :
        M_order( M_order )
        {}

    uint16_type order() const { return M_order; }
private:
    uint16_type M_order;
};


template<int IMORDER,
         int DIM,
         template<uint16_type, uint16_type, uint16_type> class Entity,
         template<class Convex, uint16_type O, typename T2> class QPS,
         typename T>
struct IMGeneric
{
    typedef typename _Q<IMORDER,QPS>::template Apply<DIM,T,Entity>::type type;

    type apply( uint16_type O ) const
        {
            return type( O );
        }

    type apply() const
        {
            return type( IMORDER );
        }
};

template <int IMORDER,
          int DIM,
          template <uint16_type, uint16_type, uint16_type> class Entity,
          template <class Convex, uint16_type O, typename T2> class QPS,
          typename T>
using imgeneric_t = typename IMGeneric<IMORDER,DIM,Entity,QPS,T>::type;

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
        mpl::identity<IMSimplex<Dim, T> >,
        mpl::identity<IM<Dim, IMORDER, T, Simplex> > >::type::type
{};
#endif

;
} // Feel


#endif /* __im_H */
