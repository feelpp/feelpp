/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-06-23

  Copyright (C) 2006 EPFL
  Copyright (C) 2011-2015 Feel++ Consortium

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
#ifndef FEELPP_IMSIMPLEX_HPP
#define FEELPP_IMSIMPLEX_HPP 1

#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelpoly/gauss.hpp>
#include <feel/feelpoly/imfactory.hpp>
#include <feel/feelpoly/geomap.hpp>


namespace Feel
{
namespace ublas = boost::numeric::ublas;

template<typename Convex,typename T>
typename mpl::if_<Feel::is_simplex<Convex>,
                  mpl::identity<std::shared_ptr<GT_Lagrange<Convex::nDim,1,Convex::nDim,Simplex,T>>>,
                  mpl::identity<std::shared_ptr<GT_Lagrange<Convex::nDim,1,Convex::nDim,Hypercube,T>>> >::type::type
    makeGeometricTransformation();

/**
 * @brief economical quadratures for simplices
 *
 * \ingroup Quadrature
 */
template<int Dim,typename T>
class IMSimplex
    :
        public PointSetQuadrature<Simplex<Dim,1> , T, index_type>
{
    typedef PointSetQuadrature<Simplex<Dim,1> , T, index_type> super;
public:

    /** @name Typedefs
     */
    //@{
    typedef Simplex<Dim,1> convex_type;
    typedef T value_type;
    typedef ublas::matrix<value_type,ublas::column_major> matrix_type;
    typedef ublas::vector<value_type> vector_type;

    using quad_type = IMBase<value_type>;

    static constexpr uint16_type Dim_m_1 = (Dim==0)?0:Dim-1;
    using face_quad_type = typename mpl::if_<mpl::bool_<(Dim_m_1>=1)>,mpl::identity<IMSimplex<Dim_m_1,T>>, mpl::identity<Gauss<Simplex<Dim_m_1,1>, invalid_uint16_type_value,T>>>::type::type;
    typedef IMSimplex<Dim,T> parent_quadrature_type;
    static const uint16_type nDim = Dim;
    static const uint16_type nRealDim = Dim;
    static const uint16_type nOrder = invalid_uint16_type_value;;

    //@}

    /** @name Constructors, destructor
     */
    //@{
    
    
    IMSimplex() = default;
    
    /**
     * build quadrature to intergrate polynomial up to degree \c order
     */
    explicit IMSimplex( uint16_type order ): super( order ) 
        {
            auto gm = makeGeometricTransformation<convex_type,T>();
            auto face_qr = std::make_shared<face_quad_type>(order);
            this->constructQROnFace( makeReferenceConvex<convex_type,nDim,1,nRealDim>(), gm, face_qr );
        }

    ~IMSimplex() = default;
    
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
    IMSimplex& operator=( IMSimplex const & i ) = default;

    /**
     * create quadrature rule that integrates exactly up to polynomial order \c p
     */
    void create( uint16_type order ) override
        {
            super::create(order);
            auto gm = makeGeometricTransformation<convex_type,T>();
            auto face_qr = std::make_shared<face_quad_type>(order);
         
            this->constructQROnFace( Reference<convex_type,nDim,1,nRealDim>(), gm, face_qr );
        }

    //@}


};

}

#endif /* FEELPP_IMSIMPLEX_HPP */
