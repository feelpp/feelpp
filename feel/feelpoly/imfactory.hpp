/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 04 mai 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_IMFACTORY_HPP
#define FEELPP_IMFACTORY_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

namespace Feel {

//!
//! base class for integration method
//!
template<typename T>
class IMBase
{
public:
    using value_type = T;

    IMBase() : M_dim(0), M_order(0), M_npoints(0) {};

    IMBase( uint16_type d, uint16_type o, uint16_type np )
        :
        M_dim(d),
        M_order(o),
        M_npoints(np)
        {}
    IMBase( IMBase const& ) = default;
    IMBase( IMBase && ) = default;
    IMBase& operator=( IMBase const& i ) = default;
    IMBase& operator=( IMBase && i ) = default;
    virtual ~IMBase() = default;
    
    uint16_type dimension() const noexcept { return M_dim; }
    uint16_type integrationDegree() const noexcept{ return M_order; }
    //!
    //! @return the degree of the quadrature
    //!
    uint16_type degree() const noexcept{ return M_order; }

    //!
    //! @return the number of quadrature points and weights
    //!
    uint16_type numberOfPoints() const noexcept { return M_npoints; }
    
        
    // structure that holds both points and weights
    std::vector<value_type> q;

    IMBase<T>* operator()() { return this; }

    //! set the quadrature to be defined
    void setDefined() { M_created = true; }
    //! @return true if the quadrature is defined, false otherwise
    bool isDefined() { return M_created; }
protected:
    bool M_created = false;
private:
    uint16_type M_dim;
    uint16_type M_order;
    uint16_type M_npoints;
    
};

template<typename T>
using IMFactory = Feel::Singleton< Feel::Factory< IMBase<T>, std::string > >;

template<typename T>
struct IMMaxOrder
{
    IMMaxOrder( uint16_type maxOrder ) : M_maxOrder( maxOrder ) {}
    IMMaxOrder( IMMaxOrder const& ) = default;
    IMMaxOrder( IMMaxOrder && ) = default;
    IMMaxOrder& operator=( IMMaxOrder const& i ) = default;
    IMMaxOrder& operator=( IMMaxOrder && i ) = default;
    uint16_type maxOrder() const { return M_maxOrder; }
private :
    uint16_type M_maxOrder;
};
template<typename T>
using IMMaxOrderFactory = Feel::Singleton< std::map< std::string, IMMaxOrder<T> > >;



#if 0
# define DIMS BOOST_PP_TUPLE_TO_LIST(2,(Triangle,Tetrahedra))
# define ORDERS BOOST_PP_TUPLE_TO_LIST(21,(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
# define SHAPES1 BOOST_PP_TUPLE_TO_LIST(1, (simplex))
//# define SHAPES1 BOOST_PP_TUPLE_TO_LIST(1, (hypercube))


#define FACTORY1NAME( LDIM, LORDER, LSHAPE )                            \
    BOOST_PP_STRINGIZE(im BOOST_PP_LPAREN() LDIM BOOST_PP_COMMA() LORDER BOOST_PP_COMMA() BOOST_PP_ARRAY_ELEM(0,LSHAPE) BOOST_PP_RPAREN())

# define FACTORY1(LDIM,LORDER,LSHAPE )                                  \
    const bool BOOST_PP_CAT( BOOST_PP_CAT( mesh, LDIM ), LORDER)  =     \
        IMFactory<double>::type::instance().registerProduct( boost::to_lower_copy(boost::algorithm::erase_all_copy( std::string( FACTORY1NAME(LDIM, LORDER, LSHAPE ) ), " " ) ), \
        *new Feel::detail::BOOST_PP_CAT(IM,BOOST_PP_ARRAY_ELEM(1,LSHAPE),Domain)(LDIM,LORDER) ); \
    const bool BOOST_PP_CAT( BOOST_PP_CAT( mesh, LDIM ), LORDER)  =        \
         IMFactory<float>::type::instance().registerProduct( boost::to_lower_copy(boost::algorithm::erase_all_copy( std::string( FACTORY1NAME(LDIM, LORDER, LSHAPE ) ), " " ) ), \
         *new Feel::detail::BOOST_PP_CAT(IM,BOOST_PP_ARRAY_ELEM(1,LSHAPE),Domain)(LDIM,LORDER) );

# define FACTORY1_OP(_, GDO) FACTORY1 GDO

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY1_OP, 3, ( DIMS, ORDERS, SHAPES1 ) )
#endif

}
#endif
