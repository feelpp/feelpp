/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2013-04-07

   Copyright (C) 2004 EPFL
   Copyright (C) 2006-2012 Université Joseph Fourier (Grenoble I)

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
   \file ReducedBasisSpace.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-07
*/
#ifndef __ReducedBasisSpace_H
#define __ReducedBasisSpace_H 1

#include <boost/shared_ptr.hpp>

#include <feel/feel.hpp>
#include <feel/feelcrb/crb.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>


namespace Feel
{

/**
 * \class ReducedBasisSpace
 * \brief Reduced Basis Space class
 *
 * @author Christophe Prud'homme
 * @see
 */

template<typename ModelType,
         typename MeshType,
         typename A1 = parameter::void_,
         typename A2 = parameter::void_,
         typename A3 = parameter::void_,
         typename A4 = parameter::void_>
class ReducedBasisSpace : public FunctionSpace<MeshType,A1,A2,A3,A4>
{

    typedef FunctionSpace<MeshType,A1,A2,A3,A4> super;

public :

    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;

    typedef typename super::element_type space_element_type;

    typedef decltype( super::context() ) super_ctx_type;

    ReducedBasisSpace( boost::shared_ptr<ModelType> model , boost::shared_ptr<MeshType> mesh )
        :
        super ( mesh )
    {
        LOG( INFO ) <<" ReducedBasisSpace constructor " ;
    }

    void setContext( super_ctx_type ctx)
    {
        LOG( INFO ) <<"add an existing context from EIM ( for example )";
    }

    class Basis
    {
        Basis() {}
        ~Basis() {}

        /*
         * add a new basis
         */
        void add( space_element_type const & e )
        {
        }

        /*
         * visualize basis of rb space
         */
        void visualize ()
        {
        }
    private :
        std::vector< space_element_type > M_basis;
    };//Basis

    /*
     *  contains coefficients in the RB expression
     *  i.e. contains coefficients u_i in
     *  u_N ( \mu ) = \sum_i^N u_i ( \mu ) \xi_i( x )
     *  where \xi_i( x ) are basis function
     */
    template<typename T = double,  typename Cont = Eigen::Matrix<T,Eigen::Dynamic,1> >
    class Element
        :
        public Cont,boost::addable<Element<T,Cont> >, boost::subtractable<Element<T,Cont> >
    {
    public:
        typedef T value_type;
    };

private :
    decltype( super::context() ) M_ctx ;

};//ReducedBasisSpace

template<int Order, typename ModelType , typename MeshType>
inline
boost::shared_ptr<ReducedBasisSpace<ModelType,MeshType,bases<Lagrange<Order,Scalar,Continuous>>>>
RbSpacePch(  boost::shared_ptr<ModelType> model , boost::shared_ptr<MeshType> mesh  )
{
    return ReducedBasisSpace<ModelType,MeshType,bases<Lagrange<Order,Scalar,Continuous>>>::New( model , mesh );
}


}//namespace Feel

#endif
