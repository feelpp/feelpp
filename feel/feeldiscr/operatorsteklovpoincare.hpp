/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake@imag.fr>
             Date: 2012-02-09

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
/**
   \file operatorsteklovpoincare.cpp
   \author Abdoulaye Samake <abdoulaye.samake@imag.fr>
   \date 2012-02-09
*/
#ifndef _OPERATORSTEKLOVPC_HPP_
#define _OPERATORSTEKLOVPC_HPP_

#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feeldiscr/operatorlift.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>
//#include<iostream>


namespace Feel
{

/**
 * \class OperatorSteklovPc
 * \brief OperatorSteklovPc made easy
 *
 * @author Abdoulaye Samake
 * @see OperatorLinear
 */
template<class fs_type>
class OperatorSteklovPc : public OperatorLinear<fs_type, fs_type>
{
    typedef OperatorLinear<fs_type, fs_type> super;

public :

    /** @name Typedefs
     */
    //@{
    typedef OperatorSteklovPc<fs_type> this_type;
    typedef OperatorLinear<fs_type, fs_type> super_type;
    typedef fs_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename super::backend_type backend_type;
    typedef typename super::backend_ptrtype backend_ptrtype;
    typedef FsFunctionalLinear<fs_type> image_element_type;
    typedef typename image_element_type::value_type value_type;
    //@}
    /** @name Constructors, destructor
     */
    //@{

    OperatorSteklovPc( space_ptrtype Xh, backend_ptrtype backend = Backend<double>::build( BACKEND_PETSC ) )
        :
        super_type( Xh, Xh, backend ),
        M_backend( backend ),
        M_Xh( Xh )
    {}

    ~OperatorSteklovPc() {}
    //@}

    /** @name  Methods
     */
    //@{
    template<typename Args,typename IntEltsDefault>
    struct integrate_type
    {
        typedef typename vf::detail::clean_type<Args,tag::expr>::type _expr_type;
        typedef typename vf::detail::clean2_type<Args,tag::range,IntEltsDefault>::type _range_type;
        typedef typename vf::detail::clean2_type<Args,tag::quad, _Q< vf::ExpressionOrder<_range_type,_expr_type>::value > >::type _quad_type;
        typedef typename vf::detail::clean2_type<Args,tag::quad1, _Q< vf::ExpressionOrder<_range_type,_expr_type>::value_1 > >::type _quad1_type;
    };

    BOOST_PARAMETER_MEMBER_FUNCTION(
        ( value_type ),
        steklovpc,
        tag,
        ( required
          ( domain,  * )
          ( image, * )
        )
        ( optional
          ( quad,   *, ( typename integrate_type<Args,decltype( elements( this->M_Xh->mesh() ) )>::_quad_type() ) )
          ( quad1,   *, ( typename integrate_type<Args,decltype( elements( this->M_Xh->mesh() ) )>::_quad1_type() ) )
          ( geomap, *, GeomapStrategyType::GEOMAP_OPT )
        ) )

    {
        using namespace vf;

        auto op_lift = operatorLift( this->M_Xh,this->M_backend );
        auto domain_lift = op_lift->lift( _range=this->M_Xh->mesh(),_expr=idv( domain ) );
        auto image_lift = op_lift->lift( _range=this->M_Xh->mesh(),_expr=idv( image ) );
        value_type steklovpcr = integrate( _range=elements( this->M_Xh->mesh() ), _expr=gradv( domain_lift )*trans( gradv( image_lift ) ), _quad=quad, _quad1=quad1 );
        return steklovpcr;
    }

    template<typename First, typename Second>
    value_type
    operator()( First const& first ,Second const& second )
    {
        return this->steklovpc( first, second );
    }

    //@}


private :

    space_ptrtype M_Xh;
    backend_ptrtype M_backend;

};//OperatorSteklovPc

/**
 * this function returns a \c OperatorSteklovPc \c shared_ptr with
 *
 * \param pace
 * \param backend
 */

template<typename space_type>
boost::shared_ptr< OperatorSteklovPc<space_type> >
operatorSteklPc( boost::shared_ptr<space_type> const& space,
                 typename OperatorSteklovPc<space_type>::backend_ptrtype const& backend = Backend<double>::build( BACKEND_PETSC ) )
{
    typedef OperatorSteklovPc<space_type> StekPc_type;
    boost::shared_ptr<StekPc_type> Steklov( new StekPc_type( space, backend ) );
    return Steklov;
}

} //namespace Feel


#endif
