/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-04-26

  Copyright (C) 2013 Feel++ Consortium

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
   \file fsfunctionallinearfree.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-26
 */
#ifndef __FSFUNCTIONALLINEARFREE_H
#define __FSFUNCTIONALLINEARFREE_H 1

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/fsfunctional.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

template<class Space, class ExprType>
class FsFunctionalLinearFree : public FsFunctionalLinear<Space>
{
public:

    typedef FsFunctionalLinearFree<Space,ExprType> this_type;
    typedef FsFunctionalLinear<Space> super_type;

    typedef Space space_type;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef typename space_type::value_type value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef ExprType expr_type;

    FsFunctionalLinearFree( space_ptrtype space , expr_type expr ) :
        super_type( space ),
        M_backend( backend_type::build( BACKEND_PETSC ) ),
        M_expr( expr )
    {}

    FsFunctionalLinearFree( space_ptrtype space, backend_ptrtype backend , expr_type expr ) :
        super_type( space ),
        M_backend( backend ),
        M_expr( expr )
    {}

    virtual ~FsFunctionalLinearFree() {}

    //return the expression
    expr_type expr()
    {
        return M_expr;
    }

    // apply the functional
    virtual value_type
    operator()( const element_type& x ) const
    {
        auto vector = M_backend->newVector( this->space() );
        form1( _test=this->space(),_vector=vector) = M_expr;
        vector->close();

        return M_backend->dot( *vector, x.container() );
    }

    //fill a vector to have the container
    virtual void containerPtr( vector_ptrtype & vector_to_fill )
    {
        auto vector = M_backend->newVector( this->space() );
        form1( _test=this->space(),_vector=vector) = M_expr;
        vector->close();
        vector_to_fill = vector;
    }

    //fill a vector to have the container
    virtual void container( vector_type & vector_to_fill )
    {
        auto vector = M_backend->newVector( this->space() );
        form1( _test=this->space(),_vector=vector) = M_expr;
        vector->close();

        vector_to_fill = *vector;
    }

    this_type& operator=( this_type const& m )
    {
        M_backend = m.M_backend;
        M_expr = m.M_expr;
        return *this;
    }


private:
    backend_ptrtype M_backend;
    expr_type M_expr;
};//FsFunctionalLinearFree


namespace detail
{

template<typename Args>
struct compute_functionalLinearFree_return
{
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::space>::type>::type::element_type space_type;
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::expr>::type>::type expr_type;

    typedef FsFunctionalLinearFree<space_type, expr_type> type;
    typedef boost::shared_ptr<FsFunctionalLinearFree<space_type,expr_type> > ptrtype;
};
}

BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::compute_functionalLinearFree_return<Args>::ptrtype ), // 1. return type
    functionalLinearFree,                        // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required
      ( space,    *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
      ( expr ,   * )
    ) // required
    ( optional
      ( backend,        *, backend() )
    ) // optionnal
)
{

    Feel::detail::ignore_unused_variable_warning( args );
    typedef typename Feel::detail::compute_functionalLinearFree_return<Args>::type functionalfree_type;
    typedef typename Feel::detail::compute_functionalLinearFree_return<Args>::ptrtype functionalfree_ptrtype;
    return functionalfree_ptrtype ( new functionalfree_type( space , backend , expr ) );

} // functionalLinearFree

}//Feel

#endif /* _FSFUNCTIONALLINEARFREE_HPP_ */
