/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-04-26

  Copyright (C) 2013-2016 Feel++ Consortium

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
#ifndef FEELPP_DISCR_FSFUNCTIONALLINEARFREE_H
#define FEELPP_DISCR_FSFUNCTIONALLINEARFREE_H 1

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/fsfunctionallinear.hpp>
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

    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef typename space_type::value_type value_type;

    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef ExprType expr_type;

    FsFunctionalLinearFree( space_ptrtype space , expr_type expr ) :
        super_type( space ),
#if FEELPP_HAS_PETSC
        M_backend( backend_type::build( BACKEND_PETSC ) ),
#else
        M_backend( backend_type::build( BACKEND_EIGEN ) ),
#endif
        M_expr( expr )
    {}

    FsFunctionalLinearFree( space_ptrtype space, backend_ptrtype backend , expr_type expr ) :
        super_type( space ),
        M_backend( backend ),
        M_expr( expr )
    {}

    ~FsFunctionalLinearFree() override {}

    //return the expression
    expr_type expr()
    {
        return M_expr;
    }

    // apply the functional
    value_type
    operator()( const element_type& x ) const override
    {
        auto vector = M_backend->newVector( this->space() );
        form1( _test=this->space(),_vector=vector) = M_expr;
        vector->close();

        return M_backend->dot( *vector, x.container() );
    }

    //fill a vector to have the container
    void containerPtr( vector_ptrtype & vector_to_fill ) override
    {
        auto vector = M_backend->newVector( this->space() );
        form1( _test=this->space(),_vector=vector) = M_expr;
        vector->close();
        vector_to_fill = vector;
    }

    //fill a vector to have the container
    void container( vector_type & vector_to_fill ) override
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


template <typename ... Ts>
auto functionalLinearFree( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && space = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_space);
    auto && expr = args.get(_expr);
    auto && backend = args.get_else_invocable( _backend, [](){ return Feel::backend(); } );
    using space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(space)>>>;
    using expr_type = std::decay_t<decltype(expr)>;
    return std::make_shared<FsFunctionalLinearFree<space_type,expr_type>>( space , backend , expr );
}


}//Feel

#endif /* _FSFUNCTIONALLINEARFREE_HPP_ */
