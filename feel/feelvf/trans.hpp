/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-13

  Copyright (C) 2014 Feel++ Consortium

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
   \file trans.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_TRANS_HPP
#define FEELPP_VF_TRANS_HPP 1

namespace Feel
{
namespace vf
{

/*!
  \ingroup vf
  \brief Transpose an scalar, vectorial or matricial expression


  @author Christophe Prud'homme
*/
template<typename ExprT>
class Trans
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ExprT::imorder;
    static const bool imIsPoly = ExprT::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef typename expression_type::evaluate_type evaluate_type;
    typedef Trans<ExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Trans( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    ~Trans()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{
    template<typename... TheExpr>
    struct Lambda
    {
        typedef Trans<typename expression_type::template Lambda<TheExpr...>::type> type;
    };
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  ) { return trans(M_expr(e...)); }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename Transpose<typename tensor_expr_type::shape>::type shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( expr.expression(), geom, fev, feu )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
        }

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_tensor_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensor_expr.update( geom, face );
        }
        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
            M_tensor_expr.updateContext( ctx );
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor_expr.evalijq( i, j, c2, c1, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return M_tensor_expr.evalijq( i, j, c2, c1, q, mpl::int_<PatternContext>() );
        }


        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor_expr.evaliq( i, c2, c1, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor_expr.evalq( c2, c1, q );
        }

        tensor_expr_type M_tensor_expr;
    };

    //@}

    /** @name Accessors
     */
    //@{

    bool isSymetric() const
    {
        return M_expr.isSymetric();
    }

    expression_type const& expression() const
    {
        return M_expr;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //@}

protected:

private:

    mutable expression_type  M_expr;
};

/**
   \ingroup vf

   Transpose an expression \p v. The expression can be scalar (no-op), vectorial
   or matricial.

   \code
   // define a column vector of two element
   auto e = vec(Px(),Py());
   // transpose it and get a row vector expression of two elements
   auto eT = trans(e);
   \endcode

   \return the transposed expression
 */
template<typename ExprT>
inline
Expr< Trans<ExprT> >
trans( ExprT v )
{
    typedef Trans<ExprT> trans_t;
    return Expr< trans_t >(  trans_t( v ) );
}


} // vf
} // Feel
#endif /* FEELPP_VF_TRANS_HPP */
