/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-02-10

  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

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
   \file trace.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-02-10
 */
#ifndef __Trace_H
#define __Trace_H 1

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class Trace
 * \brief trace of a matrix
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename ExprT>
class Trace
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
    typedef value_type evaluate_type;
    typedef Trace<ExprT> this_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Trace( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    Trace( Trace const & te )
        :
        M_expr( te.M_expr )
    {}
    ~Trace()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    expression_type const& expression() const
    {
        return M_expr;
    }

    //@}

    //template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,boost::shared_ptrvf::detail::gmc<0> > > >, typename Basis_j_t = Basis_i_t>
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape expr_shape;
        BOOST_MPL_ASSERT_MSG( ( boost::is_same<mpl::int_<expr_shape::M>,mpl::int_<expr_shape::N> >::value ), INVALID_TENSOR_SHOULD_BE_RANK_2_OR_0, ( mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N> ) );
        typedef Shape<expr_shape::nDim,Scalar,false,false> shape;


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
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev )
        {
        }

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


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            value_type res = value_type( 0 );

            for ( uint16_type l = 0; l < expr_shape::M; ++l )
                res += M_tensor_expr.evalijq( i, j, l, l, q );

            return res;
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            value_type res = value_type( 0 );

            for ( uint16_type l = 0; l < expr_shape::M; ++l )
                res += M_tensor_expr.evalijq( i, j, l, l, q, mpl::int_<PatternContext>() );

            return res;
        }

        value_type
        evaliq( uint16_type i, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            value_type res = value_type( 0 );

            for ( uint16_type l = 0; l < expr_shape::M; ++l )
                res += M_tensor_expr.evaliq( i, l, l, q );

            return res;
        }

        value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            value_type res = value_type( 0 );

            for ( uint16_type l = 0; l < expr_shape::M; ++l )
                res += M_tensor_expr.evalq( l, l, q );

            return res;
        }

        tensor_expr_type M_tensor_expr;
    };

private:
    mutable expression_type  M_expr;
};
/// \endcond

/**
 * \brief trace of the expression tensor
 */
template<typename ExprT>
inline
Expr< Trace<ExprT> >
trace( ExprT v )
{
    typedef Trace<ExprT> trace_t;
    return Expr< trace_t >(  trace_t( v ) );
}

}
}
#endif /* __Trace_H */
