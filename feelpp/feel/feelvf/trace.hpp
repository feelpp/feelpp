/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-02-10

  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2016 Feel++ Consortium

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
template <typename ExprT>
class Trace : public ExprDynamicBase
{
  public:
    using super = ExprDynamicBase;
    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    template <typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template <typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    template <typename Func>
    static const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template <typename Func>
    static const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    using evaluate_type = Eigen::Matrix<value_type,1,1>;
    typedef Trace<ExprT> this_type;

    template <typename... TheExpr>
    struct Lambda
    {
        typedef Trace<typename expression_type::template Lambda<TheExpr...>::type> type;
    };
    template <typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e ) { return typename Lambda<TheExpr...>::type( M_expr( e... ) ); }

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Trace( expression_type const& __expr )
        : super( Feel::vf::dynamicContext( __expr ) ),
          M_expr( __expr )
    {
    }
    Trace( Trace const& te )
        : M_expr( te.M_expr )
    {
    }
    ~Trace()
    {
    }

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

    //! polynomial order
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    expression_type const& expression() const
    {
        return M_expr;
    }

    //! evaluate the expression without context
    evaluate_type evaluate(bool p,  worldcomm_ptr_t const& worldcomm ) const
        {
            return evaluate_type::Constant( M_expr.evaluate(p,worldcomm).trace() );
        }

    //@}

    //template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,std::shared_ptrvf::detail::gmc<0> > > >, typename Basis_j_t = Basis_i_t>
    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape expr_shape;
        BOOST_MPL_ASSERT_MSG( ( boost::is_same<mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N>>::value ), INVALID_TENSOR_SHOULD_BE_RANK_2_OR_0, (mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N>));
        typedef Shape<expr_shape::nDim, Scalar, false, false> shape;

        template <class Args>
        struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            : M_tensor_expr( expr.expression(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            : M_tensor_expr( expr.expression(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            : M_tensor_expr( expr.expression(), geom )
        {
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
        template <int PatternContext>
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
    mutable expression_type M_expr;
};
/// \endcond

/**
 * \brief trace of the expression tensor
 */
template <typename ExprT>
inline Expr<Trace<ExprT>>
trace( ExprT v, std::enable_if_t<std::is_base_of_v<ExprBase,ExprT>>* = nullptr )
{
    typedef Trace<ExprT> trace_t;
    return Expr<trace_t>( trace_t( v ) );
}

} // namespace vf
} // namespace Feel
#endif /* __Trace_H */
