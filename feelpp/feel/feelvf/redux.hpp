/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2019-12-07

  Copyright (C) 2019 Feel++ Consortium

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
#ifndef __FEELPP_VF_REDUX_H
#define __FEELPP_VF_REDUX_H 1

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class Sum
 * \brief sum of array elements
 *
 * @author Christophe Prud'homme
 * @see
 */
template <typename ExprT, typename OpRedux>
class Redux : public ExprDynamicBase
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
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef value_type evaluate_type;
    using op_redux_t = OpRedux;
    typedef Redux<ExprT,OpRedux> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Redux( expression_type const& __expr, op_redux_t op )
        : super( Feel::vf::dynamicContext( __expr ) ),
          M_expr( __expr ),
          M_op( op )
    {
    }
    Redux( Redux const& te )
        : super( te ), M_expr( te.M_expr ), M_op( te.M_op )
    {
    }
    ~Redux() = default;

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
    uint16_type polynomialOrder() const { return 2 * M_expr.polynomialOrder(); }

    //! expression is polynomial?
    constexpr bool isPolynomial() const { return false; }

    expression_type const& expression() const
    {
        return M_expr;
    }
    op_redux_t const& op() const
    {
        return M_op;
    }
    //@}

    //template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,std::shared_ptrvf::detail::gmc<0> > > >, typename Basis_j_t = Basis_i_t>
    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape expr_shape;
        using shape = shape_rankdown_t<expr_shape>;

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
            :
            M_tensor_expr( expr.expression(), geom, fev, feu ),
            M_op( expr.op() )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            : M_tensor_expr( expr.expression(), geom, fev ),
              M_op( expr.op() )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            : M_tensor_expr( expr.expression(), geom ),
              M_op( expr.op() )
        {
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update( geom );
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
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_op( c2, q, expr_shape(), [this]( int i, int j, int q ){ return this->M_tensor_expr.evalq( i, j, q ); } );
        }
      private:
        tensor_expr_type M_tensor_expr;
        op_redux_t M_op;
    };

  private:
    mutable expression_type M_expr;
    op_redux_t M_op;
};
/// \endcond

namespace detail
{
template<typename shape_t, typename Fun, typename T>
T OpReduxSum( T init, int c2, int q, Fun f )
{
    if constexpr ( shape_t::M > 1 )
    {
        T  res = init;
        for( int i = 0; i < shape_t::M; ++i )
            res += f( i, c2, q );
        return res;
    }
    else if constexpr ( shape_t::N > 1 )
    {
        T  res = init;
        for( int i = 0; i < shape_t::N; ++i )
            res += f( 0, i, q );
        return res;
    }
    else
        return init + f( 0, 0, q );

}
template<typename shape_t, typename Fun, typename T>
T OpReduxMean( T init, int c2, int q, Fun f )
{
    if constexpr ( shape_t::M > 1 )
    {
        T  res = 0.;
        for( int i = 0; i < shape_t::M; ++i )
            res += f( i, c2, q );
        return init + res/shape_t::M;
    }
    else if constexpr ( shape_t::N > 1 )
    {
        T  res = 0.;
        for( int i = 0; i < shape_t::N; ++i )
            res += f( 0, i, q );
        return init + res/shape_t::N;
    }
    else
        return init + f( 0, 0, q );

}

template<typename shape_t, typename Fun, typename T>
T OpReduxProd( T init, int c2, int q,  Fun f )
{
    if constexpr ( shape_t::M > 1 )
    {
        T  res = init;
        for( int i = 0; i < shape_t::M; ++i )
            res *= f( i, c2, q );
        return res;
    }
    else if constexpr ( shape_t::N > 1 )
    {
        value_type  res = init;
        for( int i = 0; i < shape_t::N; ++i )
            res *= f( 0, i, q );
        return res;
    }
    else
        return init * f( 0, 0, q );

}

};
/**
 * \brief compute the sum of element array expression \p ExprT
 * \return the sum of the elements of expression v along the first array dimension whose size does not equal 1.
 */
template <typename ExprT, typename T = typename  ExprT::value_type>
inline auto
sum( ExprT v, T init = 0., std::enable_if_t<std::is_base_of_v<ExprBase,ExprT>>* = nullptr )
{
    Redux redux_sum( v, [&]( int c2, int q, auto shape, auto const& f ) { return detail::OpReduxSum<decltype(shape)>( init, c2, q, f ); } );
    return Expr{ redux_sum };
}

/**
 * \brief compute the mean of element array expression \p ExprT
 * \return the mean of the elements of expression v along the first array dimension whose size does not equal 1.
 */
template <typename ExprT, typename T = typename  ExprT::value_type>
inline auto
mean( ExprT v, T init = 0., std::enable_if_t<std::is_base_of_v<ExprBase,ExprT>>* = nullptr )
{
    Redux redux_mean( v, [&]( int c2, int q, auto shape, auto const& f ) { return detail::OpReduxMean<decltype(shape)>( init, c2, q, f ); } );
    return Expr{ redux_mean };
}

/**
 * \brief compute the productpf  element array expression \p ExprT
 * \return the product of the elements of expression v along the first array dimension whose size does not equal 1.
 */
template <typename ExprT, typename T = typename  ExprT::value_type>
inline auto
prod( ExprT v, T init = 1., std::enable_if_t<std::is_base_of_v<ExprBase,ExprT>>* = nullptr )
{
    Redux redux_prod( v, [&]( int c2, int q, auto shape, auto const& f ) { return detail::OpReduxProd<decltype(shape)>( init, c2, q, f ); } );
    return Expr{ redux_prod };
}

} // namespace vf
} // namespace Feel
#endif /* __Inv_H */
