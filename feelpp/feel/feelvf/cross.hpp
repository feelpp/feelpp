/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-09-01

  Copyright (C) 2013-2014 Université Joseph Fourier (Grenoble I)

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
#ifndef FEELPP_FEELVF_CROSS_HPP
#define FEELPP_FEELVF_CROSS_HPP 1

#include <Eigen/Core>

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class CrossProduct
 * \brief cross products
 *
 * @author Christophe Prud'homme
 * @see
 */
template <typename ExprL, typename ExprR>
class CrossProduct : public ExprDynamicBase
{
  public:
    using super = ExprDynamicBase;
    static const size_type context = ExprL::context | ExprR::context;
    static const bool is_terminal = false;

    template <typename Func>
    struct HasTestFunction
    {
        static const bool result =
            ExprL::template HasTestFunction<Func>::result ||
            ExprR::template HasTestFunction<Func>::result;
    };

    template <typename Func>
    struct HasTrialFunction
    {
        static const bool result =
            ExprL::template HasTrialFunction<Func>::result ||
            ExprR::template HasTrialFunction<Func>::result;
    };

    template <typename Func>
    static const bool has_test_basis = ExprL::template HasTestFunction<Func>::result ||
                                       ExprR::template HasTestFunction<Func>::result;
    template <typename Func>
    static const bool has_trial_basis = ExprL::template HasTrialFunction<Func>::result ||
                                        ExprR::template HasTrialFunction<Func>::result;
    using test_basis = typename ExprL::test_basis;
    using trial_basis = typename ExprL::trial_basis;

    template <typename... TheExpr>
    struct Lambda
    {
        typedef CrossProduct<typename ExprL::template Lambda<TheExpr...>::type, typename ExprR::template Lambda<TheExpr...>::type> type;
    };

    template <typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e ) { return typename Lambda<TheExpr...>::type( M_left_expr( e... ), M_right_expr( e... ) ); }

    /** @name Typedefs
     */
    //@{

    typedef ExprL left_expression_type;
    typedef ExprR right_expression_type;
    typedef typename left_expression_type::value_type value_type;
    using evaluate_type = Eigen::Matrix<value_type,
                                        (left_expression_type::evaluate_type::RowsAtCompileTime == 3)? 3 : ((left_expression_type::evaluate_type::RowsAtCompileTime == 2)? 1 : Eigen::Dynamic),
                                        1 >;
    typedef CrossProduct<ExprL, ExprR> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit CrossProduct( left_expression_type const& left_expr,
                           right_expression_type const& right_expr )
        : super( Feel::vf::dynamicContext( left_expr ) | Feel::vf::dynamicContext( right_expr ) ),
          M_left_expr( left_expr ),
          M_right_expr( right_expr )
    {
    }
    CrossProduct( CrossProduct const& te )
        : super( te ),
          M_left_expr( te.M_left_expr ),
          M_right_expr( te.M_right_expr )
    {
    }
    ~CrossProduct() = default;

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
    uint16_type polynomialOrder() const { return M_left_expr.polynomialOrder() + M_right_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_left_expr.isPolynomial() && M_right_expr.isPolynomial(); }

    left_expression_type const& left() const
    {
        return M_left_expr;
    }
    right_expression_type const& right() const
    {
        return M_right_expr;
    }

    //! evaluate the expression without context
    evaluate_type evaluate(bool p ) const
        {
            auto leval = M_left_expr.evaluate(p);
            auto reval = M_right_expr.evaluate(p);
            CHECK( leval.rows() == reval.rows() && leval.cols() == reval.cols() ) << "vector should be same dim";
            CHECK( leval.rows() == 2 || leval.rows() == 3 ) << "only vector of dim 2 or 3";
            CHECK( leval.cols() == 1 ) << "not a vector";
            if ( leval.rows() == 2 )
            {
                evaluate_type res(1,1);
                res(0,0) = leval(0,0)*reval(1,0) - leval(1,0)*reval(0,0);
                return res;
            }
            else
            {
                evaluate_type res(3,1);
                res(0,0) = leval(1,0)*reval(2,0) - leval(2,0)*reval(1,0);
                res(1,0) = leval(2,0)*reval(0,0) - leval(0,0)*reval(2,0);
                res(2,0) = leval(0,0)*reval(1,0) - leval(1,0)*reval(0,0);
                return res;
            }
        }

    //@}

    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        typedef typename left_expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_tensor_expr_type;
        typedef typename right_expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_tensor_expr_type;
        typedef typename l_tensor_expr_type::value_type value_type;

        typedef typename l_tensor_expr_type::shape left_shape;
        typedef typename r_tensor_expr_type::shape right_shape;
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<left_shape::nDim>, mpl::int_<2>>,
                                  mpl::identity<Shape<2, Scalar, false, false>>,
                                  mpl::identity<Shape<3, Vectorial, false, false>>>::type::type shape;
        static const bool l_is_terminal = left_expression_type::is_terminal;
        static const bool r_is_terminal = right_expression_type::is_terminal;

        BOOST_MPL_ASSERT_MSG( ( left_shape::nDim > 1 ),
                              CROSS_INVALID_DIMENSION,
                              (mpl::int_<left_shape::nDim>, mpl::int_<right_shape::nDim>));
        BOOST_MPL_ASSERT_MSG( left_shape::nDim == right_shape::nDim,
                              CROSS_INVALID_DIMENSION_LEFT_AND_RIGHT_SHOULD_BE_THE_SAME,
                              (mpl::int_<left_shape::nDim>, mpl::int_<right_shape::nDim>));
        BOOST_MPL_ASSERT_MSG( ( left_shape::M == right_shape::M ) && ( left_shape::N == right_shape::N ),
                              CROSS_INVALID_RANK_LEFT_AND_RIGHT_SHOULD_BE_THE_SAME,
                              (mpl::int_<left_shape::M>, mpl::int_<right_shape::M>,
                               mpl::int_<left_shape::N>, mpl::int_<right_shape::N>));

        template <class Args>
        struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = l_tensor_expr_type::is_zero::value || r_tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            : M_l_tensor_expr( expr.left(), geom, fev, feu ),
              M_r_tensor_expr( expr.right(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            : M_l_tensor_expr( expr.left(), geom, fev ),
              M_r_tensor_expr( expr.right(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            : M_l_tensor_expr( expr.left(), geom ),
              M_r_tensor_expr( expr.right(), geom )
        {
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_l_tensor_expr( std::true_type{}, exprExpanded.left(), ttse, expr.left(), geom, theInitArgs... ),
            M_r_tensor_expr( std::true_type{}, exprExpanded.right(), ttse, expr.right(), geom, theInitArgs... )
            {}

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_l_tensor_expr.update( geom, fev, feu );
            M_r_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_l_tensor_expr.update( geom, fev );
            M_r_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_l_tensor_expr.update( geom );
            M_r_tensor_expr.update( geom );
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            M_l_tensor_expr.update( std::true_type{}, exprExpanded.left(), ttse, geom, theUpdateArgs... );
            M_r_tensor_expr.update( std::true_type{}, exprExpanded.right(), ttse, geom, theUpdateArgs... );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            return evalijq( i, j, cc1, cc2, q, mpl::int_<left_shape::nDim>() );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<2> ) const
        {
            double res = M_l_tensor_expr.evalijq( i, j, 0, 0, q ) * M_r_tensor_expr.evalijq( i, j, 1, 0, q );
            res -= M_l_tensor_expr.evalijq( i, j, 1, 0, q ) * M_r_tensor_expr.evalijq( i, j, 0, 0, q );
            return res;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<3> ) const
        {
            double res = M_l_tensor_expr.evalijq( i, j, ( cc1 + 1 ) % 3, 0, q ) * M_r_tensor_expr.evalijq( i, j, ( cc1 + 2 ) % 3, 0, q );
            res -= M_l_tensor_expr.evalijq( i, j, ( cc1 + 2 ) % 3, 0, q ) * M_r_tensor_expr.evalijq( i, j, ( cc1 + 1 ) % 3, 0, q );
            return res;
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            return evaliq( i, cc1, cc2, q, mpl::int_<left_shape::nDim>() );
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<2> ) const
        {
            double res = M_l_tensor_expr.evaliq( i, 0, 0, q ) * M_r_tensor_expr.evaliq( i, 1, 0, q );
            res -= M_l_tensor_expr.evaliq( i, 1, 0, q ) * M_r_tensor_expr.evaliq( i, 0, 0, q );
            return res;
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<3> ) const
        {
            double res = M_l_tensor_expr.evaliq( i, ( cc1 + 1 ) % 3, 0, q ) * M_r_tensor_expr.evaliq( i, ( cc1 + 2 ) % 3, 0, q );
            res -= M_l_tensor_expr.evaliq( i, ( cc1 + 2 ) % 3, 0, q ) * M_r_tensor_expr.evaliq( i, ( cc1 + 1 ) % 3, 0, q );
            return res;
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            return evalq( cc1, cc2, q, mpl::int_<left_shape::nDim>() );
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<2> ) const
        {
            double res = M_l_tensor_expr.evalq( 0, 0, q ) * M_r_tensor_expr.evalq( 1, 0, q );
            res -= M_l_tensor_expr.evalq( 1, 0, q ) * M_r_tensor_expr.evalq( 0, 0, q );
            return res;
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::int_<3> ) const
        {
            double res = M_l_tensor_expr.evalq( ( cc1 + 1 ) % 3, 0, q ) * M_r_tensor_expr.evalq( ( cc1 + 2 ) % 3, 0, q );
            res -= M_l_tensor_expr.evalq( ( cc1 + 2 ) % 3, 0, q ) * M_r_tensor_expr.evalq( ( cc1 + 1 ) % 3, 0, q );
            return res;
        }

      private:
        l_tensor_expr_type M_l_tensor_expr;
        r_tensor_expr_type M_r_tensor_expr;
    };

  private:
    mutable left_expression_type M_left_expr;
    mutable right_expression_type M_right_expr;
};
/// \endcond

//!
//! @return an expression which is the cross product of two expressions \p l and \p r`
//!
template <typename ExprL, typename ExprR>
inline Expr<CrossProduct<ExprL, ExprR>>
cross( ExprL l, ExprR r )
{
    typedef CrossProduct<ExprL, ExprR> product_t;
    return Expr<product_t>( product_t( l, r ) );
}

} // namespace vf

} // namespace Feel
#endif /* __Products_H */
