/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-13

  Copyright (C) 2014-2016 Feel++ Consortium

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
#ifndef FEELPP_FEELVF_INNER_HPP
#define FEELPP_FEELVF_INNER_HPP 1

#include <Eigen/Core>
#include <feel/feelvf/expr.hpp>
#include <feel/feelcore/context.hpp>

namespace Feel
{
namespace vf
{
enum InnerProperties
{
    IS_SAME = ( 1 << 0 ),
    SQRT = ( 1 << 1 )
};

/// \cond detail
/**
 * \class Product
 * \brief inner and outer products
 *
 * @author Christophe Prud'homme
 * @see
 */
template <typename ExprL, typename ExprR, int Type = 1, int Props = NONE>
class Product : public ExprDynamicBase
{
  public:
    using super = ExprDynamicBase;
    static const size_type context = ExprL::context | ExprR::context;
    static const bool is_terminal = false;
    static const int product_type = Type;

    static const bool IsSame = has_value_v<Props,InnerProperties::IS_SAME>;//  Props & InnerProperties::IS_SAME;
    static const bool ApplySqrt =  has_value_v<Props,InnerProperties::SQRT>; //Props & InnerProperties::SQRT;

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
        typedef Product<typename ExprL::template Lambda<TheExpr...>::type, typename ExprR::template Lambda<TheExpr...>::type, Type, Props> type;
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
    using evaluate_type = Eigen::Matrix<value_type,1,1>;
    typedef Product<ExprL, ExprR, Type, Props> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Product( left_expression_type const& left_expr,
                      right_expression_type const& right_expr )
        :
        super( Feel::vf::dynamicContext( left_expr ) | Feel::vf::dynamicContext( right_expr ) ),
        M_left_expr( left_expr ),
        M_right_expr( right_expr )
    {
    }
    Product( Product const& te )
        :
        super( te ),
        M_left_expr( te.M_left_expr ),
        M_right_expr( te.M_right_expr )
    {
    }
    ~Product() = default;

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
            value_type res = 0;
            if constexpr ( Type == 1 )
            {
                auto leval = M_left_expr.evaluate(p);
                if constexpr( IsSame )
                    {
                        for ( uint16_type c2 = 0; c2 < leval.cols(); ++c2 )
                            for ( uint16_type c1 = 0; c1 < leval.rows(); ++c1 )
                            {
                                value_type val = leval( c1, c2 );
                                res += val * val;
                            }

                    }
                else
                {
                    auto reval = M_right_expr.evaluate(p);
                    for ( uint16_type c2 = 0; c2 < leval.cols(); ++c2 )
                        for ( uint16_type c1 = 0; c1 < leval.rows(); ++c1 )
                        {
                            res += leval( c1, c2 ) * reval( c1, c2 );
                        }
                }
            }

            if ( ApplySqrt )
                return evaluate_type::Constant( math::sqrt( res ) );
            else
                return evaluate_type::Constant( res );
        }

        void setParameterValues( std::map<std::string,value_type> const& mp )
        {
            M_left_expr.setParameterValues( mp );
            M_right_expr.setParameterValues( mp );
        }
        void updateParameterValues( std::map<std::string,double> & pv ) const
        {
            M_left_expr.updateParameterValues( pv );
            if constexpr( !IsSame )
                M_right_expr.updateParameterValues( pv );
        }

        template <typename SymbolsExprType>
        auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            auto newLeftExpr =  M_left_expr.applySymbolsExpr( se );
            using new_expr_left_type = std::decay_t<decltype(newLeftExpr)>;
            if constexpr( IsSame )
                {
                    using new_expr_type = Product<new_expr_left_type, new_expr_left_type, Type, Props>;
                    return new_expr_type( newLeftExpr, newLeftExpr );
                }
            else
            {
                auto newRightExpr =  M_right_expr.applySymbolsExpr( se );
                using new_expr_right_type = std::decay_t<decltype(newRightExpr)>;
                using new_expr_type = Product<new_expr_left_type, new_expr_right_type, Type, Props>;
                return new_expr_type( newLeftExpr, newRightExpr );
            }
        }

        template <typename TheSymbolExprType>
        bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            if constexpr( IsSame )
                return M_left_expr.hasSymbolDependency( symb, se );
            else
                return M_left_expr.hasSymbolDependency( symb, se ) || M_right_expr.hasSymbolDependency( symb, se );
        }

        template <typename TheSymbolExprType>
        void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            M_left_expr.dependentSymbols( symb,res,se );
            if constexpr( !IsSame )
                 M_right_expr.dependentSymbols( symb,res,se );
        }

        template <int diffOrder, typename TheSymbolExprType>
        auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
                   TheSymbolExprType const& se ) const
        {
            auto ldiffExpr = M_left_expr.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
            using expr_diff_left_type = std::decay_t<decltype(ldiffExpr)>;
            if constexpr( IsSame )
                {
                    using new_expr_nosqrt_type = Product<expr_diff_left_type,left_expression_type, Type, 0>;
                    auto resExprNoSqrt = 2.0*expr(new_expr_nosqrt_type( ldiffExpr,M_left_expr ));
                    if constexpr ( ApplySqrt )
                        return resExprNoSqrt/(2.0*expr(*this));
                    else
                        return resExprNoSqrt;
                }
            else
            {
                auto rdiffExpr = M_right_expr.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
                using expr_diff_right_type = std::decay_t<decltype(rdiffExpr)>;

                using new_expr1_nosqrt_type = Product<expr_diff_left_type,right_expression_type, Type, 0>;
                auto resExpr1 = expr(new_expr1_nosqrt_type( ldiffExpr,M_right_expr ));
                using new_expr2_nosqrt_type = Product<left_expression_type,expr_diff_right_type, Type, 0>;
                auto resExpr2 = expr(new_expr2_nosqrt_type( M_left_expr,rdiffExpr ));

                if constexpr ( ApplySqrt )
                    return (resExpr1 + resExpr2)/(2.0*expr(*this));
                else
                    return resExpr1 + resExpr2;
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
        typedef Shape<left_shape::nDim, Scalar, false, false> shape;
        static const bool l_is_terminal = left_expression_type::is_terminal;
        static const bool r_is_terminal = right_expression_type::is_terminal;
        template <class Args>
        struct sig
        {
            typedef value_type type;
        };

        BOOST_MPL_ASSERT_MSG( ( left_shape::M == right_shape::M ) && ( left_shape::N == right_shape::N ),
                              INVALID_RANK_LEFT_AND_RIGHT_SHOULD_BE_THE_SAME,
                              (mpl::int_<left_shape::M>, mpl::int_<right_shape::M>,
                               mpl::int_<left_shape::N>, mpl::int_<right_shape::N>));

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
            if constexpr ( !IsSame )
                M_r_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_l_tensor_expr.update( geom, fev );
            if constexpr ( !IsSame )
                M_r_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_l_tensor_expr.update( geom );
            if constexpr ( !IsSame )
                M_r_tensor_expr.update( geom );
        }
        template <typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            M_l_tensor_expr.updateContext( ctx... );
            if constexpr ( !IsSame )
                M_r_tensor_expr.updateContext( ctx... );
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            M_l_tensor_expr.update( std::true_type{}, exprExpanded.left(), ttse, geom, theUpdateArgs... );
            if constexpr ( !IsSame )
                M_r_tensor_expr.update( std::true_type{}, exprExpanded.right(), ttse, geom, theUpdateArgs... );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            value_type res = evalijq( i, j, cc1, cc2, q, mpl::bool_<IsSame>() );
            if constexpr ( ApplySqrt )
                return math::sqrt( res );
            else
                return res;
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_ ) const
        {
            return evalijq( i, j, cc1, cc2, q, mpl::bool_<l_is_terminal>(), mpl::bool_<r_is_terminal>(), mpl::bool_<false>() );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_, mpl::false_, mpl::false_ ) const
        {
            double res = 0;
            if ( Type == 1 )
            {
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += M_l_tensor_expr.evalijq( i, j, c1, c2, q ) * M_r_tensor_expr.evalijq( i, j, c1, c2, q );
                    }
            }
            return res;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_, mpl::true_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                //return ( M_l_tensor_expr.evalijq( i,j,q ).adjoint()*M_r_tensor_expr.evalijq( i,j,q ) ).trace();
#if 0
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++ c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++ c1 )
                    {
                        res += M_l_tensor_expr.evalijq( i,j,c1,c2,q )*M_r_tensor_expr.evalijq( i,j,c1,c2,q );
                    }
#else
                auto ltensor = M_l_tensor_expr.evalijq( i, j, q );
                auto rtensor = M_r_tensor_expr.evalijq( i, j, q );
#if 1
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += ltensor( c1, c2 ) * rtensor( c1, c2 );
                    }
#else
                //Eigen::array<Eigen::IndexPair<int>, 2> double_contraction_product_dims = { Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(1, 1) };
                //Eigen::Tensor<value_type, 0> r = ltensor.contract(rtensor, double_contraction_product_dims);
                //res = r(0);
                //res = (ltensor.transpose()*rtensor).trace();
                res += ltensor( 0, 0 ) * rtensor( 0, 0 );
                res += ltensor( 0, 1 ) * rtensor( 0, 1 );
                res += ltensor( 0, 2 ) * rtensor( 0, 2 );
#endif
#endif
            }
            return res;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_, mpl::false_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                auto const& ltensor = M_l_tensor_expr.evalijq( i, j, q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += ltensor( c1, c2 ) * M_r_tensor_expr.evalijq( i, j, c1, c2, q );
                    }
            }
            return res;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_, mpl::true_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                auto const& rtensor = M_r_tensor_expr.evalijq( i, j, q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += M_l_tensor_expr.evalijq( i, j, c1, c2, q ) * rtensor( c1, c2 );
                    }
            }
            return res;
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_ ) const
        {
            return evalijq( i, j, cc1, cc2, q, mpl::bool_<l_is_terminal>(), mpl::bool_<true>() );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_, mpl::true_ ) const
        {
            double res = 0;

            if ( Type == 1 )
            {
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        auto a = M_l_tensor_expr.evalijq( i, j, c1, c2, q );
                        res += a * a;
                    }
            }
            return res;
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_, mpl::true_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                //return ( M_l_tensor_expr.evalijq( i,j,q ).adjoint()*M_r_tensor_expr.evalijq( i,j,q ) ).trace();
#if 0
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++ c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++ c1 )
                    {
                        res += M_l_tensor_expr.evalijq( i,j,c1,c2,q )*M_r_tensor_expr.evalijq( i,j,c1,c2,q );
                    }
#else
                auto const& ltensor = M_l_tensor_expr.evalijq( i, j, q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        value_type val = ltensor( c1, c2 );
                        res += val * val;
                    }
#endif
            }
            return res;
        }

        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            value_type res = evaliq( i, cc1, cc2, q, mpl::bool_<IsSame>() );
            if ( ApplySqrt )
                return math::sqrt( res );
            else
                return res;
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_ ) const
        {
            return evaliq( i, cc1, cc2, q, mpl::bool_<l_is_terminal>(), mpl::bool_<r_is_terminal>(), mpl::bool_<false>() );
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_, mpl::false_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += M_l_tensor_expr.evaliq( i, c1, c2, q ) * M_r_tensor_expr.evaliq( i, c1, c2, q );
                    }
            }
            return res;
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_, mpl::true_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                auto const& ltensor = M_l_tensor_expr.evaliq( i, q );
                auto const& rtensor = M_r_tensor_expr.evaliq( i, q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += ltensor( c1, c2 ) * rtensor( c1, c2 );
                    }
            }
            return res;
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_, mpl::false_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                auto const& ltensor = M_l_tensor_expr.evaliq( i, q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += ltensor( c1, c2 ) * M_r_tensor_expr.evaliq( i, c1, c2, q );
                    }
            }
            return res;
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_, mpl::true_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                auto const& rtensor = M_r_tensor_expr.evaliq( i, q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += M_l_tensor_expr.evaliq( i, c1, c2, q ) * rtensor( c1, c2 );
                    }
            }
            return res;
        }

        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_ ) const
        {
            return evaliq( i, cc1, cc2, q, mpl::bool_<l_is_terminal>(), mpl::bool_<true>() );
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_, mpl::true_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        value_type val = M_l_tensor_expr.evaliq( i, c1, c2, q );
                        res += val * val;
                    }
            }
            return res;
        }
        value_type
        evaliq( uint16_type i, uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_, mpl::true_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                auto const& ltensor = M_l_tensor_expr.evaliq( i, q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        value_type val = ltensor( c1, c2 );
                        res += val * val;
                    }
            }
            return res;
        }

        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q ) const
        {
            value_type res = evalq( cc1, cc2, q, mpl::bool_<IsSame>() );
            if ( ApplySqrt )
                return math::sqrt( res );
            else
                return res;
        }

        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::bool_<false> ) const
        {
            return evalq( cc1, cc2, q, mpl::bool_<l_is_terminal>(), mpl::bool_<r_is_terminal>(), mpl::bool_<false>() );
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_, mpl::false_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += M_l_tensor_expr.evalq( c1, c2, q ) * M_r_tensor_expr.evalq( c1, c2, q );
                    }
            }
            return res;
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_, mpl::true_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                auto const& ltensor = M_l_tensor_expr.evalq( q );
                auto const& rtensor = M_r_tensor_expr.evalq( q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += ltensor( c1, c2 ) * rtensor( c1, c2 );
                    }
            }
            return res;
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_, mpl::false_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                auto const& ltensor = M_l_tensor_expr.evalq( q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += ltensor( c1, c2 ) * M_r_tensor_expr.evalq( c1, c2, q );
                    }
            }
            return res;
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_, mpl::true_, mpl::false_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                auto const& rtensor = M_r_tensor_expr.evalq( q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        res += M_l_tensor_expr.evalq( c1, c2, q ) * rtensor( c1, c2 );
                    }
            }
            return res;
        }

        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::bool_<true> ) const
        {
            return evalq( cc1, cc2, q, mpl::bool_<l_is_terminal>(), mpl::bool_<true>() );
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::false_, mpl::true_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        value_type val = M_l_tensor_expr.evalq( c1, c2, q );
                        res += val * val;
                    }
            }
            return res;
        }
        value_type
        evalq( uint16_type cc1, uint16_type cc2, uint16_type q, mpl::true_, mpl::true_ ) const
        {
            value_type res = 0;
            if ( Type == 1 )
            {
                auto const& ltensor = M_l_tensor_expr.evalq( q );
                for ( uint16_type c2 = 0; c2 < left_shape::N; ++c2 )
                    for ( uint16_type c1 = 0; c1 < left_shape::M; ++c1 )
                    {
                        value_type val = ltensor( c1, c2 );
                        res += val * val;
                    }
            }
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

/**
 * \brief symetric part of a matricial expression
 */
template <typename ExprL, typename ExprR>
inline Expr<Product<ExprL, ExprR, 1, NONE>>
inner( ExprL l, ExprR r )
{
    typedef Product<ExprL, ExprR, 1, NONE> product_t;
    return Expr<product_t>( product_t( l, r ) );
}

template <typename ExprL, typename ExprR, int Props>
inline Expr<Product<ExprL, ExprR, 1, Props>>
inner( ExprL l, ExprR r, mpl::int_<Props> )
{
    typedef Product<ExprL, ExprR, 1, Props> product_t;
    return Expr<product_t>( product_t( l, r ) );
}

/**
 * \brief symetric part of a matricial expression
 */
template <typename ExprL>
inline Expr<Product<ExprL, ExprL, 1, InnerProperties::IS_SAME>>
inner( ExprL l )
{
    typedef Product<ExprL, ExprL, 1, InnerProperties::IS_SAME> product_t;
    return Expr<product_t>( product_t( l, l ) );
}

template <typename ExprL, int Props>
inline Expr<Product<ExprL, ExprL, 1, InnerProperties::IS_SAME | Props>>
inner( ExprL l, mpl::int_<Props> )
{
    typedef Product<ExprL, ExprL, 1, InnerProperties::IS_SAME | Props> product_t;
    return Expr<product_t>( product_t( l, l ) );
}

} // namespace vf

} // namespace Feel
#endif /* FEELPP_FEELVF_INNER_HPP */
