/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-20

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
   \file matvec.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-20
 */
#ifndef FEELPP_VF_MATVEC_H
#define FEELPP_VF_MATVEC_H 1

#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/operations.hpp>

#include <type_traits>
#include <boost/mp11.hpp>
#include <utility>
#include <variant>

namespace Feel
{
namespace vf
{

/// \cond detail
namespace detail
{

/**
 * \class Mat
 * \brief class that represents a matrix in the language
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int M, int N, typename MatrixExpr>
class Mat
{
public:

    using tuple_expr_type = MatrixExpr;
    using first_expression_type = typename std::decay_t<decltype(hana::at_c<0>( tuple_expr_type{} ))>;
    static const int nExpr = std::decay_t<decltype(hana::size(tuple_expr_type{}))>::value;

    BOOST_MPL_ASSERT_MSG( ( M*N == nExpr ),
                          INVALID_MATRIX_SIZE,
                          ( mpl::int_<M>, mpl::int_<N>, mpl::int_<M*N>,
                            mpl::int_<nExpr> ) );

    /** @name Typedefs
     */
    //@{
    struct FunctorsVariadicExpr
    {
        struct GetContextExpr
        {
            template <typename R,typename T>
            constexpr auto operator()(R /*const&*/ r, T const& t) const
                {
                    return hana::integral_constant<size_type, r.value | std::decay_t<decltype(t)>::context >{};
                }
        };

        template<typename Funct>
        struct HasTestFunction
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<bool, T1::value || T2::template HasTestFunction<Funct>::result >{};
                }
        };
        template<typename Funct>
        struct HasTrialFunction
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<bool, T1::value || T2::template HasTrialFunction<Funct>::result >{};
                }
        };
        template<typename Funct>
        struct HasTestBasis
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<bool, T1::value || T2::template has_test_basis<Funct>::result >{};
                }
        };
        template<typename Funct>
        struct HasTrialBasis
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<bool, T1::value || T2::template has_trial_basis<Funct>::result >{};
                }
        };
    };

    static const size_type context = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                        hana::integral_constant<size_type, 0>{},
                                                                        typename FunctorsVariadicExpr::GetContextExpr{}
                                                                        ) )>::value;
    static inline const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                      hana::integral_constant<bool,false>{},
                                                                      typename FunctorsVariadicExpr::template HasTestFunction<Func>{}
                                                                      ) )>::value;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                      hana::integral_constant<bool,false>{},
                                                                      typename FunctorsVariadicExpr::template HasTrialFunction<Func>{}
                                                                      ) )>::value;
    };
    template<typename Func>
    static inline const bool has_test_basis = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                          hana::integral_constant<bool,false>{},
                                                                          typename FunctorsVariadicExpr::template HasTestBasis<Func>{}
                                                                          ) )>::value;
    template<typename Func>
    static inline const bool has_trial_basis = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                           hana::integral_constant<bool,false>{},
                                                                           typename FunctorsVariadicExpr::template HasTrialBasis<Func>{}
                                                                           ) )>::value;

    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    typedef MatrixExpr expression_matrix_type;
    typedef Mat<M, N, expression_matrix_type> this_type;

    static inline const uint16_type matrix_size1 = M;
    static inline const uint16_type matrix_size2 = N;
    static inline const uint16_type matrix_size  = M*N;

    //typedef double value_type;
    using value_type = typename first_expression_type::value_type;
    using evaluate_type = Eigen::Matrix<value_type,matrix_size1,matrix_size2>;


    template<typename... TheExpr>
    struct Lambda
    {
        template <typename T>
        using TransformLambdaExpr = typename T::template Lambda<TheExpr...>::type;
#if 0
        // utility to convert tuple (hana, std, ... )
        template <template <typename...> class C, typename Tuple> struct RebindImpl;
        template <template <typename...> class C, typename ... Ts>
        struct RebindImpl<C, hana::tuple<Ts...>>{
            using type = C<Ts...>;
        };
        template <template <typename...> class C, typename ... Ts>
        struct RebindImpl<C, std::tuple<Ts...>>{
            using type = C<Ts...>;
        };

        // convert hana::tuple<...> to std::tuple<...>
        using expr_tuple_as_std_tuple_type = typename RebindImpl< std::tuple, expression_matrix_type >::type;

        // transfrom std::tuple<..> with Lambda type
        using lambda_expr_tuple_as_std_tuple_type = boost::mp11::mp_transform<TransformLambdaExpr,expr_tuple_as_std_tuple_type>;

        // get type by converting std::tuple<...> to hana::tuple<...>
        using type = Mat<M,N, typename RebindImpl< hana::tuple, lambda_expr_tuple_as_std_tuple_type >::type>;
#else
        using type = Mat<M,N,  boost::mp11::mp_transform<TransformLambdaExpr,expression_matrix_type> >;
#endif
    };

    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  )
        {
            return typename Lambda<TheExpr...>::type( hana::transform( this->expression(), [&e...]( auto const& t) { return t(e...); } ) );
        }


    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Mat( MatrixExpr const& expr )
        :
        M_expr( expr )
    {
    }
    explicit Mat( MatrixExpr && expr )
        :
        M_expr( std::move( expr ) )
    {
    }
    Mat( Mat const& ) = default;
    Mat( Mat&& ) = default;
    Mat& operator=( Mat const& ) = default;
    Mat& operator=( Mat&& ) = default;
    ~Mat() = default;

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    expression_matrix_type const&  expression() const
    {
        return M_expr;
    }
    expression_matrix_type      &  expression()
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

    //! dynamic context
    size_type dynamicContext() const
        {
            size_type res = 0;
            hana::for_each( M_expr, [&res]( auto const& e )
                            {
                                res = res | Feel::vf::dynamicContext( e );
                            } );
            return res;
        }

    //! polynomial order
    uint16_type polynomialOrder() const
        {
            return hana::fold( M_expr, uint16_type(0), []( uint16_type res, auto const& e ) { return std::max( res, e.polynomialOrder() ); });
        }

    //! expression is polynomial?
    bool isPolynomial() const
        {
            return hana::fold( M_expr, true, []( bool res, auto const& e ) { return res && e.isPolynomial(); });
        }

    //! evaluate the expression without context
    evaluate_type evaluate( bool parallel ) const
        {
            evaluate_type res;
            uint16_type k = 0;
            hana::for_each( M_expr, [&parallel,&k,&res]( auto const& e )
                            {
                                uint16_type i = k / res.cols();
                                uint16_type j = k % res.cols();
                                res( i,j ) = e.evaluate(parallel)(0,0);
                                ++k;
                            } );
            return res;
        }

    void setParameterValues( std::map<std::string,value_type> const& mp )
        {
            hana::for_each( M_expr, [&mp]( auto & e ) { e.setParameterValues( mp ); } );
        }
    void updateParameterValues( std::map<std::string,double> & pv ) const
        {
            hana::for_each( M_expr, [&pv]( auto const& e ) { e.updateParameterValues( pv ); } );
        }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            auto newTupleExprs = hana::transform( M_expr, [&se](auto const& t){ return t.applySymbolsExpr( se ); });
            return Mat<M, N, std::decay_t<decltype(newTupleExprs)> >( std::move( newTupleExprs ) );
        }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            bool res = false;
            hana::for_each( M_expr, [&symb,&se,&res]( auto const& e )
                            {
                                if ( res )
                                    return;
                                res = e.hasSymbolDependency( symb, se );
                            } );
            return res;
        }
    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            hana::for_each( M_expr, [&symb,&res,&se]( auto const& e ) { e.dependentSymbols( symb,res,se ); } );
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
        {
            auto newTupleExprs = hana::transform( M_expr, [&diffVariable,&world,&dirLibExpr,&se](auto const& t){ return t.template diff<diffOrder>( diffVariable, world, dirLibExpr, se ); });
            return Mat<M, N, std::decay_t<decltype(newTupleExprs)> >( std::move( newTupleExprs ) );
        }

    //@}

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_type::expression_matrix_type expression_matrix_type;

        using shape = ShapeGeneric<gmc_t<Geo_t>::nDim,M,N>;

        using value_type = expression_type::value_type;

        template <typename T>
        using tensor_expr_type = typename T::template tensor<Geo_t, Basis_i_t, Basis_j_t>;

        struct TensorFactory
        {
            template <typename T>
            constexpr auto operator()(T const& t) const
            {
                static_assert( sizeof( T ) == 0,
                               "TensorFactory requires geometric context to construct component tensors" );
            }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) const
                {
                    using _tensor_type = tensor_expr_type<T>;
                    return _tensor_type( t, geom, fev, feu );
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom, Basis_i_t const& fev ) const
                {
                    using _tensor_type = tensor_expr_type<T>;
                    return _tensor_type( t, geom, fev );
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom ) const
                {
                    using _tensor_type = tensor_expr_type<T>;
                    return _tensor_type( t, geom );
                }

            template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename T, typename... TheArgsType>
            constexpr auto operator()(std::true_type, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                                      T const& t, Geo_t const& geom, const TheArgsType&... theInitArgs ) const
                {
                    using _tensor_type = tensor_expr_type<T>;
                    return _tensor_type( std::true_type{}, exprExpanded, ttse, t, geom, theInitArgs... );
                }
        };
        using tensor_matrix_type = boost::mp11::mp_transform<tensor_expr_type, expression_matrix_type>;

        struct is_zero
        {
            static inline const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_expr( hana::transform( expr.expression(), [&geom,&fev,&feu](auto const& t) { return TensorFactory{}(t,geom,fev,feu); } ) )
            {}

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_expr( hana::transform( expr.expression(), [&geom,&fev](auto const& t) { return TensorFactory{}(t,geom,fev); } ) )
            {}

        tensor( expression_type const& expr,
                Geo_t const& geom )
            :
            M_expr( hana::transform( expr.expression(), [&geom](auto const& t) { return TensorFactory{}(t,geom); } ) )
            {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                expression_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            //M_expr( hana::transform( hana::make_range( hana::int_c<0>, hana::int_c<nExpr> ), [&exprExpanded,ttse,&expr,&geom,&theInitArgs...](auto eId )
            M_expr( hana::transform( hana::unpack( hana::make_range( hana::int_c<0>, hana::int_c<nExpr> ), hana::make_tuple ), [&exprExpanded,&ttse,&expr,&geom,&theInitArgs...](auto eId )
                                     {
                                         return TensorFactory{}( std::true_type{}, hana::at( exprExpanded.expression(), hana::int_c<eId> ), ttse,
                                                                 hana::at( expr.expression(), hana::int_c<eId> ), geom, theInitArgs... );
                                     } ) )
            {}

        template <typename EvalFn, std::size_t... Index>
        static value_type evaluateFlatComponent( tensor_matrix_type const& exprs, uint16_type flatIndex, EvalFn&& eval,
                                                 std::index_sequence<Index...> )
        {
            value_type result( 0 );
            [[maybe_unused]] auto matched =
                ( ( flatIndex == Index ? ( result = eval( hana::at_c<Index>( exprs ) ), true ) : false ) || ... );
            return result;
        }

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            hana::for_each( M_expr, [&geom,&fev,&feu]( auto & e ) { e.update( geom, fev, feu ); } );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            hana::for_each( M_expr, [&geom,&fev]( auto & e ) { e.update( geom, fev ); } );
        }
        void update( Geo_t const& geom )
        {
            hana::for_each( M_expr, [&geom]( auto & e ) { e.update( geom ); } );
        }
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
        {
            hana::for_each( M_expr, [&ctx...]( auto & e ) { e.updateContext( ctx... ); } );
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nExpr> ), [this,&exprExpanded,&ttse,&geom,&theUpdateArgs...]( auto eId )
                            {
                                hana::at( M_expr, hana::int_c<eId> ).update( std::true_type{}, hana::at( exprExpanded.expression(), hana::int_c<eId> ), ttse, geom, theUpdateArgs... );
                            } );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            auto index = uint16_type( c1*expression_type::matrix_size2 + c2 );
            return evaluateFlatComponent( M_expr, index,
                                          [&i,&j,&q]( auto const& e ) { return e.evalijq( i, j, 0, 0, q ); },
                                          std::make_index_sequence<expression_type::matrix_size>{} );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return this->evalijq(i,j,c1,c2,q);
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            auto index = uint16_type( c1*expression_type::matrix_size2 + c2 );
            return evaluateFlatComponent( M_expr, index,
                                          [&i,&q]( auto const& e ) { return e.evaliq( i, 0, 0, q ); },
                                          std::make_index_sequence<expression_type::matrix_size>{} );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                auto index = uint16_type( c1*expression_type::matrix_size2 + c2 );
                return evaluateFlatComponent( M_expr, index,
                                              [&q]( auto const& e ) { return e.evalq( 0, 0, q ); },
                                              std::make_index_sequence<expression_type::matrix_size>{} );
        }
        tensor_matrix_type M_expr;
    };



protected:

private:
    Mat();

    expression_matrix_type M_expr;
};

template<typename VectorExpr>
using Vec = Mat<std::decay_t<decltype(hana::size(VectorExpr{}))>::value,1,VectorExpr>;

template <template <typename...> class Target, typename Tuple>
struct tuple_rebind;

template <template <typename...> class Target, typename... Ts>
struct tuple_rebind<Target, hana::tuple<Ts...>>
{
    using type = Target<Ts...>;
};

template <typename Tuple, template <typename...> class Target>
using tuple_rebind_t = typename tuple_rebind<Target, Tuple>::type;

template <typename T>
concept MatVecScalarComponent = std::is_arithmetic_v<std::remove_cvref_t<T>>;

template <typename T>
concept MatVecComponentExpr = VfExpr<std::remove_cvref_t<T>> || MatVecScalarComponent<T>;

template <typename T>
auto
normalizeMatVecComponent( T&& value )
{
    if constexpr ( VfExpr<std::remove_cvref_t<T>> )
        return std::forward<T>( value );
    else
        return cst( std::forward<T>( value ) );
}

template <typename TupleType, int Cols, int Row, int Col>
decltype(auto)
matTupleComponent( TupleType const& tuple )
{
    constexpr auto flatIndex = Row*Cols + Col;
    return hana::at_c<flatIndex>( tuple );
}

template <std::size_t Count, std::size_t I = 0, typename Tuple, typename Fn>
decltype(auto)
visitTupleAt( Tuple&& tuple, std::size_t flatIndex, Fn&& fn )
{
    if constexpr ( I+1 == Count )
    {
        CHECK( flatIndex == I ) << "invalid tuple component index " << flatIndex;
        return std::forward<Fn>( fn )( hana::at_c<I>( std::forward<Tuple>( tuple ) ) );
    }
    else
    {
        if ( flatIndex == I )
            return std::forward<Fn>( fn )( hana::at_c<I>( std::forward<Tuple>( tuple ) ) );
        return visitTupleAt<Count, I+1>( std::forward<Tuple>( tuple ), flatIndex, std::forward<Fn>( fn ) );
    }
}

template <std::size_t Count, std::size_t I = 0, typename Tuple1, typename Tuple2, typename Fn>
decltype(auto)
visitMatchedTuplesAt( Tuple1&& tuple1, Tuple2&& tuple2, std::size_t flatIndex, Fn&& fn )
{
    if constexpr ( I+1 == Count )
    {
        CHECK( flatIndex == I ) << "invalid matched tuple index " << flatIndex;
        return std::forward<Fn>( fn )( hana::at_c<I>( std::forward<Tuple1>( tuple1 ) ),
                                       hana::at_c<I>( std::forward<Tuple2>( tuple2 ) ) );
    }
    else
    {
        if ( flatIndex == I )
            return std::forward<Fn>( fn )( hana::at_c<I>( std::forward<Tuple1>( tuple1 ) ),
                                           hana::at_c<I>( std::forward<Tuple2>( tuple2 ) ) );
        return visitMatchedTuplesAt<Count, I+1>( std::forward<Tuple1>( tuple1 ),
                                                 std::forward<Tuple2>( tuple2 ),
                                                 flatIndex,
                                                 std::forward<Fn>( fn ) );
    }
}

template <std::size_t Count, std::size_t I = 0, typename Variant, typename Fn>
decltype(auto)
visitVariantAt( Variant&& variant, std::size_t flatIndex, Fn&& fn )
{
    if constexpr ( I+1 == Count )
    {
        CHECK( flatIndex == I ) << "invalid variant index " << flatIndex;
        return std::forward<Fn>( fn )( std::get<I>( std::forward<Variant>( variant ) ) );
    }
    else
    {
        if ( flatIndex == I )
            return std::forward<Fn>( fn )( std::get<I>( std::forward<Variant>( variant ) ) );
        return visitVariantAt<Count, I+1>( std::forward<Variant>( variant ), flatIndex, std::forward<Fn>( fn ) );
    }
}

template <std::size_t Count, std::size_t I = 0, typename Variant, typename Tuple, typename Fn>
decltype(auto)
visitVariantAndTupleAt( Variant&& variant, Tuple&& tuple, std::size_t flatIndex, Fn&& fn )
{
    if constexpr ( I+1 == Count )
    {
        CHECK( flatIndex == I ) << "invalid variant/tuple index " << flatIndex;
        return std::forward<Fn>( fn )( std::get<I>( std::forward<Variant>( variant ) ),
                                       hana::at_c<I>( std::forward<Tuple>( tuple ) ) );
    }
    else
    {
        if ( flatIndex == I )
            return std::forward<Fn>( fn )( std::get<I>( std::forward<Variant>( variant ) ),
                                           hana::at_c<I>( std::forward<Tuple>( tuple ) ) );
        return visitVariantAndTupleAt<Count, I+1>( std::forward<Variant>( variant ),
                                                   std::forward<Tuple>( tuple ),
                                                   flatIndex,
                                                   std::forward<Fn>( fn ) );
    }
}

template <int SharedDim, int RightCols, int Row, int Col, typename LeftTuple, typename RightTuple, std::size_t... SumIndex>
auto
makeMatProductEntry( LeftTuple const& leftTuple, RightTuple const& rightTuple, std::index_sequence<SumIndex...> )
{
    return ( ( matTupleComponent<LeftTuple, SharedDim, Row, static_cast<int>( SumIndex )>( leftTuple ) *
               matTupleComponent<RightTuple, RightCols, static_cast<int>( SumIndex ), Col>( rightTuple ) ) + ... );
}

template <int LeftRows, int SharedDim, int RightCols, typename LeftTuple, typename RightTuple, std::size_t... FlatIndex>
auto
makeMatProductTuple( LeftTuple const& leftTuple, RightTuple const& rightTuple, std::index_sequence<FlatIndex...> )
{
    return hana::make_tuple(
        makeMatProductEntry<SharedDim,
                            RightCols,
                            static_cast<int>( FlatIndex/RightCols ),
                            static_cast<int>( FlatIndex%RightCols )>( leftTuple, rightTuple, std::make_index_sequence<SharedDim>{} )... );
}

} // detail
/// \endcond

template <int M, int N, typename MatrixExpr>
class ComponentsExpr<Expr<detail::Mat<M, N, MatrixExpr>>>
{
public:
    using matrix_expression_type = detail::Mat<M, N, MatrixExpr>;
    using expression_type = Expr<matrix_expression_type>;
    using expression_matrix_type = MatrixExpr;
    using value_type = typename expression_type::value_type;
    using evaluate_type = Eigen::Matrix<value_type,1,1>;
    using this_type = ComponentsExpr<expression_type>;

    static const size_type context = expression_type::context;
    static inline const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = expression_type::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = expression_type::template HasTrialFunction<Func>::result;
    };

    template<typename Func>
    static inline const bool has_test_basis = expression_type::template has_test_basis<Func>;
    template<typename Func>
    static inline const bool has_trial_basis = expression_type::template has_trial_basis<Func>;
    using test_basis = typename expression_type::test_basis;
    using trial_basis = typename expression_type::trial_basis;

    ComponentsExpr()
        :
        M_expr(),
        M_c1( 0 ),
        M_c2( 0 ),
        M_flatIndex( 0 )
    {}

    explicit ComponentsExpr( expression_type const& expr, int c1, int c2 )
        :
        M_expr( expr ),
        M_c1( c1 ),
        M_c2( c2 ),
        M_flatIndex( this->flatIndexFromComponents( c1, c2 ) )
    {}

    size_type dynamicContext() const
    {
        return detail::visitTupleAt<M*N>( M_expr.expression().expression(), M_flatIndex,
                                          []( auto const& e ) { return Feel::vf::dynamicContext( e ); } );
    }

    uint16_type polynomialOrder() const
    {
        return detail::visitTupleAt<M*N>( M_expr.expression().expression(), M_flatIndex,
                                          []( auto const& e ) { return e.polynomialOrder(); } );
    }

    bool isPolynomial() const
    {
        return detail::visitTupleAt<M*N>( M_expr.expression().expression(), M_flatIndex,
                                          []( auto const& e ) { return e.isPolynomial(); } );
    }

    expression_type const& expression() const
    {
        return M_expr;
    }

    std::size_t flatIndex() const
    {
        return M_flatIndex;
    }

    evaluate_type evaluate( bool p ) const
    {
        return evaluate_type::Constant( detail::visitTupleAt<M*N>( M_expr.expression().expression(), M_flatIndex,
                                                                   [p]( auto const& e ) { return e.evaluate( p )( 0, 0 ); } ) );
    }

    void setParameterValues( std::map<std::string,double> const& mp )
    {
        M_expr.setParameterValues( mp );
    }
    void updateParameterValues( std::map<std::string,double> & pv ) const
    {
        detail::visitTupleAt<M*N>( M_expr.expression().expression(), M_flatIndex,
                                   [&pv]( auto const& e ) { e.updateParameterValues( pv ); } );
    }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto newExpr = M_expr.applySymbolsExpr( se );
        using new_expr_type = std::decay_t<decltype(newExpr)>;
        return ComponentsExpr<new_expr_type>( newExpr, M_c1, M_c2 );
    }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        return detail::visitTupleAt<M*N>( M_expr.expression().expression(), M_flatIndex,
                                          [&symb,&se]( auto const& e ) { return e.hasSymbolDependency( symb, se ); } );
    }
    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
    {
        detail::visitTupleAt<M*N>( M_expr.expression().expression(), M_flatIndex,
                                   [&symb,&res,&se]( auto const& e ) { e.dependentSymbols( symb, res, se ); } );
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        auto theDiffExpr = M_expr.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
        using new_expr_type = std::decay_t<decltype(theDiffExpr)>;
        return ComponentsExpr<new_expr_type>( theDiffExpr, M_c1, M_c2 );
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        template <typename T>
        using tensor_expr_type = typename std::remove_cvref_t<T>::template tensor<Geo_t, Basis_i_t, Basis_j_t>;

        using tensor_tuple_type = boost::mp11::mp_transform<tensor_expr_type, expression_matrix_type>;
        using tensor_variant_type = detail::tuple_rebind_t<tensor_tuple_type, std::variant>;
        using value_type = typename this_type::value_type;
        using key_type = key_t<Geo_t>;
        using gmc_type = typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type;
        using shape = Shape<gmc_type::NDim, Scalar, false>;

        template <class Args> struct sig
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = false;
        };

        template <std::size_t Index = 0>
        static tensor_variant_type makeTensorVariantImpl( this_type const& expr, std::size_t flatIndex,
                                                          Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            if constexpr ( Index + 1 == M*N )
            {
                CHECK( flatIndex == Index ) << "invalid tensor component index " << flatIndex;
                auto const& e = hana::at_c<Index>( expr.expression().expression().expression() );
                return tensor_variant_type( std::in_place_index<Index>, e, geom, fev, feu );
            }
            else
            {
                if ( flatIndex == Index )
                {
                    auto const& e = hana::at_c<Index>( expr.expression().expression().expression() );
                    return tensor_variant_type( std::in_place_index<Index>, e, geom, fev, feu );
                }
                return makeTensorVariantImpl<Index + 1>( expr, flatIndex, geom, fev, feu );
            }
        }

        template <std::size_t Index = 0>
        static tensor_variant_type makeTensorVariantImpl( this_type const& expr, std::size_t flatIndex,
                                                          Geo_t const& geom, Basis_i_t const& fev )
        {
            if constexpr ( Index + 1 == M*N )
            {
                CHECK( flatIndex == Index ) << "invalid tensor component index " << flatIndex;
                auto const& e = hana::at_c<Index>( expr.expression().expression().expression() );
                return tensor_variant_type( std::in_place_index<Index>, e, geom, fev );
            }
            else
            {
                if ( flatIndex == Index )
                {
                    auto const& e = hana::at_c<Index>( expr.expression().expression().expression() );
                    return tensor_variant_type( std::in_place_index<Index>, e, geom, fev );
                }
                return makeTensorVariantImpl<Index + 1>( expr, flatIndex, geom, fev );
            }
        }

        template <std::size_t Index = 0>
        static tensor_variant_type makeTensorVariantImpl( this_type const& expr, std::size_t flatIndex,
                                                          Geo_t const& geom )
        {
            if constexpr ( Index + 1 == M*N )
            {
                CHECK( flatIndex == Index ) << "invalid tensor component index " << flatIndex;
                auto const& e = hana::at_c<Index>( expr.expression().expression().expression() );
                return tensor_variant_type( std::in_place_index<Index>, e, geom );
            }
            else
            {
                if ( flatIndex == Index )
                {
                    auto const& e = hana::at_c<Index>( expr.expression().expression().expression() );
                    return tensor_variant_type( std::in_place_index<Index>, e, geom );
                }
                return makeTensorVariantImpl<Index + 1>( expr, flatIndex, geom );
            }
        }

        template <std::size_t Index = 0, typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        static tensor_variant_type makeExpandedTensorVariantImpl( TheExprExpandedType const& exprExpanded, this_type const& expr,
                                                                  std::size_t flatIndex, TupleTensorSymbolsExprType& ttse,
                                                                  Geo_t const& geom, const TheArgsType&... theInitArgs )
        {
            if constexpr ( Index + 1 == M*N )
            {
                CHECK( flatIndex == Index ) << "invalid expanded tensor component index " << flatIndex;
                auto const& expandedComponent = hana::at_c<Index>( exprExpanded.expression().expression().expression() );
                auto const& originalComponent = hana::at_c<Index>( expr.expression().expression().expression() );
                return tensor_variant_type( std::in_place_index<Index>, std::true_type{}, expandedComponent, ttse,
                                            originalComponent, geom, theInitArgs... );
            }
            else
            {
                if ( flatIndex == Index )
                {
                    auto const& expandedComponent = hana::at_c<Index>( exprExpanded.expression().expression().expression() );
                    auto const& originalComponent = hana::at_c<Index>( expr.expression().expression().expression() );
                    return tensor_variant_type( std::in_place_index<Index>, std::true_type{}, expandedComponent, ttse,
                                                originalComponent, geom, theInitArgs... );
                }
                return makeExpandedTensorVariantImpl<Index + 1>( exprExpanded, expr, flatIndex, ttse, geom, theInitArgs... );
            }
        }

        static tensor_variant_type makeTensorVariant( this_type const& expr, std::size_t flatIndex,
                                                      Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            return makeTensorVariantImpl( expr, flatIndex, geom, fev, feu );
        }

        static tensor_variant_type makeTensorVariant( this_type const& expr, std::size_t flatIndex,
                                                      Geo_t const& geom, Basis_i_t const& fev )
        {
            return makeTensorVariantImpl( expr, flatIndex, geom, fev );
        }

        static tensor_variant_type makeTensorVariant( this_type const& expr, std::size_t flatIndex,
                                                      Geo_t const& geom )
        {
            return makeTensorVariantImpl( expr, flatIndex, geom );
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        static tensor_variant_type makeExpandedTensorVariant( TheExprExpandedType const& exprExpanded, this_type const& expr,
                                                              std::size_t flatIndex, TupleTensorSymbolsExprType & ttse,
                                                              Geo_t const& geom, const TheArgsType&... theInitArgs )
        {
            return makeExpandedTensorVariantImpl( exprExpanded, expr, flatIndex, ttse, geom, theInitArgs... );
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( makeTensorVariant( expr, expr.flatIndex(), geom, fev, feu ) ),
            M_flatIndex( expr.flatIndex() )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( makeTensorVariant( expr, expr.flatIndex(), geom, fev ) ),
            M_flatIndex( expr.flatIndex() )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( makeTensorVariant( expr, expr.flatIndex(), geom ) ),
            M_flatIndex( expr.flatIndex() )
        {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_tensor_expr( makeExpandedTensorVariant( exprExpanded, expr, exprExpanded.flatIndex(), ttse, geom, theInitArgs... ) ),
            M_flatIndex( exprExpanded.flatIndex() )
        {}

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            detail::visitVariantAt<M*N>( M_tensor_expr, M_flatIndex, [&geom,&fev,&feu]( auto& e ) { e.update( geom, fev, feu ); } );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            detail::visitVariantAt<M*N>( M_tensor_expr, M_flatIndex, [&geom,&fev]( auto& e ) { e.update( geom, fev ); } );
        }
        void update( Geo_t const& geom )
        {
            detail::visitVariantAt<M*N>( M_tensor_expr, M_flatIndex, [&geom]( auto& e ) { e.update( geom ); } );
        }
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
        {
            detail::visitVariantAt<M*N>( M_tensor_expr, M_flatIndex, [&ctx...]( auto& e ) { e.updateContext( ctx... ); } );
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            detail::visitVariantAndTupleAt<M*N>( M_tensor_expr, exprExpanded.expression().expression().expression(), M_flatIndex,
                                                 [&ttse,&geom,&theUpdateArgs...]( auto& tensorExpr, auto const& componentExpr )
                                                 {
                                                     tensorExpr.update( std::true_type{}, componentExpr, ttse, geom, theUpdateArgs... );
                                                 } );
        }

        value_type evalij( uint16_type i, uint16_type j ) const
        {
            return detail::visitVariantAt<M*N>( M_tensor_expr, M_flatIndex, [&i,&j]( auto const& e ) { return e.evalij( i, j ); } );
        }

        value_type evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return detail::visitVariantAt<M*N>( M_tensor_expr, M_flatIndex,
                                                [&i,&j,&q]( auto const& e ) { return e.evalijq( i, j, 0, 0, q ); } );
        }

        template<int PatternContext>
        value_type evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q,
                            mpl::int_<PatternContext> ) const
        {
            return detail::visitVariantAt<M*N>( M_tensor_expr, M_flatIndex,
                                                [&i,&j,&q]( auto const& e ) { return e.evalijq( i, j, 0, 0, q, mpl::int_<PatternContext>() ); } );
        }

        value_type evaliq( uint16_type i, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return detail::visitVariantAt<M*N>( M_tensor_expr, M_flatIndex,
                                                [&i,&q]( auto const& e ) { return e.evaliq( i, 0, 0, q ); } );
        }

        value_type evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return detail::visitVariantAt<M*N>( M_tensor_expr, M_flatIndex,
                                                [&q]( auto const& e ) { return e.evalq( 0, 0, q ); } );
        }

        tensor_variant_type M_tensor_expr;
        std::size_t M_flatIndex;
    };

private:
    std::size_t flatIndexFromComponents( int c1, int c2 ) const
    {
        CHECK( c1 >= 0 && c1 < M ) << "row component index out of range: " << c1;
        CHECK( c2 >= 0 && c2 < N ) << "column component index out of range: " << c2;
        return static_cast<std::size_t>( c1*N + c2 );
    }

    expression_type M_expr;
    int M_c1, M_c2;
    std::size_t M_flatIndex;
};

template <int M, int N, typename MatrixExpr, int C1, int C2>
class StaticComponentsExpr<Expr<detail::Mat<M, N, MatrixExpr>>, C1, C2>
{
public:
    static_assert( C1 >= 0 && C1 < M, "Static matrix row component index is out of range" );
    static_assert( C2 >= 0 && C2 < N, "Static matrix column component index is out of range" );

    using matrix_expression_type = detail::Mat<M, N, MatrixExpr>;
    using expression_type = Expr<matrix_expression_type>;
    using expression_matrix_type = MatrixExpr;
    using component_expression_type = std::remove_cvref_t<decltype( detail::matTupleComponent<expression_matrix_type, N, C1, C2>( std::declval<expression_matrix_type const&>() ) )>;
    using value_type = typename expression_type::value_type;
    using evaluate_type = Eigen::Matrix<value_type,1,1>;
    using this_type = StaticComponentsExpr<expression_type, C1, C2>;

    static const size_type context = expression_type::context;
    static inline const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = expression_type::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = expression_type::template HasTrialFunction<Func>::result;
    };

    template<typename Func>
    static inline const bool has_test_basis = expression_type::template has_test_basis<Func>;
    template<typename Func>
    static inline const bool has_trial_basis = expression_type::template has_trial_basis<Func>;
    using test_basis = typename expression_type::test_basis;
    using trial_basis = typename expression_type::trial_basis;

    StaticComponentsExpr()
        :
        M_expr()
    {}

    explicit StaticComponentsExpr( expression_type const& expr )
        :
        M_expr( expr )
    {}

    size_type dynamicContext() const
    {
        return Feel::vf::dynamicContext( componentExpr( M_expr ) );
    }

    uint16_type polynomialOrder() const
    {
        return componentExpr( M_expr ).polynomialOrder();
    }

    bool isPolynomial() const
    {
        return componentExpr( M_expr ).isPolynomial();
    }

    expression_type const& expression() const
    {
        return M_expr;
    }

    evaluate_type evaluate( bool p ) const
    {
        return evaluate_type::Constant( componentExpr( M_expr ).evaluate( p )( 0, 0 ) );
    }

    void setParameterValues( std::map<std::string,double> const& mp )
    {
        M_expr.setParameterValues( mp );
    }
    void updateParameterValues( std::map<std::string,double> & pv ) const
    {
        componentExpr( M_expr ).updateParameterValues( pv );
    }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto newExpr = M_expr.applySymbolsExpr( se );
        using new_expr_type = std::decay_t<decltype(newExpr)>;
        return StaticComponentsExpr<new_expr_type, C1, C2>( newExpr );
    }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        return componentExpr( M_expr ).hasSymbolDependency( symb, se );
    }
    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
    {
        componentExpr( M_expr ).dependentSymbols( symb, res, se );
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        auto theDiffExpr = M_expr.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
        using new_expr_type = std::decay_t<decltype(theDiffExpr)>;
        return StaticComponentsExpr<new_expr_type, C1, C2>( theDiffExpr );
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using tensor_expr_type = typename component_expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using value_type = typename this_type::value_type;
        using key_type = key_t<Geo_t>;
        using gmc_type = typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type;
        using shape = Shape<gmc_type::NDim, Scalar, false>;

        template <class Args> struct sig
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( componentExpr( expr.expression() ), geom, fev, feu )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( componentExpr( expr.expression() ), geom, fev )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( componentExpr( expr.expression() ), geom )
        {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_tensor_expr( std::true_type{}, componentExpr( exprExpanded.expression() ), ttse, componentExpr( expr.expression() ), geom, theInitArgs... )
            {}

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
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
        {
            M_tensor_expr.updateContext( ctx... );
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            M_tensor_expr.update( std::true_type{}, componentExpr( exprExpanded.expression() ), ttse, geom, theUpdateArgs... );
        }

        value_type evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }

        value_type evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_tensor_expr.evalijq( i, j, 0, 0, q );
        }

        template<int PatternContext>
        value_type evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q,
                            mpl::int_<PatternContext> ) const
        {
            return M_tensor_expr.evalijq( i, j, 0, 0, q, mpl::int_<PatternContext>() );
        }

        value_type evaliq( uint16_type i, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_tensor_expr.evaliq( i, 0, 0, q );
        }

        value_type evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_tensor_expr.evalq( 0, 0, q );
        }

        tensor_expr_type M_tensor_expr;
    };

private:
    static decltype(auto) componentExpr( expression_type const& expr )
    {
        return detail::matTupleComponent<expression_matrix_type, N, C1, C2>( expr.expression().expression() );
    }

    expression_type M_expr;
};

/**
 * \brief vector definition
 */
template<typename Expr1,typename ... ExprT>
inline
auto
vec( Expr1&& expr1, ExprT&&... expr2 )
requires detail::MatVecComponentExpr<Expr1> && ( detail::MatVecComponentExpr<ExprT> && ... )
{
    auto normalizedExpr = hana::make_tuple( detail::normalizeMatVecComponent( std::forward<Expr1>( expr1 ) ),
                                            detail::normalizeMatVecComponent( std::forward<ExprT>( expr2 ) )... );
    using expr_t = vf::detail::Vec<std::decay_t<decltype(normalizedExpr)>>;
    return Expr<expr_t>( expr_t( std::move( normalizedExpr ) ) );
}


/**
 * \brief matrix definition
 */
template<int M, int N, typename Expr1,typename ... ExprT>
inline
auto
mat( Expr1&& expr1, ExprT&&... expr2 )
requires detail::MatVecComponentExpr<Expr1> &&
         ( detail::MatVecComponentExpr<ExprT> && ... ) &&
         ( 1 + sizeof...( ExprT ) == M*N )
{
    auto normalizedExpr = hana::make_tuple( detail::normalizeMatVecComponent( std::forward<Expr1>( expr1 ) ),
                                            detail::normalizeMatVecComponent( std::forward<ExprT>( expr2 ) )... );
    using expr_t = vf::detail::Mat<M, N, std::decay_t<decltype(normalizedExpr)> >;
    return Expr<expr_t>( expr_t( std::move( normalizedExpr ) ) );
}

template <int M, int K, typename LeftTuple, int N, typename RightTuple>
inline
auto
operator*( Expr<vf::detail::Mat<M, K, LeftTuple>> const& left, Expr<vf::detail::Mat<K, N, RightTuple>> const& right )
{
    auto const& leftTuple = left.expression().expression();
    auto const& rightTuple = right.expression().expression();
    auto productTuple = vf::detail::makeMatProductTuple<M, K, N>( leftTuple, rightTuple, std::make_index_sequence<M*N>{} );

    using product_tuple_type = std::decay_t<decltype( productTuple )>;
    using expr_t = vf::detail::Mat<M, N, product_tuple_type>;

    return Expr<expr_t>( expr_t( std::move( productTuple ) ) );
}


} // vf
} // Feel
#endif /* __VFVec_H */
