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

#include <boost/mp11.hpp>

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
    static const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                      hana::integral_constant<bool,false>{},
                                                                      typename FunctorsVariadicExpr::template HasTestFunction<Func>{}
                                                                      ) )>::value;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                      hana::integral_constant<bool,false>{},
                                                                      typename FunctorsVariadicExpr::template HasTrialFunction<Func>{}
                                                                      ) )>::value;
    };
    template<typename Func>
    static const bool has_test_basis = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                          hana::integral_constant<bool,false>{},
                                                                          typename FunctorsVariadicExpr::template HasTestBasis<Func>{}
                                                                          ) )>::value;
    template<typename Func>
    static const bool has_trial_basis = std::decay_t<decltype( hana::fold( tuple_expr_type{},
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
        M_expr( expr )
    {
    }
    Mat( Mat const & expr )
        :
        M_expr( expr.M_expr )
    {
    }
    ~Mat()
    {}

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
    evaluate_type evaluate( bool parallel, worldcomm_ptr_t const& worldcomm ) const
        {
            evaluate_type res;
            uint16_type k = 0;
            hana::for_each( M_expr, [&parallel,&worldcomm,&k,&res]( auto const& e )
                            {
                                uint16_type i = k / res.cols();
                                uint16_type j = k % res.cols();
                                res( i,j ) = e.evaluate(parallel,worldcomm)(0,0);
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

        struct TransformExprToTensor
        {
            template <typename T>
            struct apply {
                using type = typename T::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
            };
            template <typename T>
            constexpr auto operator()(T const& t) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    //return _tensor_type( t, Geo_t{} );
                    std::unique_ptr<_tensor_type> result;
                    return *result;
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    return _tensor_type( t, geom, fev, feu );
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom, Basis_i_t const& fev ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    return _tensor_type( t, geom, fev );
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    return _tensor_type( t, geom );
                }

            template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename T, typename... TheArgsType>
            constexpr auto operator()(std::true_type, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                                      T const& t, Geo_t const& geom, const TheArgsType&... theInitArgs ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    return _tensor_type( std::true_type{}, exprExpanded, ttse, t, geom, theInitArgs... );
                }
        };
        using tensor_matrix_type = std::decay_t<decltype( hana::transform( expression_matrix_type{}, TransformExprToTensor{} ) ) >;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_expr( hana::transform( expr.expression(), [&geom,&fev,&feu](auto const& t) { return TransformExprToTensor{}(t,geom,fev,feu); } ) )
            {}

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_expr( hana::transform( expr.expression(), [&geom,&fev](auto const& t) { return TransformExprToTensor{}(t,geom,fev); } ) )
            {}

        tensor( expression_type const& expr,
                Geo_t const& geom )
            :
            M_expr( hana::transform( expr.expression(), [&geom](auto const& t) { return TransformExprToTensor{}(t,geom); } ) )
            {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                expression_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            //M_expr( hana::transform( hana::make_range( hana::int_c<0>, hana::int_c<nExpr> ), [&exprExpanded,ttse,&expr,&geom,&theInitArgs...](auto eId )
            M_expr( hana::transform( hana::unpack( hana::make_range( hana::int_c<0>, hana::int_c<nExpr> ), hana::make_tuple ), [&exprExpanded,&ttse,&expr,&geom,&theInitArgs...](auto eId )
                                     {
                                         return TransformExprToTensor{}(std::true_type{},hana::at( exprExpanded.expression(), hana::int_c<eId> ), ttse,
                                                                        hana::at( expr.expression(), hana::int_c<eId> ), geom, theInitArgs...);
                                     } ) )
            {}

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
            uint16_type index = c1*expression_type::matrix_size2+c2;
            uint16_type k = 0;
            value_type res(0);
            hana::for_each( M_expr, [&i,&j,&q,&index,&k,&res]( auto const& e )
                            {
                                if ( k == index )
                                    res = e.evalijq( i,j,0,0,q );
                                ++k;
                            } );
            return res;
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
            uint16_type index = c1*expression_type::matrix_size2+c2;
            uint16_type k = 0;
            value_type res(0);
            hana::for_each( M_expr, [&i,&q,&index,&k,&res]( auto const& e )
                            {
                                if ( k == index )
                                    res = e.evaliq( i,0,0,q );
                                ++k;
                            } );
            return res;
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                uint16_type index = c1*expression_type::matrix_size2+c2;
                uint16_type k = 0;
                value_type res(0);
                hana::for_each( M_expr, [&q,&index,&k,&res]( auto const& e )
                                {
                                    if ( k == index )
                                        res = e.evalq( 0,0,q );
                                    ++k;
                                } );
                return res;
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

} // detail
/// \endcond

/**
 * \brief vector definition
 */
template<typename Expr1,typename ... ExprT>
inline
auto
vec( Expr1 const& expr1, const ExprT&... expr2 )
{
    using expr_t = vf::detail::Vec< hana::tuple<Expr1,ExprT...> >;
    return Expr<expr_t>( expr_t( hana::make_tuple( expr1, expr2... ) ) );
}


/**
 * \brief matrix definition
 */
template<int M, int N, typename Expr1,typename ... ExprT>
inline
Expr<vf::detail::Mat<M, N, hana::tuple<Expr1,ExprT...> > >
mat( Expr1 const& expr1, const ExprT&... expr2 )
{
    using expr_t = vf::detail::Mat<M, N, hana::tuple<Expr1,ExprT...> > ;
    return Expr<expr_t>( expr_t( hana::make_tuple( expr1, expr2... ) ) );
}


} // vf
} // Feel
#endif /* __VFVec_H */
