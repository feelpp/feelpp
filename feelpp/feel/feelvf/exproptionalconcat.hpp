/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2020-03-25

  Copyright (C) 2020 Feel++ Consortium

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
#ifndef FEELPP_MODELS_VF_ExprOptionalConcat_H
#define FEELPP_MODELS_VF_ExprOptionalConcat_H 1

namespace Feel
{
namespace vf
{

template <typename TupleExprType>
class ExprOptionalConcat : public Feel::vf::ExprDynamicBase
{
    struct TransformOptionalExpr
    {
        template <typename T>
        constexpr auto operator()(T const& t) const
            {
                return std::vector<T>{};
            }
    };

     struct TransformOptionalExprInverse
     {
         template <typename T>
         typename std::decay_t<T>::value_type operator()( T && t) const
             {
                 return t.front();
             }
     };
public :
    using this_type = ExprOptionalConcat<TupleExprType>;

    using tuple_expr_type = TupleExprType;

    static constexpr int nExpr = std::decay_t<decltype(hana::size(tuple_expr_type{}))>::value;

    using tuple_optional_expr_type =  std::decay_t<decltype( hana::transform( tuple_expr_type{}, TransformOptionalExpr{} ) )>;

    using first_expression_type = typename std::decay_t<decltype(hana::at_c<0>( tuple_expr_type{} ))>;

    static const size_type context = vm::DYNAMIC|vm::DYNAMIC_BASIS_FUNCTION;
    static const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                      hana::integral_constant<bool,false>{},
                                                                      []( auto const& res, auto const& e ) { hana::integral_constant<bool, res.value || std::decay_t<decltype(e)>::template HasTestFunction<Func>::result >{}; }
                                                                      ) )>::value;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                      hana::integral_constant<bool,false>{},
                                                                      []( auto const& res, auto const& e ) { hana::integral_constant<bool, res.value || std::decay_t<decltype(e)>::template HasTrialFunction<Func>::result >{}; }
                                                                      ) )>::value;
    };

    template<typename Func>
    static const bool has_test_basis = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                          hana::integral_constant<bool,false>{},
                                                                          []( auto const& res, auto const& e ) { hana::integral_constant<bool, res.value || std::decay_t<decltype(e)>::template has_test_basis<Func>::result >{}; }
                                                                          ) )>::value;
    template<typename Func>
    static const bool has_trial_basis = std::decay_t<decltype( hana::fold( tuple_expr_type{},
                                                                           hana::integral_constant<bool,false>{},
                                                                           []( auto const& res, auto const& e ) { hana::integral_constant<bool, res.value || std::decay_t<decltype(e)>::template has_trial_basis<Func>::result >{}; }
                                                                           ) )>::value;

    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;
    typedef typename first_expression_type::value_type value_type;
    typedef typename first_expression_type::evaluate_type evaluate_type;

    ExprOptionalConcat() = default;
    ExprOptionalConcat( ExprOptionalConcat const& ) = default;
    ExprOptionalConcat( ExprOptionalConcat && ) = default;

    template <typename TheArgTupleExprType>
    ExprOptionalConcat( std::true_type, TheArgTupleExprType && t )
        :
        M_exprs( std::forward<TheArgTupleExprType>( t ) )
        {}

    //! polynomial order
    uint16_type polynomialOrder() const
        {
            uint16_type res = 0;
            hana::for_each( M_exprs, [&res]( auto const& e )
                            {
                                for ( auto const& e2 : e )
                                    res = std::max( res, e2.polynomialOrder() );
                            });
            return res;
        }

    //! expression is polynomial?
    bool isPolynomial() const
        {
            bool res = true;
            hana::for_each( M_exprs, [&res]( auto const& e )
                            {
                                for ( auto const& e2 : e )
                                    res = res && e2.isPolynomial();
                            });
            return res;
        }

    size_type dynamicContext() const
        {
            size_type res = 0;
            hana::for_each( M_exprs, [&res]( auto const& e )
                            {
                                for ( auto const& e2 : e )
                                    res = res | Feel::vf::dynamicContext( e2 );
                            });
            return res;
        }

    tuple_optional_expr_type const& tupleExpressions() const { return M_exprs; }

    bool hasExpr() const
        {
            bool res = false;
            hana::for_each( M_exprs, [&res]( auto const& e )
                            {
                                if ( !e.empty() )
                                    res = true;
                            });
            return res;
        }

    template <int ExprId>
    void set( std::decay_t<decltype(hana::at_c<ExprId>(tuple_expr_type{}))> const& theExpr )
        {
            hana::at_c<ExprId>( M_exprs ).clear();
            hana::at_c<ExprId>( M_exprs ).push_back( theExpr );
        }

    template< typename TheExprType>
    void add( TheExprType const& theExpr )
        {
            bool applyAdd = false;
            hana::for_each( M_exprs, [&theExpr,&applyAdd]( auto & e )
                            {
                                if constexpr ( std::is_same_v<TheExprType, typename std::decay_t<decltype(e)>::value_type> )
                                    {
                                        if ( !applyAdd )
                                        {
                                            e.push_back( theExpr );
                                            applyAdd = true;
                                        }
                                    }
                            });
            CHECK( applyAdd ) << "cannot add the expression, the type is not register";
        }

    void setParameterValues( std::map<std::string,double> const& mp )
        {
            hana::for_each( M_exprs, [&mp]( auto & e )
                            {
                                for ( auto & e2 : e )
                                    e2.setParameterValues( mp );
                            });
        }
    void updateParameterValues( std::map<std::string,double> & pv ) const
        {
            hana::for_each( M_exprs, [&pv]( auto & e )
                            {
                                for ( auto const& e2 : e )
                                    e2.updateParameterValues( pv );
                            });
        }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            auto newTuple = hana::transform( M_exprs, [&se]( auto const& e ) {
                    using _new_expr_type = std::decay_t<decltype(e.front().applySymbolsExpr(se))>;
                    std::vector<_new_expr_type> curentRes;
                    curentRes.reserve(e.size());
                    for ( auto & e2 : e )
                        curentRes.push_back( e2.applySymbolsExpr(se) );
                    return curentRes;
                });

            using new_tuple_optional_expr_type = std::decay_t<decltype( hana::transform( newTuple, TransformOptionalExprInverse{} ) )>;
            return ExprOptionalConcat<new_tuple_optional_expr_type>( std::true_type{},std::move( newTuple ) );
        }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            bool res = false;
            hana::for_each( M_exprs, [&res,&symb,&se]( auto const& e )
                            {
                                if ( res )
                                    return;
                                for ( auto const& e2 : e )
                                {
                                    res = e2.hasSymbolDependency(symb,se);
                                    if ( res )
                                        break;
                                }
                            });
            return res;
        }

    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            hana::for_each( M_exprs, [&symb,&res,&se]( auto const& e )
                            {
                                for ( auto const& e2 : e )
                                    e2.dependentSymbols(symb,res,se);
                            });
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
        {
            auto newTuple = hana::transform( M_exprs, [&diffVariable,&world,&dirLibExpr,&se]( auto const& e ) {
                    using _new_expr_type = std::decay_t<decltype(e.front().template diff<diffOrder>(diffVariable,world,dirLibExpr,se))>;
                    std::vector<_new_expr_type> curentRes;
                    curentRes.reserve(e.size());
                    for ( auto & e2 : e )
                        curentRes.push_back( e2.template diff<diffOrder>(diffVariable,world,dirLibExpr,se) );
                    return curentRes;
                });
            using new_tuple_optional_expr_type = std::decay_t<decltype( hana::transform( newTuple, TransformOptionalExprInverse{} ) )>;
            return ExprOptionalConcat<new_tuple_optional_expr_type>( std::true_type{},std::move( newTuple ) );
        }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        struct TransformExprToTensor
        {
            template <typename T>
            struct apply {
                using type = typename T::value_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
            };

            template <typename T>
            constexpr auto operator()(T const& t) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    return std::vector<_tensor_type>{};
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    std::vector<_tensor_type> res;
                    res.reserve( t.size() );
                    for ( auto const& e : t )
                        res.push_back( _tensor_type( e,geom,fev,feu ) );
                    return res;
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom, Basis_i_t const& fev ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    std::vector<_tensor_type> res;
                    res.reserve( t.size() );
                    for ( auto const& e : t )
                        res.push_back( _tensor_type( e,geom,fev ) );
                    return res;
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    std::vector<_tensor_type> res;
                    res.reserve( t.size() );
                    for ( auto const& e : t )
                        res.push_back( _tensor_type( e,geom ) );
                    return res;
                }
            template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename T, typename... TheArgsType>
            constexpr auto operator()(std::true_type, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                                      T const& t, Geo_t const& geom, const TheArgsType&... theInitArgs ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    std::vector<_tensor_type> res;
                    res.reserve( t.size() );
                    int k=0;
                    for ( auto const& e : t )
                        res.push_back( _tensor_type( std::true_type{}, exprExpanded[k++], ttse, e, geom, theInitArgs... ) );
                    return res;
                }

        };

        using tuple_tensor_expr_type = std::decay_t<decltype( hana::transform( tuple_optional_expr_type{}, TransformExprToTensor{} ) ) >;

        typedef typename first_expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> first_tensor_expr_type;
        typedef typename first_tensor_expr_type::value_type value_type;
        typedef typename first_tensor_expr_type::shape expr_shape;
        typedef expr_shape shape;

        struct is_zero
        {
            static const bool value = false;//tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tupleTensorExprs( hana::transform( expr.tupleExpressions(), [this,&geom,&fev,&feu](auto const& t) { return TransformExprToTensor{}(t,geom,fev,feu); } ) )
            {}
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tupleTensorExprs( hana::transform( expr.tupleExpressions(), [this,&geom,&fev](auto const& t) { return TransformExprToTensor{}(t,geom,fev); } ) )
            {}
        tensor( this_type const& expr,
                Geo_t const& geom )
            :
            M_tupleTensorExprs( hana::transform( expr.tupleExpressions(), [this,&geom](auto const& t) { return TransformExprToTensor{}(t,geom); } ) )
            {}
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_tupleTensorExprs( hana::transform( hana::unpack( hana::make_range( hana::int_c<0>, hana::int_c<nExpr> ), hana::make_tuple ), [this,&exprExpanded,&ttse,&expr,&geom,&theInitArgs...](auto eId )
                                                 {
                                                     return TransformExprToTensor{}(std::true_type{},hana::at( exprExpanded.tupleExpressions(), hana::int_c<eId> ), ttse,
                                                                                    hana::at( expr.tupleExpressions(), hana::int_c<eId> ), geom, theInitArgs...);
                                                 } ) )
            {}
        template<typename IM>
        void init( IM const& im )
        {
            //M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            hana::for_each( M_tupleTensorExprs, [&geom,&fev,&feu]( auto & e )
                            {
                                for ( auto & e2 : e )
                                    e2.update( geom, fev, feu );
                            });
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            hana::for_each( M_tupleTensorExprs, [&geom,&fev]( auto & e )
                            {
                                for ( auto & e2 : e )
                                    e2.update( geom, fev );
                            });
        }
        void update( Geo_t const& geom )
        {
            hana::for_each( M_tupleTensorExprs, [&geom]( auto & e )
                            {
                                for ( auto & e2 : e )
                                    e2.update( geom );
                            });
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            hana::for_each( M_tupleTensorExprs, [&geom,&face]( auto & e )
                            {
                                for ( auto & e2 : e )
                                    e2.update( geom, face );
                            });
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nExpr> ), [this,&exprExpanded,&ttse,&geom,&theUpdateArgs...]( auto eId )
                            {
                                    auto const& subExprExpanded = hana::at( exprExpanded.tupleExpressions(), hana::int_c<eId> );
                                    int k=0;
                                    for ( auto & e2 :  hana::at( M_tupleTensorExprs, hana::int_c<eId> ) )
                                        e2.update( std::true_type{}, subExprExpanded[k++], ttse, geom, theUpdateArgs... );
                            } );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res(0);
            hana::for_each( M_tupleTensorExprs, [&i,&j,&c1,&c2,&q,&res]( auto const& e )
                            {
                                for ( auto & e2 : e )
                                    res += e2.evalijq( i, j, c1, c2, q );
                            });
            return res;
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res(0);
            hana::for_each( M_tupleTensorExprs, [&i,&c1,&c2,&q,&res]( auto const& e )
                            {
                                for ( auto & e2 : e )
                                    res += e2.evaliq( i, c1, c2, q );
                            });
            return res;
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res(0);
            hana::for_each( M_tupleTensorExprs, [&c1,&c2,&q,&res]( auto const& e )
                            {
                                for ( auto & e2 : e )
                                    res += e2.evalq( c1, c2, q );
                            });
            return res;
        }

    private :
        tuple_tensor_expr_type M_tupleTensorExprs;
    };

private :
    tuple_optional_expr_type M_exprs;
};


template <typename Expr1Type, typename ... ExprOtherType>
inline
Expr<ExprOptionalConcat<hana::tuple<Expr1Type,ExprOtherType...>>>
exprOptionalConcat()
{
    typedef ExprOptionalConcat<hana::tuple<Expr1Type,ExprOtherType...>> eoc_t;
    return Expr< eoc_t >( eoc_t{} );
}

} // namespace vf
} // namespace Feel

#endif
