/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2018-06-03

  Copyright (C) 2018 Feel++ Consortium

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
#ifndef FEELPP_VF_SYMBOLSEXPR_HPP
#define FEELPP_VF_SYMBOLSEXPR_HPP 1

namespace Feel
{
namespace vf
{

struct SymbolExprTag {};

//! attach a symbol (string) with a feel++ expression
//! ex : auto se = SymbolExpr( "u", cst(3.)*idv(u) );
template <typename ExprT>
struct SymbolExpr : public std::vector<std::pair<std::string,ExprT>>
{
    using super_type = std::vector<std::pair<std::string,ExprT>>;
    using feelpp_tag = SymbolExprTag;
    SymbolExpr() = default;
    SymbolExpr( SymbolExpr const& ) = default;
    SymbolExpr( SymbolExpr && ) = default;
    SymbolExpr( std::pair<std::string,ExprT> const& e ) : super_type( 1,e ) {}
    SymbolExpr( std::pair<std::string,ExprT> && e ) : super_type( 1,e ) {}
    SymbolExpr( std::initializer_list<std::pair<std::string,ExprT>> const& e ) : super_type( e ) {}
};
//! build a SymbolExpr object
template <typename T>
SymbolExpr<Expr<T>>
symbolExpr( std::string const& s,Expr<T> const& e ) { return SymbolExpr<Expr<T>>( std::make_pair(s,e) ); }

template <typename T>
SymbolExpr<Expr<T>>
symbolExpr( std::initializer_list<std::pair<std::string,Expr<T>>> const& e ) { return SymbolExpr<Expr<T>>( e ); }


struct SymbolsExprTag {};

//! defined type from input args (variadic expression)
struct SymbolsExprTraits
{
    template <typename Tag, typename T>
    struct is_a_t :
        hana::integral_constant<bool, std::is_same<Tag, typename T::feelpp_tag >::value >
    {};

    template<typename... ExprT2>
    static constexpr auto apply2( hana::integral_constant<bool,false> )
        {
            return hana::type_c<hana::tuple<>>;
        }
    template<typename T1, typename... ExprT2>
    static constexpr auto apply2NonEmpty( hana::integral_constant<bool,true> )
        {
            return hana::type_c<hana::tuple<T1,ExprT2...>>;
        }
    template<typename T1, typename... ExprT2>
    static constexpr auto apply2NonEmpty( hana::integral_constant<bool,false> )
        {
            return hana::type_c<T1>;
        }
    template<typename T1, typename... ExprT2>
    static constexpr auto apply2( hana::integral_constant<bool,true> )
        {
            return apply2NonEmpty<T1,ExprT2...>( is_a_t<SymbolExprTag, T1 >() );
        }
    template<typename... ExprT>
    static constexpr auto apply2()
        {
            constexpr int nExpr = decltype(hana::length( hana::tuple<ExprT...>{} ) )::value;
            constexpr bool hasExpr = nExpr > 0;
            return apply2<ExprT...>( hana::integral_constant<bool, hasExpr >() );
        }


    template<typename T1, typename... ExprT2>
    static constexpr auto applyNonEmpty( hana::integral_constant<bool,true>, T1 const& t1, const ExprT2&... exprs )
        {
            return hana::tuple<T1,ExprT2...>( t1,exprs... );
        }
    template<typename T1, typename... ExprT2>
    static constexpr auto applyNonEmpty( hana::integral_constant<bool,false>, T1 const& t1, const ExprT2&... exprs )
        {
            return t1;
        }
    template<typename... ExprT>
    static constexpr auto apply( hana::integral_constant<bool,false> )
        {
            return hana::tuple<>{};
        }
    template<typename T1, typename... ExprT2 >
    static constexpr auto apply( hana::integral_constant<bool,true>, T1 const& t1, const ExprT2&... exprs )
        {
            return applyNonEmpty( is_a_t<SymbolExprTag, T1 >(), t1, exprs... );
        }
    template<typename... ExprT>
    static constexpr auto apply( const ExprT&... exprs )
        {
            constexpr int nExpr = decltype(hana::length( hana::tuple<ExprT...>{} ) )::value;
            constexpr bool hasExpr = nExpr > 0;
            return apply( hana::integral_constant<bool, hasExpr >(), exprs... );
        }
};

//! collect SymbolExpr object and store it an hana::tuple
template<typename... ExprT>
struct SymbolsExpr : public decltype(SymbolsExprTraits::template apply2<ExprT...>() )::type
{
    using super_type = typename decltype(SymbolsExprTraits::template apply2<ExprT...>() )::type;
    using feelpp_tag = SymbolsExprTag;

    //SymbolsExpr() = default;
    template<typename ...dummy,typename = typename std::enable_if< decltype(hana::length( hana::tuple<ExprT...,dummy...>{} ) )::value != 0 >::type >
    SymbolsExpr()
        :
        super_type()
    {}
    SymbolsExpr( const ExprT&... exprs )
        :
        super_type( SymbolsExprTraits::template apply( exprs... ) )
        {}
};
template <>
struct SymbolsExpr<> : public hana::tuple<>
{
    hana::basic_tuple<> storage_; // fix empty hana::tuple used hana::fold
};

//! build a SymbolsExpr object
template<typename... ExprT>
SymbolsExpr<ExprT...>
symbolsExpr( const ExprT&... exprs ) { return SymbolsExpr<ExprT...>( exprs... ); }

}

} // namespace Feel

#endif
