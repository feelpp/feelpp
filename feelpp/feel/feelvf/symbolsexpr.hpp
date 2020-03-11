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

struct SymbolExprComponentSuffix : public std::vector< std::tuple<std::string,std::array<uint16_type,2>>>
{
    SymbolExprComponentSuffix() = default;
    SymbolExprComponentSuffix( SymbolExprComponentSuffix const& ) = default;
    SymbolExprComponentSuffix( SymbolExprComponentSuffix && ) = default;

    SymbolExprComponentSuffix( uint16_type nComp1, uint16_type nComp2, bool useXYZ = false )
    {
        if ( nComp1 > 1 && nComp2 > 1 )
        {
            for (uint16_type c1=0;c1<nComp1;++c1 )
                for (uint16_type c2=0;c2<nComp2;++c2 )
                    this->push_back( std::make_tuple( this->generateSuffix( c1,c2,useXYZ ), std::array<uint16_type,2>{ { c1,c2 } } ) );
        }
        else if ( nComp2 == 1 && nComp1 > 1 )
        {
            for (uint16_type c1=0;c1<nComp1;++c1 )
                this->push_back( std::make_tuple( this->generateSuffix( c1,useXYZ ), std::array<uint16_type,2>{ { c1,0 } } ) );
        }
        else if ( nComp1 == 1 && nComp2 > 1 )
        {
            for (uint16_type c2=0;c2<nComp2;++c2 )
                this->push_back( std::make_tuple( this->generateSuffix( c2,useXYZ ), std::array<uint16_type,2>{ { 0,c2 } } ) );
        }
        if ( !this->empty() )
            this->shrink_to_fit();
    }

    void print() const
    {
        for ( auto const& d : *this )
            std::cout << "suffix " << std::get<0>( d ) << " c1=" << std::get<1>( d )[0] << " c1=" << std::get<1>( d )[1] << std::endl;
    }
  private :
    static std::string generateSuffix( uint16_type c1, bool useXYZ )
    {
        return "_" + convertComponent( c1, useXYZ );
    }
    static std::string generateSuffix( uint16_type c1, uint16_type c2, bool useXYZ )
    {
        return "_" + convertComponent( c1, useXYZ ) + convertComponent( c2, useXYZ );
    }
   static std::string convertComponent( uint16_type c1, bool useXYZ )
    {
        if ( useXYZ )
        {
            switch ( c1 )
            {
            case 0 : return "x";
            case 1 : return "y";
            case 2 : return "z";
            default : CHECK( false ) << "can't use xyz suffix if comp > 2";
            }
            return "";
        }
        else
        {
            return (boost::format("%1%")%c1).str();
        }
    }
};


struct SymbolExprTag {};

//! attach a symbol (string) with a feel++ expression
//! ex : auto se = SymbolExpr( "u", cst(3.)*idv(u) );
template <typename ExprT>
struct SymbolExpr : public std::vector<std::tuple<std::string,ExprT,SymbolExprComponentSuffix>>
{
    using super_type = std::vector<std::tuple<std::string,ExprT,SymbolExprComponentSuffix>>;
    using feelpp_tag = SymbolExprTag;
    SymbolExpr() = default;
    SymbolExpr( SymbolExpr const& ) = default;
    SymbolExpr( SymbolExpr && ) = default;

    explicit SymbolExpr( std::tuple<std::string,ExprT,SymbolExprComponentSuffix> const& e ) : super_type( 1,e ) {}
    explicit SymbolExpr( std::tuple<std::string,ExprT,SymbolExprComponentSuffix> && e ) : super_type( 1,e ) {}
    explicit SymbolExpr( std::initializer_list<std::tuple<std::string,ExprT,SymbolExprComponentSuffix>> const& e ) : super_type( e ) {}

    explicit SymbolExpr( std::vector<std::tuple<std::string,ExprT>> const& e )
        :
        super_type()
    {
        this->reserve( e.size() );
        SymbolExprComponentSuffix emptySuffix;
        for ( int k=0;k<e.size();++k )
            this->push_back( std::make_tuple( std::get<0>( e[k] ), std::get<1>( e[k] ), emptySuffix ) );
    }
    explicit SymbolExpr( std::vector<std::pair<std::string,ExprT>> const& e )
        :
        super_type()
    {
        this->reserve( e.size() );
        SymbolExprComponentSuffix emptySuffix;
        for ( int k=0;k<e.size();++k )
            this->push_back( std::make_tuple( e[k].first, e[k].second, emptySuffix ) );
    }
    SymbolExpr( super_type const& e ) : super_type( e ) {}
};
//! build a SymbolExpr object
template <typename T>
SymbolExpr<Expr<T>>
symbolExpr( std::string const& s,Expr<T> const& e, SymbolExprComponentSuffix secs = SymbolExprComponentSuffix() ) { return SymbolExpr<Expr<T>>( std::make_tuple(s,e,secs) ); }

template <typename T>
SymbolExpr<Expr<T>>
symbolExpr( std::initializer_list<std::pair<std::string,Expr<T>>> const& e ) { return SymbolExpr<Expr<T>>( std::vector<std::pair<std::string,Expr<T>>>( e ) ); }

template <typename T>
SymbolExpr<Expr<T>>
symbolExpr( std::initializer_list<std::tuple<std::string,Expr<T>>> const& e ) { return SymbolExpr<Expr<T>>( std::vector<std::tuple<std::string,Expr<T>>>( e ) ); }

template <typename T>
SymbolExpr<Expr<T>>
symbolExpr( std::initializer_list<std::tuple<std::string,Expr<T>,SymbolExprComponentSuffix>> const& e ) { return SymbolExpr<Expr<T>>( e ); }

template <typename T>
SymbolExpr<Expr<T>>
symbolExpr( std::vector<std::pair<std::string,Expr<T>>> const& e ) { return SymbolExpr<Expr<T>>( e ); }

template <typename T>
SymbolExpr<Expr<T>>
symbolExpr( std::vector<std::tuple<std::string,Expr<T>>> const& e ) { return SymbolExpr<Expr<T>>( e ); }

template <typename T>
SymbolExpr<Expr<T>>
symbolExpr( std::vector<std::tuple<std::string,Expr<T>,SymbolExprComponentSuffix>> const& e ) { return SymbolExpr<Expr<T>>( e ); }


struct SymbolsExprTag {};

//! defined type from input args (variadic expression)
struct SymbolsExprTraits
{
    template <typename Tag, typename T>
    struct is_a_t :
        hana::integral_constant<bool, std::is_same<Tag, typename T::feelpp_tag >::value >
    {};

    template<typename... ExprT>
    static constexpr auto apply( const ExprT&... exprs )
        {
            return applyImpl( hana::tuple<>{}, exprs... );
        }
private :

    template<typename ResType >
    static constexpr auto applyImpl( ResType && res )
        {
            return std::move( res );
        }

    template<typename ResType >
    static constexpr auto applyImpl2( ResType && res,  hana::tuple<> const& t )
        {
            return std::move( res );
        }
    template<typename ResType, typename T1, typename... ExprT >
    static constexpr auto applyImpl2( ResType && res,  hana::tuple<T1,ExprT...> const& t )
        {
            return applyImpl2( applyImpl( std::forward<ResType>( res ), hana::at( t, 0_c ) ), hana::remove_at( t, 0_c ) );
        }

    template < typename T1, typename... ExprT >
    static constexpr auto applySymbolExpr( T1 const& t1, hana::tuple<ExprT...> && res )
        {
            if constexpr ( hana::find( hana::to_tuple(hana::tuple_t<ExprT...> ), hana::type_c<T1>) == hana::nothing )
                         {
                             return hana::append( res, t1 );
                         }
            else
            {
                hana::for_each( res, [&t1]( auto & e )
                                {
                                    if constexpr ( std::is_same_v<std::decay_t<decltype(e)>, T1> )
                                        {
                                            for ( auto const& se : t1 )
                                                e.push_back( se );
                                        }
                                });
                return std::move( res );
            }
        }
    template<typename ResType, typename T1, typename... ExprT2 >
    static constexpr auto applyImpl( ResType && res, T1 const& t1, const ExprT2&... exprs )
        {
            if constexpr ( is_a_t<SymbolExprTag, T1 >::value )
                {
                    return applyImpl( applySymbolExpr(t1,std::forward<ResType>(res) ), exprs... );
                }
            else if constexpr ( is_a_t<SymbolsExprTag, T1 >::value )
                {
                    if constexpr ( std::decay_t<decltype(hana::size(t1.tupleExpr))>::value == 0 )
                                     return applyImpl( std::forward<ResType>( res ), exprs... );
                        else
                            return applyImpl( applyImpl2( std::forward<ResType>( res ), t1.tupleExpr ), exprs... );
                }
            else
                return std::move( res );
        }

};


//! store set of SymbolExpr object into a hana::tuple
template<typename TupleExprType>
struct SymbolsExpr
{
    using tuple_type = TupleExprType;
    using feelpp_tag = SymbolsExprTag;

    SymbolsExpr() = default;

    explicit SymbolsExpr( tuple_type const& tse )
        :
        tupleExpr( tse )
        {}

    explicit SymbolsExpr( tuple_type && tse )
        :
        tupleExpr( tse )
        {}

    std::map<std::string,std::set<std::string>> names() const
        {
            std::map<std::string,std::set<std::string>> res;
            hana::for_each( tupleExpr, [&res]( auto const& e )
                            {
                                for ( auto const& se : e )
                                {
                                    std::string const& symbolNameBase = std::get<0>( se );
                                    SymbolExprComponentSuffix const& symbolSuffix = std::get<2>( se );
                                    if ( symbolSuffix.empty() )
                                        res[symbolNameBase].insert( symbolNameBase );
                                    else
                                    {
                                        for ( auto const& [_suffix,compArray] : symbolSuffix )
                                            res[symbolNameBase].insert( symbolNameBase+ _suffix );
                                    }
                                }
                            });
            return res;
        }

    tuple_type tupleExpr;
};
template <>
struct SymbolsExpr<hana::tuple<>>
{
    using tuple_type = hana::tuple<>;
    using feelpp_tag = SymbolsExprTag;
    tuple_type tupleExpr;

    SymbolsExpr() = default;
    explicit SymbolsExpr( tuple_type const& tse ) {}
    explicit SymbolsExpr( tuple_type && tse ) {}

    std::map<std::string,std::set<std::string>> names() const { return std::map<std::string,std::set<std::string>>{}; }

};

template<typename... ExprT>
struct SymbolsExprTraits2
{
    static constexpr auto callApply = [](const auto& ...exprs) { return SymbolsExprTraits::template apply( exprs... ); };
    using tuple_type = std::decay_t<decltype( hana::unpack( hana::tuple<ExprT...>{},  callApply ) )>;
    using type = SymbolsExpr<tuple_type>;
};

template<typename... ExprT>
using symbols_expression_t = typename SymbolsExprTraits2<ExprT...>::type;

using symbols_expression_empty_t = SymbolsExpr<hana::tuple<>>;

//! build a SymbolsExpr object
template<typename... ExprT>
symbols_expression_t<ExprT...>
symbolsExpr( const ExprT&... exprs )
{
    return symbols_expression_t<ExprT...>(SymbolsExprTraits::template apply( exprs... ) );
}

} // namespace vf

} // namespace Feel

#endif
