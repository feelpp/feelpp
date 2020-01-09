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


    template<typename T1/*, typename = typename std::enable_if< is_a_t<SymbolExprTag, T1 >::value >::type*/>
    static constexpr auto applyOneElement( T1 const& t1, std::enable_if_t< is_a_t<SymbolExprTag, T1 >::value >* = nullptr )
        {
            return hana::tuple<T1>( t1 );
        }
    template<typename T1/*, typename = typename std::enable_if< is_a_t<SymbolsExprTag, T1 >::value >::type*/>
    static constexpr auto applyOneElement( T1 const& t1, std::enable_if_t< is_a_t<SymbolsExprTag, T1 >::value >* = nullptr )
        {
            return t1.tupleExpr;
        }

    template<typename T1, typename... ExprT2>
    static constexpr auto applyNonEmpty( T1 const& t1, const ExprT2&... exprs, std::enable_if_t< decltype(hana::length( hana::tuple<ExprT2...>{} ) )::value != 0 >* = nullptr )
        {
            return hana::insert_range( applyNonEmpty<ExprT2...>( exprs... ), 0_c,  applyOneElement<T1>( t1 ) );
        }
    template<typename T1>
    static constexpr auto applyNonEmpty( T1 const& t1 )
        {
            return applyOneElement<T1>( t1 );
        }

    template<typename... ExprT>
    static constexpr auto apply( hana::integral_constant<bool,false> )
        {
            return hana::tuple<>{};
        }
    template<typename T1, typename... ExprT2 >
    static constexpr auto apply( hana::integral_constant<bool,true>, T1 const& t1, const ExprT2&... exprs )
        {
            return applyNonEmpty<T1,ExprT2...>( t1, exprs... );
        }
    template<typename... ExprT>
    static constexpr auto apply( const ExprT&... exprs )
        {
            constexpr int nExpr = decltype(hana::length( hana::tuple<ExprT...>{} ) )::value;
            constexpr bool hasExpr = nExpr > 0;
            return apply( hana::integral_constant<bool, hasExpr >(), exprs... );
        }

#if 1 // maybe the code below can be remove in the future, see comments below
    template<typename... ExprT2>
    static constexpr auto apply2( hana::integral_constant<bool,false> )
        {
            return hana::type_c<hana::tuple<>>;
        }

    template<typename T1/*, typename = typename std::enable_if< is_a_t<SymbolExprTag, T1 >::value >::type*/>
    static constexpr auto apply2OneElement( std::enable_if_t< is_a_t<SymbolExprTag, T1 >::value >* = nullptr )
        {
            return hana::type_c<hana::tuple<T1>>;
        }
    template<typename T1/*, typename = typename std::enable_if< is_a_t<SymbolsExprTag, T1 >::value >::type*/>
    static constexpr auto apply2OneElement( std::enable_if_t< is_a_t<SymbolsExprTag, T1 >::value >* = nullptr )
        {
            return hana::type_c<typename T1::tuple_type>;
        }

    template<typename T1, typename... ExprT2>
    static constexpr auto apply2NonEmpty( std::enable_if_t< decltype(hana::length( hana::tuple<ExprT2...>{} ) )::value != 0 >* = nullptr )
        {
            return hana::type_c< decltype(hana::insert_range( typename decltype(apply2NonEmpty<ExprT2...>())::type{}  ,0_c, typename  decltype(apply2OneElement<T1>())::type{} ) ) >;
        }
    template<typename T1>
    static constexpr auto apply2NonEmpty()
        {
            return apply2OneElement<T1>();
        }
    template<typename T1, typename... ExprT2>
    static constexpr auto apply2( hana::integral_constant<bool,true> )
        {
            return apply2NonEmpty<T1,ExprT2...>();
        }
    template<typename... ExprT>
    static constexpr auto apply2()
        {
            constexpr int nExpr = decltype(hana::length( hana::tuple<ExprT...>{} ) )::value;
            constexpr bool hasExpr = nExpr > 0;
            return apply2<ExprT...>( hana::integral_constant<bool, hasExpr >() );
        }
#endif


};

//! collect SymbolExpr object and store it an hana::tuple
template<typename... ExprT>
struct SymbolsExpr
{
#if 0
    // NOTE : a lot of code work with this 2 lines below but some codes (as test_modelproperties) does not compile.
    // The error is : constexpr variable cannot have non-literal type ....
    // We have also this message at the end : lambda closure types are non-literal types before C++17
    // So, maybe it's ok C++17
    static constexpr auto callApply = [](const auto& ...exprs) { return SymbolsExprTraits::template apply( exprs... ); };
    using tuple_type = decltype( hana::unpack( hana::tuple<ExprT...>{},  callApply ) );
#else
    using tuple_type = typename decltype(SymbolsExprTraits::template apply2<ExprT...>() )::type;
#endif

    using feelpp_tag = SymbolsExprTag;

    //SymbolsExpr() = default;
    template<typename ...dummy,typename = typename std::enable_if< decltype(hana::length( hana::tuple<ExprT...,dummy...>{} ) )::value != 0 >::type >
    SymbolsExpr()
        :
        tupleExpr()
    {}
    SymbolsExpr( const ExprT&... exprs )
        :
        tupleExpr( SymbolsExprTraits::template apply( exprs... ) )
        {}

    tuple_type tupleExpr;
};
template <>
struct SymbolsExpr<>
{
    using tuple_type = hana::tuple<>;
    tuple_type tupleExpr;
};

//! build a SymbolsExpr object
template<typename... ExprT>
SymbolsExpr<ExprT...>
symbolsExpr( const ExprT&... exprs ) { return SymbolsExpr<ExprT...>( exprs... ); }

} // namespace vf

} // namespace Feel

#endif
