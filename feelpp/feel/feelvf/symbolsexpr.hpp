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

#include <feel/feelcore/tuple_utils.hpp>

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

using SymbolExprUpdateFunction = std::function<void()>;

struct SymbolExprTag {};

template <typename ExprT>
struct SymbolExpr1
{
    using expr_type = ExprT;
    SymbolExpr1( std::string const& s, ExprT const& e, SymbolExprComponentSuffix const& secs = SymbolExprComponentSuffix(), SymbolExprUpdateFunction const& seuf = SymbolExprUpdateFunction{} )
        :
        M_symbol( s ),
        M_expr( e ),
        M_secs( secs ),
        M_seuf( seuf )
        {}
    SymbolExpr1( SymbolExpr1 const& ) = default;
    SymbolExpr1( SymbolExpr1 && ) = default;

    std::string const& symbol() const { return M_symbol; }
    expr_type const& expr() const { return M_expr; }
    expr_type & expr() { return M_expr; }
    SymbolExprComponentSuffix const& componentSuffix() const { return M_secs; }
    SymbolExprUpdateFunction const& updateFunction() const { return M_seuf; }

private :
    std::string M_symbol;
    expr_type M_expr;
    SymbolExprComponentSuffix M_secs;
    SymbolExprUpdateFunction M_seuf;
};

//! attach a symbol (string) with a feel++ expression
//! ex : auto se = SymbolExpr( "u", cst(3.)*idv(u) );
template <typename ExprT>
struct SymbolExpr : public std::vector<SymbolExpr1<ExprT>>
{
    using super_type = std::vector<SymbolExpr1<ExprT>>;
    using symbolexpr1_type = SymbolExpr1<ExprT>;
    using update_function_type = SymbolExprUpdateFunction;
    using feelpp_tag = SymbolExprTag;
    SymbolExpr() = default;
    SymbolExpr( SymbolExpr const& ) = default;
    SymbolExpr( SymbolExpr && ) = default;
    explicit SymbolExpr( typename super_type::value_type const& e ) : super_type( 1,e ) {}
    explicit SymbolExpr( typename super_type::value_type && e ) : super_type( 1,e ) {}
    explicit SymbolExpr( std::initializer_list<typename super_type::value_type> const& e ) : super_type( e ) {}
    explicit SymbolExpr( std::vector<std::tuple<std::string,ExprT>> const& e )
        :
        super_type()
    {
        this->reserve( e.size() );
        SymbolExprComponentSuffix emptySuffix;
        update_function_type emptyUpdateFunc;
        for ( int k=0;k<e.size();++k )
            this->push_back( symbolexpr1_type( std::get<0>( e[k] ), std::get<1>( e[k] ), emptySuffix, emptyUpdateFunc ) );
    }
    explicit SymbolExpr( std::vector<std::pair<std::string,ExprT>> const& e )
        :
        super_type()
    {
        this->reserve( e.size() );
        SymbolExprComponentSuffix emptySuffix;
        update_function_type emptyUpdateFunc;
        for ( int k=0;k<e.size();++k )
            this->push_back( symbolexpr1_type( e[k].first, e[k].second, emptySuffix, emptyUpdateFunc ) );
    }
    SymbolExpr( super_type const& e ) : super_type( e ) {}


    void add( std::string const& s, ExprT const& e, SymbolExprComponentSuffix const& secs = SymbolExprComponentSuffix(), SymbolExprUpdateFunction const& seuf = SymbolExprUpdateFunction{} )
    {
        this->push_back( symbolexpr1_type( s,e,secs,seuf ) );
    }
};

template<typename ExprT>
using symbol_expression_t = SymbolExpr<ExprT>;

//! build a SymbolExpr object
template <typename ExprT>
SymbolExpr<ExprT>
symbolExpr( std::string const& s,ExprT const& e, SymbolExprComponentSuffix const& secs = SymbolExprComponentSuffix(), SymbolExprUpdateFunction const& seuf = SymbolExprUpdateFunction{} )
{
    return SymbolExpr<ExprT>( typename SymbolExpr<ExprT>::symbolexpr1_type(s,e,secs,seuf) );
}

template <typename ExprT>
SymbolExpr<ExprT>
symbolExpr( std::initializer_list<std::pair<std::string,ExprT>> const& e ) { return SymbolExpr<ExprT>( std::vector<std::pair<std::string,ExprT>>( e ) ); }

template <typename ExprT>
SymbolExpr<ExprT>
symbolExpr( std::initializer_list<std::tuple<std::string,ExprT>> const& e ) { return SymbolExpr<ExprT>( std::vector<std::tuple<std::string,ExprT>>( e ) ); }

template <typename ExprT>
SymbolExpr<ExprT>
symbolExpr( std::initializer_list<std::tuple<std::string,ExprT,SymbolExprComponentSuffix,SymbolExprUpdateFunction>> const& e ) { return SymbolExpr<ExprT>( e ); }

template <typename ExprT>
SymbolExpr<ExprT>
symbolExpr( std::vector<std::pair<std::string,ExprT>> const& e ) { return SymbolExpr<ExprT>( e ); }

template <typename ExprT>
SymbolExpr<ExprT>
symbolExpr( std::vector<std::tuple<std::string,ExprT>> const& e ) { return SymbolExpr<ExprT>( e ); }

template <typename ExprT>
SymbolExpr<ExprT>
symbolExpr( std::vector<std::tuple<std::string,ExprT,SymbolExprComponentSuffix,SymbolExprUpdateFunction>> const& e ) { return SymbolExpr<ExprT>( e ); }


struct SymbolsExprTag {};

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
                                    std::string const& symbolNameBase = se.symbol();
                                    SymbolExprComponentSuffix const& symbolSuffix = se.componentSuffix();
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

    tuple_type const& tuple() const { return tupleExpr; }
    tuple_type & tuple() { return tupleExpr; }

    tuple_type tupleExpr;
};
#if 0
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
#endif
template<typename... ExprT>
struct SymbolsExprTraits
{
    static constexpr auto callApply = [](const auto& ...exprs) { return Feel::detail::AdvancedConcatOfTupleContainerType<SymbolsExprTag,SymbolExprTag>::template apply( exprs... ); };
    using tuple_type = std::decay_t<decltype( hana::unpack( hana::tuple<ExprT...>{},  callApply ) )>;
    using type = SymbolsExpr<tuple_type>;
};

template<typename... ExprT>
using symbols_expression_t = typename SymbolsExprTraits<ExprT...>::type;

using symbols_expression_empty_t = SymbolsExpr<hana::tuple<>>;

//! build a SymbolsExpr object
template<typename... ExprT>
symbols_expression_t<ExprT...>
symbolsExpr( const ExprT&... exprs )
{
    return symbols_expression_t<ExprT...>(Feel::detail::AdvancedConcatOfTupleContainerType<SymbolsExprTag,SymbolExprTag>::template apply( exprs... ) );
}

} // namespace vf

} // namespace Feel

#endif
