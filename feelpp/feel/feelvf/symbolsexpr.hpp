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

#include <any>
#include <feel/feelcore/tuple_utils.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelvf/detail/gmc.hpp>


namespace Feel
{
namespace vf
{

template<typename Geo_t/*, typename Basis_i_t, typename Basis_j_t*/>
struct ExprTensorsFromSymbolsExpr;


template<typename Geo_t/*, typename Basis_i_t, typename Basis_j_t*/>
struct tensorGen
{
    using value_type = double;

    virtual ~tensorGen() {}

    virtual void update( std::true_type /**/, ExprTensorsFromSymbolsExpr<Geo_t/*,Basis_i_t,Basis_j_t*/> & ttse, Geo_t const& geom ) = 0;
    virtual value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const = 0;

};


template<typename Geo_t/*, typename Basis_i_t, typename Basis_j_t*/, typename ExprType, typename ExprExpandType>
struct tensorFromExpr : public tensorGen<Geo_t/*,Basis_i_t,Basis_j_t*/>
{
    using super_type = tensorGen<Geo_t/*,Basis_i_t,Basis_j_t*/>;
    using value_type = typename super_type::value_type;
    using expr_type = ExprType;
    using expr_expand_type = ExprExpandType;
    using tensor_type = typename ExprType::template tensor<Geo_t/*,Basis_i_t,Basis_j_t*/>;

    tensorFromExpr( std::true_type /**/, expr_expand_type const& exprExpanded, ExprTensorsFromSymbolsExpr<Geo_t/*,Basis_i_t,Basis_j_t*/> & ttse,
                    expr_type const& expr, Geo_t const& geom/*, const TheArgsType&... theInitArgs*/ )
        :
        M_expr( std::make_shared<expr_type>( expr ) ),
        M_exprExpand( exprExpanded ),
        M_tensor( std::true_type{}, exprExpanded, ttse, *M_expr, geom/*, theInitArgs...*/ )
        {}

    ~tensorFromExpr() override {}

    void update( std::true_type /**/, ExprTensorsFromSymbolsExpr<Geo_t/*,Basis_i_t,Basis_j_t*/> & ttse, Geo_t const& geom ) override
        {
            M_tensor.update( std::true_type{}, M_exprExpand, ttse, geom );
        }

    value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const override { return M_tensor.evalq(c1,c2,q); }

private :
    std::shared_ptr<expr_type> M_expr;
    expr_expand_type const& M_exprExpand;
    tensor_type M_tensor;
};

template<typename Geo_t/*, typename Basis_i_t, typename Basis_j_t*/, typename ExprType>
struct tensorFromExprClassic : public tensorGen<Geo_t/*,Basis_i_t,Basis_j_t*/>
{
    using super_type = tensorGen<Geo_t/*,Basis_i_t,Basis_j_t*/>;
    using value_type = typename super_type::value_type;
    using expr_type = ExprType;
    using tensor_type = typename expr_type::template tensor<Geo_t/*,Basis_i_t,Basis_j_t*/>;

    tensorFromExprClassic( expr_type const& expr, Geo_t const& geom/*, const TheArgsType&... theInitArgs*/ )
        :
        M_tensor( expr, geom/*, theInitArgs...*/ )
        {}

    ~tensorFromExprClassic() override {}

    void update( std::true_type /**/, ExprTensorsFromSymbolsExpr<Geo_t/*,Basis_i_t,Basis_j_t*/> & ttse, Geo_t const& geom ) override
        {
            M_tensor.update( geom );
        }

    value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const override { return M_tensor.evalq(c1,c2,q); }

private :
    tensor_type M_tensor;
};


template<typename Geo_t/*, typename Basis_i_t, typename Basis_j_t*/>
struct ExprTensorsFromSymbolsExpr
{
    using self_type = ExprTensorsFromSymbolsExpr<Geo_t/*,Basis_i_t,Basis_j_t*/>;
    using tensor_type = tensorGen<Geo_t/*,Basis_i_t,Basis_j_t*/>;
    using tensor_ptrtype = std::shared_ptr<tensor_type>;
    using value_type = typename tensor_type::value_type;

    virtual ~ExprTensorsFromSymbolsExpr() {}

    virtual std::shared_ptr<self_type> newObjectPtr() const = 0;

    virtual std::map<uint16_type,uint16_type> init( std::vector<std::any> const& evecExpand,
                                                     std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type>>> const& M_t_expr_index,
                                                     Geo_t const& geom ) = 0;

    virtual void update( std::vector<std::any> const& evecExpand,
                         std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type>>> const& M_t_expr_index,
                         std::map<uint16_type,uint16_type> const& M_subexprIdToSubtensorId,
                         std::vector<std::vector<std::pair<uint16_type,value_type>>> & evalSubtensor,
                         Geo_t const& geom ) = 0;


    tensor_type & tensor( uint16_type id )
        {
            CHECK( id < M_tensor.size() ) << "invalid id";
            CHECK( M_tensor[id] ) << "not init";
            return *M_tensor[id];
        }

protected:
    std::vector<tensor_ptrtype> M_tensor;

};

template<typename Geo_t/*, typename Basis_i_t, typename Basis_j_t*/, typename SymbolsExprType>
struct ExprTensorsFromSymbolsExprImpl : public ExprTensorsFromSymbolsExpr<Geo_t/*,Basis_i_t,Basis_j_t*/>
{
    using super_type = ExprTensorsFromSymbolsExpr<Geo_t/*,Basis_i_t,Basis_j_t*/>;
    using self_type = ExprTensorsFromSymbolsExprImpl<Geo_t/*,Basis_i_t,Basis_j_t*/,SymbolsExprType>;
    using symbols_expr_type = SymbolsExprType;
    using value_type = typename super_type::value_type;
private :

    template<template<class...>class Target>
    struct DeferToExpandIsSameType
    {
        template<class Ts1/*,class Ts2,class Ts3*/,class Ts4,class Ts5  >
        using execute = Target<Ts1/*,Ts2,Ts3*/,Ts4>;
    };
    template<template<class...>class Target>
    struct DeferToExpandIsNotSameType
    {
        template<class Ts1/*,class Ts2,class Ts3*/,class Ts4,class Ts5 >
        using execute = Target<Ts1/*,Ts2,Ts3*/,Ts4,Ts5>;
    };

    template< class Prog, class... Ts >
    using runRRR = typename Prog::template execute<Ts...>;

public :

    ExprTensorsFromSymbolsExprImpl( symbols_expr_type const& se ) : M_se( se ) {}

    ~ExprTensorsFromSymbolsExprImpl() override {}

    std::shared_ptr<super_type> newObjectPtr() const override
        {
            return std::make_shared<self_type>( M_se );
        }

    std::map<uint16_type,uint16_type>
    init( std::vector<std::any> const& evecExpand,
          std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type>>> const& M_t_expr_index,
          Geo_t const& geom ) override
        {

            std::map<uint16_type,uint16_type> subexprIdToSubtensorId;

            auto & baseObject = *static_cast<super_type*>(this);
            uint16_type k=0;
            hana::for_each( M_se.tuple(), [this,&geom,&evecExpand,&M_t_expr_index,&k,&subexprIdToSubtensorId,&baseObject]( auto const& evec ) {
                    int nSubExpr = evec.size();
                    using the_expr_type = typename std::decay_t<decltype(evec)>::expr_type;
                    using the_expr_expand_type = std::decay_t<decltype(evec.front().expr().applySymbolsExpr( M_se ))>;
                    //using tensor_from_expr_type = tensorFromExpr<Geo_t,Basis_i_t,Basis_j_t,the_expr_type,the_expr_expand_type>;

                    static const bool is_same_expr = std::is_same_v<the_expr_type,the_expr_expand_type>;
                    using choice = typename std::conditional< is_same_expr, DeferToExpandIsSameType<tensorFromExprClassic>, DeferToExpandIsNotSameType<tensorFromExpr>  >::type;
                    using tensor_from_expr_type = runRRR< choice, Geo_t/*,Basis_i_t,Basis_j_t*/,the_expr_type,the_expr_expand_type >;

                    for (int l=0;l<nSubExpr;++l,++k)
                    {
                        if ( M_t_expr_index[k].empty() )
                            continue;

                        auto const& e = evec[l];
                        auto const& theexprBase = e.expr();
                        auto const& theexpr = std::any_cast<the_expr_expand_type const&>(evecExpand[k]);
                        typename super_type::tensor_ptrtype tPtr;
                        if constexpr ( is_same_expr )
                            tPtr = std::make_shared<tensor_from_expr_type>( theexpr, geom );
                        else
                            tPtr = std::make_shared<tensor_from_expr_type>( std::true_type{}, theexpr, baseObject, theexprBase, geom );

                        uint16_type id = this->M_tensor.size();
                        this->M_tensor.push_back( std::move( tPtr ) );
                        subexprIdToSubtensorId[k] = id;
                    }
                });
            return subexprIdToSubtensorId;
        }

    void update( std::vector<std::any> const& evecExpand,
                 std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type>>> const& M_t_expr_index,
                 std::map<uint16_type,uint16_type> const& M_subexprIdToSubtensorId,
                 std::vector<std::vector<std::pair<uint16_type,value_type>>> & evalSubtensor,
                 Geo_t const& geom ) override
        {
            using key_type = key_t<Geo_t>;
            //using key_type = key_t<Geo_t>;
            auto M_gmc = fusion::at_key<key_type>( geom ).get();

            auto & baseObject = *static_cast<super_type*>(this);
            uint16_type k=0;
            hana::for_each( M_se.tuple(), [this,&geom,&evecExpand,&M_t_expr_index,&M_subexprIdToSubtensorId,&M_gmc,&evalSubtensor,&k,&baseObject]( auto const& evec )
                            {
                                int nSubExpr = evec.size();
                                for (int l=0;l<nSubExpr;++l,++k)
                                {
                                    if ( M_t_expr_index[k].empty() )
                                        continue;

                                    auto itFindId = M_subexprIdToSubtensorId.find( k );
                                    CHECK( itFindId != M_subexprIdToSubtensorId.end() ) << "invalid id";
                                    auto & subTensor = this->tensor( itFindId->second );

                                    subTensor.update( std::true_type{}, baseObject, geom/*, theUpdateArgs...*/ );

                                    for ( auto const& [idx,c1,c2] : M_t_expr_index[k] )
                                    {
                                        for(int q = 0; q < M_gmc->nPoints();++q )
                                            evalSubtensor[q].push_back( std::make_pair(idx, subTensor.evalq( c1, c2, q ) ) );
                                    }
                                }
                            });

        }
private :
    symbols_expr_type const& M_se;
};



template<typename Geo_t/*, typename Basis_i_t, typename Basis_j_t*/>
struct ExprTensorsFromSymbolsExprDummyImpl : public ExprTensorsFromSymbolsExpr<Geo_t/*,Basis_i_t,Basis_j_t*/>
{
    using self_type = ExprTensorsFromSymbolsExprDummyImpl<Geo_t/*,Basis_i_t,Basis_j_t*/>;
    using super_type =  ExprTensorsFromSymbolsExpr<Geo_t/*,Basis_i_t,Basis_j_t*/>;
    using value_type = typename super_type::value_type;

    std::shared_ptr<super_type> newObjectPtr() const override
        {
            return std::make_shared<self_type>();
        }

    ~ExprTensorsFromSymbolsExprDummyImpl() override {}

    std::map<uint16_type,uint16_type>
    init( std::vector<std::any> const& evecExpand,
          std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type>>> const& M_t_expr_index,
          Geo_t const& geom ) override
        {
            return std::map<uint16_type,uint16_type>{};
        }

    void update( std::vector<std::any> const& evecExpand,
                 std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type>>> const& M_t_expr_index,
                 std::map<uint16_type,uint16_type> const& M_subexprIdToSubtensorId,
                 std::vector<std::vector<std::pair<uint16_type,value_type>>> & evalSubtensor,
                 Geo_t const& geom ) override
        {}
};

















             

struct SymbolExprComponentSuffix : public std::vector< std::tuple<std::string,std::array<uint16_type,2>>>
{
    SymbolExprComponentSuffix() : M_nComp1(0), M_nComp2(0) {}
    SymbolExprComponentSuffix( SymbolExprComponentSuffix const& ) = default;
    SymbolExprComponentSuffix( SymbolExprComponentSuffix && ) = default;
    SymbolExprComponentSuffix& operator=( SymbolExprComponentSuffix const& ) = default;
    SymbolExprComponentSuffix& operator=( SymbolExprComponentSuffix&& ) = default;

    SymbolExprComponentSuffix( uint16_type nComp1, uint16_type nComp2, bool useXYZ = false )
        :
        M_nComp1( nComp1 ),
        M_nComp2( nComp2 )
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

    uint16_type nComp1() const { return M_nComp1; }
    uint16_type nComp2() const { return M_nComp2; }

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
  private :
   uint16_type M_nComp1, M_nComp2;
};

//! return informations about a symbol in json format
nl::json symbolExprInformations( std::string const& symbol, std::string const& expr, SymbolExprComponentSuffix const& symbolSuffix, std::string const& name = "" );


using SymbolExprUpdateFunction = std::function<void()>;

struct SymbolExprTag {};

template <typename ExprT>
struct SymbolExpr1
{
    using expr_type = ExprT;

    template <typename TheExprType>
    SymbolExpr1( std::string const& s, TheExprType && e, SymbolExprComponentSuffix const& secs = SymbolExprComponentSuffix(), SymbolExprUpdateFunction const& seuf = SymbolExprUpdateFunction{} )
        :
        M_symbol( s ),
        M_expr( std::forward<TheExprType>( e ) ),
        M_secs( secs ),
        M_seuf( seuf )
        {}
    SymbolExpr1( SymbolExpr1 const& ) = default;
    SymbolExpr1( SymbolExpr1 && ) = default;
    SymbolExpr1& operator=( SymbolExpr1 const& ) = default;
    SymbolExpr1& operator=( SymbolExpr1&& ) = default;

    std::string const& symbol() const { return M_symbol; }
    expr_type const& expr() const { return M_expr; }
    expr_type & expr() { return M_expr; }
    SymbolExprComponentSuffix const& componentSuffix() const { return M_secs; }
    SymbolExprUpdateFunction const& updateFunction() const { return M_seuf; }

    template<typename... TheExpr>
    auto
    applyLambda( TheExpr... e ) const
    {
        using new_expr_type = typename ExprT::template Lambda<TheExpr...>::type;
        return SymbolExpr1<new_expr_type>( M_symbol, M_expr(e...), M_secs, M_seuf );
    }
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
    using symbolexpr1_type = typename super_type::value_type;
    using expr_type = ExprT;

    using update_function_type = SymbolExprUpdateFunction;
    using feelpp_tag = SymbolExprTag;
    SymbolExpr() = default;
    SymbolExpr( SymbolExpr const& ) = default;
    SymbolExpr( SymbolExpr && ) = default;
    SymbolExpr& operator=( SymbolExpr const& ) = default;
    SymbolExpr& operator=( SymbolExpr&& ) = default;
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


    template <typename TheExprType>
    void add( std::string const& s, TheExprType && e, SymbolExprComponentSuffix const& secs = SymbolExprComponentSuffix(), SymbolExprUpdateFunction const& seuf = SymbolExprUpdateFunction{} )
    {
        this->push_back( symbolexpr1_type( s,std::forward<TheExprType>( e ),secs,seuf ) );
    }

    template<typename... TheExpr>
    struct Lambda
    {
        using type = SymbolExpr<typename ExprT::template Lambda<TheExpr...>::type>;
    };

    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    applyLambda( TheExpr... e ) const
    {
        typename Lambda<TheExpr...>::type res;
        res.reserve( this->size() );
        for ( auto const& se1 : *this )
            res.push_back( se1.applyLambda( e... ) );
        return res;
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

struct SymbolsExprBase {};
struct SymbolsExprTensorContextBase {};

struct SymbolsExprTag {};

//! store set of SymbolExpr object into a hana::tuple
template<typename TupleExprType >
struct SymbolsExpr : public SymbolsExprBase
{
    using self_type = SymbolsExpr<TupleExprType>;
    using symbols_expr_type = self_type;
    using tuple_type = TupleExprType;
    using feelpp_tag = SymbolsExprTag;

    SymbolsExpr() = default;
    SymbolsExpr( SymbolsExpr const& ) = default;
    SymbolsExpr( SymbolsExpr && ) = default;
    SymbolsExpr& operator=( SymbolsExpr const& ) = default;
    SymbolsExpr& operator=( SymbolsExpr&& ) = default;

    explicit SymbolsExpr( tuple_type const& tse )
        :
        M_tupleExpr( tse )
        {}

    explicit SymbolsExpr( tuple_type && tse )
        :
        M_tupleExpr( tse )
        {}


    std::map<std::string,std::set<std::string>> names() const
        {
            std::map<std::string,std::set<std::string>> res;
            hana::for_each( this->tuple(), [&res]( auto const& e )
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

    void updateInformationObject( nl::json & p ) const
        {
            nl::json::array_t ja;
            hana::for_each( this->tuple(), [&ja]( auto const& e )
                            {
                                for ( auto const& se : e )
                                {
                                    ja.push_back( symbolExprInformations( se.symbol(), "", se.componentSuffix() ) );
                                }
                            });
            p = ja;
        }

    tuple_type const& tuple() const { return M_tupleExpr; }
    tuple_type & tuple() { return M_tupleExpr; }

    template<typename... TheExpr>
    struct Lambda
    {
        struct TransformLambdaExpr
        {
            template <typename T>
            constexpr auto operator()(T const& t) const
                {
                    return typename T::template Lambda<TheExpr...>::type{};
                }
        };

        using lambda_tuple_type = std::decay_t<decltype( hana::transform( tuple_type{}, TransformLambdaExpr{} ) ) >;
        using type = SymbolsExpr<lambda_tuple_type>;
    };

    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    applyLambda( TheExpr... e ) const
        {
            return typename Lambda<TheExpr...>::type( hana::transform( M_tupleExpr, [&e...]( auto const& t) { return t.applyLambda(e...); } ) );
        }

    template <typename MapExprTensorType>
    struct TensorContext : public SymbolsExprTensorContextBase
    {
        using map_expr_tensor_type = MapExprTensorType;
        using symbols_expr_type = symbols_expr_type;

        TensorContext() = default;

        template <typename TheMetType, std::enable_if_t<std::is_same_v<std::decay_t<TheMetType>, map_expr_tensor_type>,bool> = true >
        TensorContext( std::shared_ptr<symbols_expr_type> const& se, TheMetType && met )
            :
            M_se( se ),
            M_mapExprTensor( std::forward<TheMetType>( met ) )
            {
#if 0
                hana::for_each( M_mapExprTensor, [this]( auto& pair )
                                {
                                    using Geo_t = typename std::decay_t<decltype(hana::first(pair))>::type;
                                    if (!hana::second(pair))
                                    {
                                        hana::second(pair) = std::make_shared<ExprTensorsFromSymbolsExprImpl<Geo_t/*,Basis_i_t,Basis_j_t*/,symbols_expr_type>>( *M_se );
                                    }
                                });
#endif
            }
        TensorContext( TensorContext const& ) = default;
        TensorContext( TensorContext && ) = default;
        TensorContext& operator=( TensorContext const& ) = default;
        TensorContext& operator=( TensorContext&& ) = default;

        symbols_expr_type const& symbolsExpression() const { return *M_se; }
        symbols_expr_type & symbolsExpression() { return *M_se; }
        map_expr_tensor_type const& mapExprTensor() const { return M_mapExprTensor; }

        tuple_type const& tuple() const { return M_se->tuple(); }
        tuple_type & tuple() { return M_se->tuple(); }

    private :
        std::shared_ptr<symbols_expr_type> M_se;
        map_expr_tensor_type M_mapExprTensor;
    };

    template <typename... MeshesType>
    auto createTensorContext() const
        {
            auto sePtr = std::make_shared<self_type>( *this );
            auto met = hana::fold( hana::tuple_t<MeshesType...>, hana::make_map(), [this,&sePtr]( auto && state, auto const& e )
                                   {
                                       using mesh_type = typename std::decay_t<decltype(e)>::type;
                                       using gm_type = typename mesh_type::gm_type;
                                       using geoelement_type = typename mesh_type::element_type;

                                       using state_type = std::decay_t<decltype(state)>;

                                       if constexpr ( geoelement_type::nOrder > 1 )
                                       {
                                           using gm1_type = typename mesh_type::gm1_type;
                                           auto met_gN = createTensorContextImpl<gm_type,geoelement_type>( sePtr );
                                           auto met_g1 = createTensorContextImpl<gm1_type,geoelement_type>( sePtr );
                                           return hana::union_( std::forward<state_type>( state ), hana::union_( std::move(met_gN), std::move(met_g1) ) );
                                       }
                                       else
                                       {
                                           auto met = createTensorContextImpl<gm_type,geoelement_type>( sePtr );
                                           return hana::union_( std::forward<state_type>( state ), std::move(met) );
                                       }
                                   } );
            using MapExprTensorType = std::decay_t<decltype(met)>;
            return TensorContext<MapExprTensorType>( sePtr, std::move( met ) );
        }

private :
    template <typename GmType,typename GeoElementType>
    auto
    createTensorContextImpl( std::shared_ptr<self_type> const& sePtr ) const
        {
            using gm_type = GmType;
            using geoelement_type = GeoElementType;

            auto buildExprTensor = [&sePtr]( auto const& theMapGmc )
            {
                using Geo_t = typename std::decay_t<decltype(theMapGmc)>::type;
                using etfse_type = ExprTensorsFromSymbolsExpr<Geo_t/*,Basis_i_t,Basis_j_t*/>;
                using etfse_ptrtype = std::shared_ptr<etfse_type>;
                etfse_ptrtype etfsei = std::make_shared<ExprTensorsFromSymbolsExprImpl<Geo_t/*,Basis_i_t,Basis_j_t*/,symbols_expr_type>>( *sePtr );
                return etfsei;
            };

            using gmc_type = typename gm_type::template Context<geoelement_type>;
            using Geo_element_t = map_gmc_type<gmc_type>;
            using Geo_face_t = map_gmc_type<typename gm_type::template Context<geoelement_type,1>>;
            using Geo_face2_t = map2_gmc_type<typename gm_type::template Context<geoelement_type,1>>;

            if constexpr ( geoelement_type::nDim == 3 )
            {
                using Geo_point_t = map_gmc_type<typename gm_type::template Context<geoelement_type,3>>;
                using Geo_edge_t = map_gmc_type<typename gm_type::template Context<geoelement_type,2>>;
                auto met = hana::make_map( hana::make_pair(hana::type_c<Geo_element_t>, buildExprTensor( hana::type_c<Geo_element_t> ) ),
                                           hana::make_pair(hana::type_c<Geo_face_t>, buildExprTensor( hana::type_c<Geo_face_t> ) ),
                                           hana::make_pair(hana::type_c<Geo_face2_t>, buildExprTensor( hana::type_c<Geo_face2_t> ) ),
                                           hana::make_pair(hana::type_c<Geo_edge_t>, buildExprTensor( hana::type_c<Geo_edge_t> ) ),
                                           hana::make_pair(hana::type_c<Geo_point_t>, buildExprTensor( hana::type_c<Geo_point_t> ) )
                                           );
                return met;
            }
            else if constexpr ( geoelement_type::nDim == 2 )
            {
                using Geo_point_t = map_gmc_type<typename gm_type::template Context<geoelement_type,2>>;
                auto met = hana::make_map( hana::make_pair(hana::type_c<Geo_element_t>, buildExprTensor( hana::type_c<Geo_element_t> ) ),
                                           hana::make_pair(hana::type_c<Geo_face_t>, buildExprTensor( hana::type_c<Geo_face_t> ) ),
                                           hana::make_pair(hana::type_c<Geo_face2_t>, buildExprTensor( hana::type_c<Geo_face2_t> ) ),
                                           hana::make_pair(hana::type_c<Geo_point_t>, buildExprTensor( hana::type_c<Geo_point_t> ) )
                                           );
                return met;
            }
            else
            {
                auto met = hana::make_map( hana::make_pair(hana::type_c<Geo_element_t>, buildExprTensor( hana::type_c<Geo_element_t> ) ),
                                           hana::make_pair(hana::type_c<Geo_face_t>, buildExprTensor( hana::type_c<Geo_face_t> ) ),
                                           hana::make_pair(hana::type_c<Geo_face2_t>, buildExprTensor( hana::type_c<Geo_face2_t> ) )
                                           );
                return met;
            }
        }

private :
    tuple_type M_tupleExpr;
};

template<typename SymbolsExprType>
using is_symbols_expression = typename std::is_base_of<SymbolsExprBase,SymbolsExprType>::type;

template<typename SymbolsExprType>
constexpr bool is_symbols_expression_v = is_symbols_expression<SymbolsExprType>::value;

template<typename SymbolsExprType>
using is_symbols_expression_tensor_context = typename std::is_base_of<SymbolsExprTensorContextBase,SymbolsExprType>::type;

template<typename SymbolsExprType>
constexpr bool is_symbols_expression_tensor_context_v = is_symbols_expression_tensor_context<SymbolsExprType>::value;



template<typename T1,typename... ExprT>
struct SymbolsExprTraits
{
    static constexpr auto callApply = [](const auto& ...exprs) { return Feel::detail::AdvancedConcatOfTupleContainerType<SymbolsExprTag,SymbolExprTag>::template apply( exprs... ); };
    using tuple_type = std::decay_t<decltype( hana::unpack( hana::tuple<T1,ExprT...>{},  callApply ) )>;
    using type = SymbolsExpr<tuple_type>;
};


template<typename T>
struct SymbolsExprTraits<T>
{
    struct DeferToSymbolsExpr
    {
        template<class Ts>
        using execute = Ts;
    };
    struct DeferToSymbolExpr
    {
        template<class Ts>
        using execute = SymbolsExpr< hana::tuple<Ts> >;
    };
    template< class Prog, class... Ts >
    using run = typename Prog::template execute<Ts...>;

    using choice = typename std::conditional< is_symbols_expression_v<T> || is_symbols_expression_tensor_context_v<T>,
                                              DeferToSymbolsExpr,
                                              DeferToSymbolExpr>::type;
    using type = run< choice, T >;
};

template<typename... ExprT>
using symbols_expression_t = typename SymbolsExprTraits<ExprT...>::type;

using symbols_expression_empty_t = SymbolsExpr<hana::tuple<>>;

//template<typename SymbolsExprType>
//using is_symbols_expression_empty = typename std::is_same_vbase_of<SymbolsExprBase,SymbolsExprType>::type;

template<typename SymbolsExprType>
constexpr bool is_symbols_expression_empty_v = std::is_same_v<SymbolsExprType,symbols_expression_empty_t>;

//! build a SymbolsExpr object
template<typename... ExprT>
symbols_expression_t<ExprT...>
symbolsExpr( const ExprT&... exprs )
{
    return symbols_expression_t<ExprT...>(Feel::detail::AdvancedConcatOfTupleContainerType<SymbolsExprTag,SymbolExprTag>::template apply( exprs... ) );
}
template<typename T>
symbols_expression_t<T> const&
symbolsExpr( T const& se, std::enable_if_t< (is_symbols_expression_v<T> || is_symbols_expression_tensor_context_v<T>) >* = nullptr )
{
    return se;
}

} // namespace vf

} // namespace Feel

#endif
