/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2019-11-16

  Copyright (C) 2019 Feel++ Consortium

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
#ifndef FEELPP_MODELS_VF_ExprSelectorByMeshElement_H
#define FEELPP_MODELS_VF_ExprSelectorByMeshElement_H 1

namespace Feel
{
namespace vf
{

template <typename IndexType>
class ExprSelectorByMeshElementMapping
{
public :
    using tag_type = uint16_type;
    ExprSelectorByMeshElementMapping() = default;

    template <typename MeshType>
    void updateForUse( std::map<std::string, std::tuple<elements_reference_wrapper_t<MeshType>,elements_reference_wrapper_t<MeshType>> > const& data )
        {
            uint16_type cpt=0;
            for ( auto const& [name,pairRangeElt] : data )
            {
                for ( auto const& eltWrap : std::get<0>( pairRangeElt ) ) // active elements
                    M_eltIdToTag[ unwrap_ref( eltWrap ).id() ] = cpt;
                for ( auto const& eltWrap : std::get<1>( pairRangeElt ) ) // ghost elements
                    M_eltIdToTag[ unwrap_ref( eltWrap ).id() ] = cpt;
                M_nameToTag[name] = cpt++;
            }
        }

    tag_type idToTag( IndexType id ) const
        {
            auto itFindId = M_eltIdToTag.find( id );
            if ( itFindId != M_eltIdToTag.end() )
                return itFindId->second;
            return invalid_v<tag_type>;
        }

    tag_type nameToTag( std::string const& name ) const
        {
            auto itFindName = M_nameToTag.find( name );
            if ( itFindName != M_nameToTag.end() )
                return itFindName->second;
            return invalid_v<tag_type>;
        }

    void clear()
        {
            M_eltIdToTag.clear();
            M_nameToTag.clear();
        }
private :
    std::unordered_map<index_type,tag_type> M_eltIdToTag;
    std::map<std::string,tag_type> M_nameToTag;
};

template <typename IndexType, typename TupleVectorExprType>
class ExprSelectorByMeshElement : public Feel::vf::ExprDynamicBase
{
public :
    using this_type = ExprSelectorByMeshElement<IndexType,TupleVectorExprType>;
    using mapping_type = ExprSelectorByMeshElementMapping<IndexType>;
    using mapping_ptrtype = std::shared_ptr<mapping_type>;
    using mapping_tag_type = typename mapping_type::tag_type;

    using tuple_vector_expr_type = TupleVectorExprType;
    static const int nExpr = std::decay_t<decltype(hana::size(tuple_vector_expr_type{}))>::value;

    using first_expression_type = typename std::decay_t<decltype(hana::at_c<0>( tuple_vector_expr_type{} ))>::value_type::second_type;

    static const size_type context = vm::DYNAMIC;
    static const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = std::decay_t<decltype( hana::fold( tuple_vector_expr_type{},
                                                                      hana::integral_constant<bool,false>{},
                                                                      typename GinacExVF<>::FunctorsVariadicExpr::template HasTestFunction<Func>{} ) )>::value;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result =  std::decay_t<decltype( hana::fold( tuple_vector_expr_type{},
                                                                       hana::integral_constant<bool,false>{},
                                                                       typename GinacExVF<>::FunctorsVariadicExpr::template HasTrialFunction<Func>{} ) )>::value;
    };

    template<typename Funct>
    static const bool has_test_basis = std::decay_t<decltype( hana::fold( tuple_vector_expr_type{},
                                                                          hana::integral_constant<bool,false>{},
                                                                          typename GinacExVF<>::FunctorsVariadicExpr::template HasTestBasis<Funct>{} ) )>::value;
    template<typename Funct>
    static const bool has_trial_basis = std::decay_t<decltype( hana::fold( tuple_vector_expr_type{},
                                                                           hana::integral_constant<bool,false>{},
                                                                           typename GinacExVF<>::FunctorsVariadicExpr::template HasTrialBasis<Funct>{} ) )>::value;
    using test_basis = std::nullptr_t;//TODO//typename expression_type::test_basis;
    using trial_basis = std::nullptr_t;//TODO//typename expression_type::trial_basis;
    typedef typename first_expression_type::value_type value_type;
    typedef typename first_expression_type::evaluate_type evaluate_type;

    ExprSelectorByMeshElement( mapping_ptrtype mapping, tuple_vector_expr_type const& exprs )
        :
        M_mapping( mapping ),
        M_exprs( exprs )
        {}

    ExprSelectorByMeshElement( ExprSelectorByMeshElement const& ) = default;

    //! polynomial order
    uint16_type polynomialOrder() const
        {
            uint16_type res = 0;
            hana::for_each( M_exprs, [&res]( auto const& e )
                            {
                                for ( auto const& [name,expr] : e )
                                    res = std::max( res, expr.polynomialOrder() );
                            });
            return res;
        }

    //! expression is polynomial?
    bool isPolynomial() const
        {
            bool res = true;
            hana::for_each( M_exprs, [&res]( auto const& e )
                            {
                                for ( auto const& [name,expr] : e )
                                    res = res && expr.isPolynomial();
                            });
            return res;
        }

    size_type dynamicContext() const
        {
            size_type res = 0;
            hana::for_each( M_exprs, [&res]( auto const& e )
                            {
                                for ( auto const& [name,expr] : e )
                                    res = res | Feel::vf::dynamicContext( expr );
                            });
            return res;
        }

    mapping_type const& mapping() const { return *M_mapping; }
    tuple_vector_expr_type const& tupleExpressions() const { return M_exprs; }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            auto newTupleTensorExprs = hana::transform( this->tupleExpressions(), [&se](auto const& t)
                                                        {
                                                            using new_expr_type = std::decay_t<decltype( std::get<1>( t.front() ).applySymbolsExpr( se ) )>;
                                                            std::vector<std::pair<std::string,new_expr_type>> res;
                                                            for ( auto const&  [name,theexpr] : t )
                                                                res.push_back( std::make_pair(name, theexpr.applySymbolsExpr( se ) ) );
                                                            return res;
                                                        } );
            return ExprSelectorByMeshElement<IndexType, std::decay_t<decltype(newTupleTensorExprs)> >( M_mapping,newTupleTensorExprs );
        }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            bool res = false;
            hana::for_each( M_exprs, [&res,&symb,&se]( auto const& e )
                            {
                                if ( res )
                                    return;
                                for ( auto const& [name,expr] : e )
                                {
                                    res = expr.hasSymbolDependency(symb,se);
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
                                for ( auto const& [name,expr] : e )
                                    expr.dependentSymbols(symb,res,se);
                            });
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        auto newTupleTensorExprs = hana::transform( this->tupleExpressions(), [&diffVariable,&world,&dirLibExpr,&se](auto const& t)
                                                    {
                                                        using new_expr_type = std::decay_t<decltype( std::get<1>( t.front() ).template diff<diffOrder>( diffVariable,world,dirLibExpr,se ) )>;
                                                        std::vector<std::pair<std::string,new_expr_type>> res;
                                                        for ( auto const&  [name,theexpr] : t )
                                                            res.push_back( std::make_pair(name, theexpr.template diff<diffOrder>( diffVariable,world,dirLibExpr,se ) ) );
                                                        return res;
                                                    } );
        return ExprSelectorByMeshElement<IndexType, std::decay_t<decltype(newTupleTensorExprs)> >( M_mapping,newTupleTensorExprs );
    }


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        struct TransformExprToTensor
        {
            explicit TransformExprToTensor( mapping_type const& mapping ) : M_mapping( mapping ) {}

            template <typename T>
            struct apply {
                using type = typename T::value_type::second_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
            };

            template <typename T>
            constexpr auto operator()(T const& t) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    std::map<uint16_type,_tensor_type> res;
                    for ( auto const&  [name,theexpr]  : t )
                        res.emplace( M_mapping.nameToTag( name ), _tensor_type( theexpr,Geo_t{} ) );
                    _tensor_type * res2 = nullptr;
                    return std::make_pair(res,res2);
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    std::map<uint16_type,_tensor_type> res;
                    for ( auto const&  [name,theexpr]  : t )
                        res.emplace( M_mapping.nameToTag( name ), _tensor_type( theexpr,geom,fev,feu  ) );
                    _tensor_type * res2 = nullptr;
                    return std::make_pair(res,res2);
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom, Basis_i_t const& fev ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    std::map<uint16_type,_tensor_type> res;
                    for ( auto const&  [name,theexpr]  : t )
                        res.emplace( M_mapping.nameToTag( name ), _tensor_type( theexpr,geom,fev ) );
                    _tensor_type * res2 = nullptr;
                    return std::make_pair(res,res2);
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    std::map<uint16_type,_tensor_type> res;
                    for ( auto const&  [name,theexpr]  : t )
                        res.emplace( M_mapping.nameToTag( name ), _tensor_type( theexpr,geom ) );
                    _tensor_type * res2 = nullptr;
                    return std::make_pair(res,res2);
                }
            template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename T, typename... TheArgsType>
            constexpr auto operator()(std::true_type, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                                      T const& t, Geo_t const& geom, const TheArgsType&... theInitArgs ) const
                {
                    using _tensor_type = typename TransformExprToTensor::template apply<T>::type;
                    std::map<uint16_type,_tensor_type> res;
                    int k=0;
                    for ( auto const&  [name,theexpr]  : t )
                    {
                        CHECK( name == exprExpanded[k].first ) << "incompatible name";
                        res.emplace( M_mapping.nameToTag( name ), _tensor_type( std::true_type{}, exprExpanded[k++].second, ttse, theexpr, geom, theInitArgs... ) );
                    }

                    _tensor_type * res2 = nullptr;
                    return std::make_pair(res,res2);
                }


            mapping_type const& M_mapping;
        };

        using tuple_tensor_expr_type = std::decay_t<decltype( hana::transform( tuple_vector_expr_type{}, TransformExprToTensor{mapping_type{}} ) ) >;

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
            M_mapping( expr.mapping() ),
            M_tupleTensorExprs( hana::transform( expr.tupleExpressions(), [this,&geom,&fev,&feu](auto const& t) { return TransformExprToTensor{M_mapping}(t,geom,fev,feu); } ) )
            {}
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_mapping( expr.mapping() ),
            M_tupleTensorExprs( hana::transform( expr.tupleExpressions(), [this,&geom,&fev](auto const& t) { return TransformExprToTensor{M_mapping}(t,geom,fev); } ) )
            {}
        tensor( this_type const& expr,
                Geo_t const& geom )
            :
            M_mapping( expr.mapping() ),
            M_tupleTensorExprs( hana::transform( expr.tupleExpressions(), [this,&geom](auto const& t) { return TransformExprToTensor{M_mapping}(t,geom); } ) )
            {}
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_mapping( exprExpanded.mapping() ),
            M_tupleTensorExprs( hana::transform( hana::unpack( hana::make_range( hana::int_c<0>, hana::int_c<nExpr> ), hana::make_tuple ), [this,&exprExpanded,&ttse,&expr,&geom,&theInitArgs...](auto eId )
                                                 {
                                                     return TransformExprToTensor{M_mapping}(std::true_type{},hana::at( exprExpanded.tupleExpressions(), hana::int_c<eId> ), ttse,
                                                                                             hana::at( expr.tupleExpressions(), hana::int_c<eId> ), geom, theInitArgs...);
                                                 } ) )
            {}



        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            this->selectSubTensor( geom );
            hana::for_each( M_tupleTensorExprs, [&geom,&fev,&feu]( auto & e )
                            {
                                if ( e.second )
                                    e.second->update( geom, fev, feu );
                            });
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            this->selectSubTensor( geom );
            hana::for_each( M_tupleTensorExprs, [&geom,&fev]( auto & e )
                            {
                                if ( e.second )
                                    e.second->update( geom, fev );
                            });
        }
        void update( Geo_t const& geom )
        {
            this->selectSubTensor( geom );
            hana::for_each( M_tupleTensorExprs, [&geom]( auto & e )
                            {
                                if ( e.second )
                                    e.second->update( geom );
                            });
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            mapping_tag_type tag = this->selectSubTensor( geom );
            hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nExpr> ), [this,&tag,&exprExpanded,&ttse,&geom,&theUpdateArgs...]( auto eId )
                            {
                                auto & subtensorPtr = hana::at( M_tupleTensorExprs, hana::int_c<eId> ).second;
                                if ( subtensorPtr )
                                {
                                    auto const& exprExpandedCurrent = hana::at( exprExpanded.tupleExpressions(), hana::int_c<eId> );
                                    auto itFindExpr = std::find_if( exprExpandedCurrent.begin(), exprExpandedCurrent.end(), [this,&tag]( auto const& ee ) { return M_mapping.nameToTag( ee.first ) == tag; } );
                                    CHECK( itFindExpr != exprExpandedCurrent.end() ) << "expr from tag not find : " << tag;
                                    subtensorPtr->update( std::true_type{}, itFindExpr->second, ttse, geom, theUpdateArgs... );
                                }
                            } );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res(0);
            hana::for_each( M_tupleTensorExprs, [&i,&j,&c1,&c2,&q,&res]( auto const& e )
                            {
                                if ( e.second )
                                    res = e.second->evalijq( i, j, c1, c2, q );
                            });
            return res;
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res(0);
            hana::for_each( M_tupleTensorExprs, [&i,&c1,&c2,&q,&res]( auto const& e )
                            {
                                if ( e.second )
                                    res = e.second->evaliq( i, c1, c2, q );
                            });
            return res;
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            value_type res(0);
            hana::for_each( M_tupleTensorExprs, [&c1,&c2,&q,&res]( auto const& e )
                            {
                                if ( e.second )
                                    res = e.second->evalq( c1, c2, q );
                            });
            return res;
        }

    private :
        mapping_tag_type selectSubTensor( Geo_t const& geom )
            {
                // first, reset all current tensors
                hana::for_each( M_tupleTensorExprs, []( auto & e ) { e.second = nullptr; } );

                IndexType eid = vf::detail::ExtractGm<Geo_t>::get( geom )->id();
                mapping_tag_type tag = M_mapping.idToTag( eid );
                if ( tag != invalid_v<mapping_tag_type> )
                {
                    hana::for_each( M_tupleTensorExprs, [&tag]( auto & e )
                                    {
                                        //e.second = nullptr;
                                        auto & tensorExprs = e.first;
                                        auto itFindTensorExpr = tensorExprs.find( tag );
                                        if ( itFindTensorExpr !=  tensorExprs.end() )
                                            e.second = &itFindTensorExpr->second;
                                    });
                }
                return tag;
            }
    private :
        mapping_type const& M_mapping;
        tuple_tensor_expr_type M_tupleTensorExprs;
    };

private :
    mapping_ptrtype M_mapping;
    TupleVectorExprType M_exprs;
};


template <typename IndexType, typename VectorExpr1Type, typename ... VectorExprOtherType>
inline
Expr<ExprSelectorByMeshElement<IndexType,hana::tuple<VectorExpr1Type,VectorExprOtherType...>>>
    expr( typename ExprSelectorByMeshElement<IndexType,hana::tuple<VectorExpr1Type, VectorExprOtherType...>>::mapping_ptrtype mapping, VectorExpr1Type const& exprs1, VectorExprOtherType ... exprsOther )
{
    typedef ExprSelectorByMeshElement<IndexType,hana::tuple<VectorExpr1Type, VectorExprOtherType...>> esbme_t;
    return Expr< esbme_t >( esbme_t( mapping, hana::make_tuple( exprs1, exprsOther... ) ) );
}

} // namespace vf
} // namespace Feel

#endif
