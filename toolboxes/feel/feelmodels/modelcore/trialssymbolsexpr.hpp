/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_MODELCORE_TRIALSSYMBOLSEXPR_H
#define FEELPP_TOOLBOXES_MODELCORE_TRIALSSYMBOLSEXPR_H 1

#include <feel/feelcore/tuple_utils.hpp>
#include <feel/feelmodels/modelcore/traits.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename SpaceType,typename ExprType>
class TrialSymbolExpr1
{
public :
    using space_type = SpaceType;
    using space_ptrtype = std::shared_ptr<space_type>;
    using expr_type = ExprType;
    using expr_shape_type = typename ExprTraitsFromMeshElement<typename space_type::mesh_type::element_type,expr_type>::shape;

    TrialSymbolExpr1( space_ptrtype const& space, std::string const& symbol, expr_type const& expr, size_type blockSpaceIndex )
        :
        M_space( space ),
        M_symbol( symbol ),
        M_expr( expr ),
        M_blockSpaceIndex( blockSpaceIndex ),
        M_secs( expr_shape_type::M, expr_shape_type::N )
        {}
    TrialSymbolExpr1( TrialSymbolExpr1 const& ) = default;
    TrialSymbolExpr1( TrialSymbolExpr1 && ) = default;

    space_ptrtype const& space() const { return M_space; }
    std::string const& symbol() const { return M_symbol; }
    expr_type const& expr() const { return M_expr; }
    size_type blockSpaceIndex() const { return M_blockSpaceIndex; }
    SymbolExprComponentSuffix const& componentSuffix() const { return M_secs; }

    bool isDefinedOn( space_ptrtype const& space, size_type blockSpaceIndex ) const
        {
            if ( M_space != space )
                return false;
            if ( M_blockSpaceIndex != blockSpaceIndex )
                return false;
            return true;
        }
private :
    space_ptrtype M_space;
    std::string M_symbol;
    expr_type M_expr;
    size_type M_blockSpaceIndex;
    SymbolExprComponentSuffix M_secs;
};


struct TrialSymbolExprFeelppTag {};
struct TrialSymbolsExprFeelppTag {};
struct TrialsSymbolsExprFeelppTag {};

template <typename SpaceType,typename ExprType>
class TrialSymbolExpr : public std::vector<TrialSymbolExpr1<SpaceType,ExprType>>
{
    using super_type = std::vector<TrialSymbolExpr1<SpaceType,ExprType>>;
public :
    using trial_symbol_expr1_type = TrialSymbolExpr1<SpaceType,ExprType>;
    using space_type = typename trial_symbol_expr1_type::space_type;
    using space_ptrtype = typename trial_symbol_expr1_type::space_ptrtype;
    using expr_type = typename trial_symbol_expr1_type::expr_type;

    using feelpp_tag = TrialSymbolExprFeelppTag;

    TrialSymbolExpr() = default;
    TrialSymbolExpr( space_ptrtype const& space, std::string const& symbol, expr_type const& expr, size_type blockSpaceIndex ) : super_type( 1, trial_symbol_expr1_type( space, symbol, expr, blockSpaceIndex ) ) {}
#if 1
    TrialSymbolExpr( TrialSymbolExpr const& ) = default;
#else
    TrialSymbolExpr( TrialSymbolExpr const& tse ) : super_type( tse ) { std::cout << "TrialSymbolExpr copy constructor with size " << this->size() << std::endl; }
#endif
    TrialSymbolExpr( TrialSymbolExpr && ) = default;

    void add( space_ptrtype const& space, std::string const& symbol, expr_type const& expr, size_type blockSpaceIndex )
    {
        this->push_back( trial_symbol_expr1_type( space, symbol, expr, blockSpaceIndex ) );
    }
};

template <typename SpaceType,typename ExprType>
using trial_symbol_expr_t = TrialSymbolExpr<SpaceType,ExprType>;

template <typename SpaceType,typename ExprType>
trial_symbol_expr_t<SpaceType,ExprType>
trialSymbolExpr( std::shared_ptr<SpaceType> const& space, std::string const& symbol, ExprType const& expr, size_type blockSpaceIndex )
{
    return trial_symbol_expr_t<SpaceType,ExprType>( space, symbol, expr, blockSpaceIndex );
}


template <typename SpaceType,typename TupleTrialSymbolExprType>
class TrialSymbolsExpr
{
public :
    using feelpp_tag = TrialSymbolsExprFeelppTag;
    using tuple_type = TupleTrialSymbolExprType;
    using space_type = SpaceType;
    using space_ptrtype = std::shared_ptr<space_type>;

    TrialSymbolsExpr() = default;
    TrialSymbolsExpr( TrialSymbolsExpr const& ) = default;
    TrialSymbolsExpr( TrialSymbolsExpr && ) = default;
    TrialSymbolsExpr& operator=( TrialSymbolsExpr&& ) = default;

    explicit TrialSymbolsExpr( tuple_type const& t )
        :
        M_tuple( t )
        {
            this->updateForUse();
        }

    explicit TrialSymbolsExpr( tuple_type && t )
        :
        M_tuple( t )
        {
            this->updateForUse();
        }

    tuple_type const& tuple() const { return M_tuple; }

    std::map<size_type,space_ptrtype> const& blockSpaceIndex() const { return M_blockSpaceIndex; }
private :
    void updateForUse()
        {
            hana::for_each( M_tuple, [this]( auto const& e ) {
                    for ( auto const & e2 : e )
                    {
                        auto [it,isInserted] = M_blockSpaceIndex.insert( std::make_pair( e2.blockSpaceIndex(), e2.space() ) );
                        if ( !isInserted )
                            CHECK( it->second == e2.space() ) << "not same space at block index" << e2.blockSpaceIndex();
                    }
                });
        }
private :
    tuple_type M_tuple;
    std::map<size_type,space_ptrtype> M_blockSpaceIndex;
};

template<typename SpaceType,typename... TrialSymbolsExprType>
struct TrialSymbolsExprTraits
{
    static constexpr auto callApply = [](const auto& ...args) { return Feel::detail::AdvancedConcatOfTupleContainerType<TrialSymbolsExprFeelppTag,TrialSymbolExprFeelppTag>::template apply( args... ); };
    using tuple_type = std::decay_t<decltype( hana::unpack( hana::tuple<TrialSymbolsExprType...>{},  callApply ) )>;
    // TODO check that all space in tuple_type is equal to SpaceType
    using type = TrialSymbolsExpr<SpaceType,tuple_type>;
};

template<typename SpaceType,typename... TrialSymbolsExprType>
using trial_symbols_expr_t = typename TrialSymbolsExprTraits<SpaceType,TrialSymbolsExprType...>::type;

template<typename SpaceType>
using trial_symbols_expr_empty_t = TrialSymbolsExpr<SpaceType,hana::tuple<>>;

template<typename SpaceType,typename... TrialSymbolsExprType>
trial_symbols_expr_t<SpaceType,TrialSymbolsExprType...>
trialSymbolsExpr( const TrialSymbolsExprType&... tse )
{
    return trial_symbols_expr_t<SpaceType,TrialSymbolsExprType...>( Feel::detail::AdvancedConcatOfTupleContainerType<TrialSymbolsExprFeelppTag,TrialSymbolExprFeelppTag>::template apply( tse... ) );
}

template <typename FeelppTagOfMapType,typename FeelppTagOfContainerType>
struct AdvancedConcatOfMapContainerType
{
    template <typename Tag, typename T>
    struct is_a_t :
        hana::integral_constant<bool, std::is_same<Tag, typename T::feelpp_tag >::value >
    {};

    template<typename... SomeType>
    static constexpr auto apply( const SomeType&... v )
        {
            return applyImpl( hana::map<>{}, v... );
        }
private :

    template < typename T1, typename ResType, typename... SomeType >
    static constexpr auto applyFromContainerType( T1 const& t1, ResType && res, hana::basic_tuple<SomeType...> && keysInRes )
        {
            auto constexpr the_key = hana::type_c<typename T1::space_type>;
            if constexpr ( hana::find( hana::to_tuple(hana::tuple/*_t*/<SomeType...>{} ), hana::type_c<typename T1::space_type>) == hana::nothing )
            {
                return hana::insert( std::forward<ResType>( res ), hana::make_pair(the_key, t1 ) );
            }
            else
            {
                auto newEntry = trialSymbolsExpr<typename T1::space_type>( hana::at_key( res, hana::type_c<typename T1::space_type> ), t1 );
                return hana::insert( hana::erase_key( std::forward<ResType>( res ), the_key ),  hana::make_pair( the_key, std::move( newEntry ) ) );
            }
        }

    template<int Index, typename ResType, typename T1, typename... SomeType >
    static constexpr auto applyImpl2( ResType && res, T1 const& t1, hana::basic_tuple<SomeType...> && keysInT1 )
        {
            constexpr int nKeys = std::decay_t<decltype(hana::size( hana::basic_tuple<SomeType...> {} ))>::value;
            if constexpr ( Index < nKeys )
                return applyImpl2<Index+1>( applyImpl( std::forward<ResType>( res ), hana::at_key( t1, hana::at( keysInT1, hana::int_c<Index> ) ) ),
                                            t1,  std::forward< hana::basic_tuple<SomeType...> >( keysInT1 ) );
            else
                return std::move( res );
        }

    template<typename ResType >
    static constexpr auto applyImpl( ResType && res )
        {
            return std::move( res );
        }
    template<typename ResType, typename T1, typename... SomeType2 >
    static constexpr auto applyImpl( ResType && res, T1 const& t1, const SomeType2&... tothers )
        {
            if constexpr ( is_a_t<FeelppTagOfContainerType, T1 >::value )
                {
                    return applyImpl( applyFromContainerType( t1,std::forward<ResType>(res), hana::keys(res) ), tothers... );
                }
            else if constexpr ( is_a_t<FeelppTagOfMapType, T1 >::value )
                {
                    return applyImpl( applyImpl2<0>( std::forward<ResType>( res ), t1.map(), hana::keys(t1.map()) ), tothers... );
                }
            else
                return std::move( res );
        }
};

template <typename MapTrialSymbolsExprType>
class TrialsSymbolsExpr
{
public :
    using feelpp_tag = TrialsSymbolsExprFeelppTag;
    using map_type = MapTrialSymbolsExprType;

    TrialsSymbolsExpr() = default;
    TrialsSymbolsExpr( TrialsSymbolsExpr const& ) = default;
    TrialsSymbolsExpr( TrialsSymbolsExpr && ) = default;
    TrialsSymbolsExpr& operator=( TrialsSymbolsExpr&& ) = default;

    explicit TrialsSymbolsExpr( map_type const& m )
        :
        M_map( m )
        {}

    explicit TrialsSymbolsExpr( map_type && m )
        :
        M_map( m )
        {}

    map_type const& map() const { return M_map; }


    std::set<std::string> names() const
        {
            std::set<std::string> res;
            hana::for_each( M_map, [&res]( auto const& e )
                            {
                                hana::for_each( hana::second(e).tuple(), [&res]( auto const& e2 ) {
                                        for ( auto const& e3 : e2 )
                                        {
                                            if constexpr ( std::decay_t<decltype( e3 )>::expr_shape_type::is_scalar )
                                            {
                                                res.insert( e3.symbol() );
                                            }
                                            else
                                            {
                                                for ( auto const& [_suffix,compArray] : e3.componentSuffix() )
                                                    res.insert( e3.symbol()+_suffix );
                                            }
                                        }
                                    });
                            });
            return res;
        }

private :
    map_type M_map;
};

template<typename... TrialSymbolsExprType>
struct TrialsSymbolsExprTraits
{
    static constexpr auto callApply = [](const auto& ...args) { return Feel::FeelModels::AdvancedConcatOfMapContainerType<TrialsSymbolsExprFeelppTag,TrialSymbolsExprFeelppTag>::template apply( args... ); };
    using map_type = std::decay_t<decltype( hana::unpack( hana::tuple<TrialSymbolsExprType...>{},  callApply ) )>;
    using type = TrialsSymbolsExpr<map_type>;
};

template<typename... TrialSymbolsExprType>
using trials_symbols_expr_t = typename TrialsSymbolsExprTraits<TrialSymbolsExprType...>::type;

using trials_symbols_expr_empty_t = TrialsSymbolsExpr<hana::map<>>;

template<typename... TrialSymbolsExprType>
trials_symbols_expr_t<TrialSymbolsExprType...>
trialsSymbolsExpr( const TrialSymbolsExprType&... tse )
{
    return trials_symbols_expr_t<TrialSymbolsExprType...>( Feel::FeelModels::AdvancedConcatOfMapContainerType<TrialsSymbolsExprFeelppTag,TrialSymbolsExprFeelppTag>::template apply( tse... ) );
}


} // namespace FeelModels
} // namespace Feel

#endif
