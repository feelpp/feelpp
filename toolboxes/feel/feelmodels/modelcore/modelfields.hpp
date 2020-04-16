/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_MODELCORE_MODELFIELDS_H
#define FEELPP_TOOLBOXES_MODELCORE_MODELFIELDS_H 1

#include <type_traits>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <feel/feelcore/feeltypes.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelvf/symbolsexpr.hpp>
#include <feel/feelvf/inner.hpp>

namespace Feel
{
namespace FeelModels
{

namespace FieldCtx
{

const size_type ID           = ( 1<<0 );
const size_type MAGNITUDE    = ( 1<<1 );
const size_type GRAD         = ( 1<<2 );
const size_type GRAD_NORMAL  = ( 1<<3 );

}




template <typename SpaceType,typename ExprType>
class TrialSymbolExpr1
{
public :
    using space_type = SpaceType;
    using space_ptrtype = std::shared_ptr<space_type>;
    using expr_type = ExprType;
    TrialSymbolExpr1( space_ptrtype const& space, std::string const& symbol, expr_type const& expr, size_type startBlockSpaceIndex )
        :
        M_space( space ),
        M_symbol( symbol ),
        M_expr( expr ),
        M_startBlockSpaceIndex( startBlockSpaceIndex )
        {}
    TrialSymbolExpr1( TrialSymbolExpr1 const& ) = default;
    TrialSymbolExpr1( TrialSymbolExpr1 && ) = default;

    space_ptrtype const& space() const { return M_space; }
    std::string const& symbol() const { return M_symbol; }
    expr_type const& expr() const { return M_expr; }
    size_type startBlockSpaceIndex() const { return M_startBlockSpaceIndex; }
private :
    space_ptrtype M_space;
    std::string M_symbol;
    expr_type M_expr;
    size_type M_startBlockSpaceIndex;
};


struct TrialSymbolExprFeelppTag {};
struct TrialSymbolsExprFeelppTag {};

template <typename SpaceType,typename ExprType>
class TrialSymbolExpr : public std::vector<TrialSymbolExpr1<SpaceType,ExprType>>
{
    using super_type = std::vector<TrialSymbolExpr1<SpaceType,ExprType>>;
    using trial_symbol_expr1_type = TrialSymbolExpr1<SpaceType,ExprType>;
    using space_ptrtype = typename trial_symbol_expr1_type::space_ptrtype;
    using expr_type = typename trial_symbol_expr1_type::expr_type;
public :
    using feelpp_tag = TrialSymbolExprFeelppTag;

    TrialSymbolExpr() = default;
    TrialSymbolExpr( space_ptrtype const& space, std::string const& symbol, expr_type const& expr, size_type startBlockSpaceIndex ) : super_type( 1, trial_symbol_expr1_type( space, symbol, expr, startBlockSpaceIndex ) ) {}
    TrialSymbolExpr( TrialSymbolExpr const& ) = default;
    TrialSymbolExpr( TrialSymbolExpr && ) = default;

    void add( space_ptrtype const& space, std::string const& symbol, expr_type const& expr, size_type startBlockSpaceIndex )
    {
        this->push_back( trial_symbol_expr1_type( space, symbol, expr, startBlockSpaceIndex ) );
    }
};

template <typename SpaceType,typename ExprType>
using trial_symbol_expr_t = TrialSymbolExpr<SpaceType,ExprType>;

template <typename SpaceType,typename ExprType>
trial_symbol_expr_t<SpaceType,ExprType>
trialSymbolExpr( std::shared_ptr<SpaceType> const& space, std::string const& symbol, ExprType const& expr, size_type startBlockSpaceIndex )
{
    return trial_symbol_expr_t<SpaceType,ExprType>( space, symbol, expr, startBlockSpaceIndex );
}


template <typename TupleTrialSymbolExprType>
class TrialSymbolsExpr
{
public :
    using feelpp_tag = TrialSymbolsExprFeelppTag;
    using tuple_type = TupleTrialSymbolExprType;

    TrialSymbolsExpr() = default;
    TrialSymbolsExpr( TrialSymbolsExpr const& ) = default;
    TrialSymbolsExpr( TrialSymbolsExpr && ) = default;
    TrialSymbolsExpr& operator=( TrialSymbolsExpr&& ) = default;

    explicit TrialSymbolsExpr( tuple_type const& t )
        :
        M_tuple( t )
        {}

    explicit TrialSymbolsExpr( tuple_type && t )
        :
        M_tuple( t )
        {}

    tuple_type const& tuple() const { return M_tuple; }
private :
    tuple_type M_tuple;
};

template<typename... TrialSymbolsExprType>
struct TrialSymbolsExprTraits
{
    static constexpr auto callApply = [](const auto& ...args) { return Feel::detail::AdvancedConcatOfTupleContainerType<TrialSymbolsExprFeelppTag,TrialSymbolExprFeelppTag>::template apply( args... ); };
    using tuple_type = std::decay_t<decltype( hana::unpack( hana::tuple<TrialSymbolsExprType...>{},  callApply ) )>;
    using type = TrialSymbolsExpr<tuple_type>;
};

template<typename... TrialSymbolsExprType>
using trial_symbols_expr_t = typename TrialSymbolsExprTraits<TrialSymbolsExprType...>::type;

using trial_symbols_expr_empty_t = TrialSymbolsExpr<hana::tuple<>>;

template<typename... TrialSymbolsExprType>
trial_symbols_expr_t<TrialSymbolsExprType...>
trialSymbolExprs( const TrialSymbolsExprType&... tse )
{
    return trial_symbols_expr_t<TrialSymbolsExprType...>( Feel::detail::AdvancedConcatOfTupleContainerType<TrialSymbolsExprFeelppTag,TrialSymbolExprFeelppTag>::template apply( tse... ) );
}








struct ModelFieldFeelppTag {};

template <typename ToolboxType,uint16_type TheFieldTag>
struct ModelFieldTag
{
    using toolbox_type = ToolboxType;
    static constexpr uint16_type field_tag = TheFieldTag;

    explicit ModelFieldTag( toolbox_type const* t ) : M_toolbox( t ) {}
    ModelFieldTag( ModelFieldTag const& ) = default;
    ModelFieldTag( ModelFieldTag && ) = default;

    //bool operator==( ModelFieldTag const& o ) const { return ( dynamic_cast<void const*>( this->M_toolbox ) == dynamic_cast<void const*>( o.M_toolbox ) ); }
    bool operator==( ModelFieldTag const& o ) const { return ( dynamic_cast<toolbox_type const*>( this->M_toolbox ) == dynamic_cast<toolbox_type const*>( o.M_toolbox ) ); }
    bool operator!=( ModelFieldTag const& o ) const { return !this->operator==( o ); }
private :
    toolbox_type const* M_toolbox;
};

template <typename ModelFieldTagType,typename FieldType>
struct ModelField1
{
    using self_type = ModelField1<ModelFieldTagType,FieldType>;
    using tag_type = ModelFieldTagType;
    using field_type = FieldType;
    using update_function_type = std::function<void(field_type &)>;

    ModelField1( tag_type const& thetag, std::string const& prefix, std::string const& name, field_type const& u, std::string const& symbol, std::string const& prefix_symbol, update_function_type const& updateFunction )
        :
        M_tag( thetag ),
        M_prefix( prefix ),
        M_name( name ),
        M_field( u ),
        M_symbol( symbol ),
        M_prefixSymbol( prefix_symbol ),
        M_updateFunction( updateFunction )
        {}

    ModelField1( ModelField1 const& ) = default;
    ModelField1( ModelField1 && ) = default;

    tag_type const& tag() const { return M_tag; }
    std::string const& name() const { return M_name; }
    std::string const& prefix() const { return M_prefix; }
    field_type const& field() const { return M_field; }
    std::string const& symbol() const { return M_symbol; }
    std::string const& prefixSymbol() const { return M_prefixSymbol; }
    SymbolExprUpdateFunction updateFunctionSymbolExpr() const
        {
            if ( M_updateFunction )
                return std::bind( &self_type::applyUpdateFunction, this );
            else
                return SymbolExprUpdateFunction{};
        }

    std::string nameWithPrefix() const { return prefixvm( M_prefix,M_name ); }

    void applyUpdateFunction() const
        {
            if ( M_updateFunction )
                M_updateFunction(const_cast<field_type&>(M_field));
        }
private :
    tag_type M_tag;
    std::string M_prefix;
    std::string M_name;
    field_type M_field;
    std::string M_symbol, M_prefixSymbol;
    update_function_type M_updateFunction;
};

template <size_type Ctx,typename ModelFieldTagType,typename FieldType>
class ModelField : public std::vector<ModelField1<ModelFieldTagType,FieldType> >
{
    static constexpr uint16_type nComponents1 = Feel::remove_shared_ptr_type<FieldType>::nComponents1;
    static constexpr uint16_type nComponents2 = Feel::remove_shared_ptr_type<FieldType>::nComponents2;
    static constexpr uint16_type nRealDim = Feel::remove_shared_ptr_type<FieldType>::nRealDim;
    using functionspace_type = typename Feel::remove_shared_ptr_type<FieldType>::functionspace_type;

    using model_field1_type = ModelField1<ModelFieldTagType,FieldType>;
  public :
    using tag_type = ModelFieldTagType;
    using type = ModelField<Ctx,ModelFieldTagType,FieldType>; // used by boost::hana and find_if
    using field_type = FieldType;
    using feelpp_tag = ModelFieldFeelppTag;
    using update_function_type = typename model_field1_type::update_function_type;

    ModelField() = default;
#if 1
    ModelField( ModelField const& ) = default;
#else
    ModelField( ModelField const& e ) : std::vector<ModelField1<ModelFieldTagType,FieldType> >( e ) { std::cout << "copy constructor"<< std::endl; }
#endif
    ModelField( ModelField && ) = default;
    ModelField( tag_type const& thetag, std::string const& prefix, std::string const& name, field_type const& u, std::string const& symbol = "", std::string const& prefix_symbol = "", update_function_type const& updateFunction = update_function_type{} )
    {
        this->add( thetag, prefix, name, u, symbol, prefix_symbol, updateFunction );
    }

    ModelField& operator=( ModelField&& ) = default;

    void add( tag_type const& thetag, std::string const& prefix, std::string const& name, field_type const& u, std::string const& symbol = "", std::string const& prefix_symbol = "", update_function_type const& updateFunction = update_function_type{} )
    {
        std::string const& symbolNameUsed = symbol.empty()? name : symbol;
        if constexpr ( is_shared_ptr_v<field_type> )
            {
                if ( u )
                    this->push_back( model_field1_type( thetag,prefix, name, u, symbolNameUsed, prefix_symbol, updateFunction ) );
            }
        else
            this->push_back( model_field1_type( thetag,prefix, name, u, symbolNameUsed, prefix_symbol, updateFunction ) );
    }

    auto symbolsExpr() const
    {
        return Feel::vf::symbolsExpr( this->symbolsExpr_ID(),
                                      this->symbolsExpr_MAGNITUDE(),
                                      this->symbolsExpr_GRAD(),
                                      this->symbolsExpr_GRAD_NORMAL()
                                      );

    }

    template <typename SelectorModelFieldType>
    auto trialSymbolsExpr( SelectorModelFieldType const& smf ) const
        {
            return Feel::FeelModels::trialSymbolExprs( this->trialSymbolsExpr_ID( smf ) );
        }
  private :

    auto symbolsExpr_ID() const
    {
        if constexpr ( has_value_v<Ctx,FieldCtx::ID> )
        {
            SymbolExprComponentSuffix secs( nComponents1, nComponents2 );
            using _expr_type = std::decay_t<decltype( idv(this->front().field()) )>;
            symbol_expression_t<_expr_type> se;
            for ( auto const& mfield : *this )
                se.add( prefixvm( mfield.prefixSymbol(),mfield.symbol(),"_" ), idv(mfield.field()), secs, mfield.updateFunctionSymbolExpr() );
            return Feel::vf::symbolsExpr( se );
        }
        else
            return symbols_expression_empty_t{};
    }

    auto symbolsExpr_MAGNITUDE() const
    {
        if constexpr ( has_value_v<Ctx,FieldCtx::MAGNITUDE> )
            {
                SymbolExprComponentSuffix secs( 1, 1 );
                using _expr_type = std::decay_t<decltype( inner(idv(this->front().field()), mpl::int_<InnerProperties::SQRT>() ) )>;
                symbol_expression_t<_expr_type> se;
                for ( auto const& mfield : *this )
                    se.add( prefixvm( mfield.prefixSymbol(), mfield.symbol()+"_magnitude","_" ), inner(idv(mfield.field()),mpl::int_<InnerProperties::SQRT>()) , secs, mfield.updateFunctionSymbolExpr() );
                return Feel::vf::symbolsExpr( se );
            }
        else
            return symbols_expression_empty_t{};
    }

    auto symbolsExpr_GRAD() const
    {
        if constexpr ( has_value_v<Ctx,FieldCtx::GRAD> && nComponents2 == 1 )
        {
            SymbolExprComponentSuffix secs( nComponents1, nRealDim );
            using _expr_type = std::decay_t<decltype( gradv(this->front().field()) )>;
            symbol_expression_t<_expr_type> se;
            for ( auto const& mfield : *this )
                se.add( prefixvm( mfield.prefixSymbol(), "grad_"+mfield.symbol(),"_" ), gradv(mfield.field()), secs, mfield.updateFunctionSymbolExpr() );
            return Feel::vf::symbolsExpr( se );
        }
        else
            return symbols_expression_empty_t{};
    }

    auto symbolsExpr_GRAD_NORMAL() const
    {
        if constexpr ( has_value_v<Ctx,FieldCtx::GRAD_NORMAL> && nComponents2 == 1 )
        {
            SymbolExprComponentSuffix secs( nComponents1, 1 );
            using _expr_type = std::decay_t<decltype(dnv(this->front().field()) )>;
            symbol_expression_t<_expr_type> se;
            for ( auto const& mfield : *this )
                se.add( prefixvm( mfield.prefixSymbol(), "dn_"+mfield.symbol(),"_" ), dnv(mfield.field()), secs, mfield.updateFunctionSymbolExpr() );
            return Feel::vf::symbolsExpr( se );
        }
        else
            return symbols_expression_empty_t{};
    }



    template <typename SelectorModelFieldType>
    auto trialSymbolsExpr_ID( SelectorModelFieldType const& smfs ) const
    {
        if constexpr ( has_value_v<Ctx,FieldCtx::ID> )
        {
            SymbolExprComponentSuffix secs( nComponents1, nComponents2 );
            using _expr_type = std::decay_t<decltype( idt(this->front().field()) )>;
            trial_symbol_expr_t<functionspace_type,_expr_type> tse;
            for ( auto const& mfield : *this )
            {
                for ( auto const& smf : smfs )
                {
                    if ( smf.name() != mfield.name() )
                        continue;
                    if ( smf.tag() != mfield.tag() )
                        continue;

                    size_type startBlockSpaceIndex =3;
                    tse.add( unwrap_ptr( mfield.field() ).functionSpace(), prefixvm( mfield.prefixSymbol(),mfield.symbol(),"_" ), idt(mfield.field()),
                             //secs, mfield.updateFunctionSymbolExpr()
                             startBlockSpaceIndex
                             );
                    break;
                }
            }
            return Feel::FeelModels::trialSymbolExprs( tse );
        }
        else
            return trial_symbols_expr_empty_t{};
    }

};

template <size_type Ctx,typename TagType,typename FieldType>
auto modelField()
{
    return ModelField<Ctx,TagType,FieldType>{};
}
template <size_type Ctx,typename FieldType,typename TagType>
auto modelField( TagType const& )
{
    return ModelField<Ctx,TagType,FieldType>{};
}

template <size_type Ctx,typename TagType,typename FieldType>
auto modelField( TagType const& thetag, std::string const& prefix, std::string const& name, FieldType const& u, std::string const& symbol = "", std::string const& prefix_symbol = "",
                 typename ModelField<Ctx,TagType,FieldType>::update_function_type const& updateFunction = typename ModelField<Ctx,TagType,FieldType>::update_function_type{} )
{
    return ModelField<Ctx,TagType,FieldType>( thetag,prefix,name,u,symbol,prefix_symbol,updateFunction );
}




template <typename FieldTagType>
class SelectorModelField1
{
public :
    using tag_type = FieldTagType;
    SelectorModelField1( tag_type const& tag, std::string const& name )
        :
        M_tag( tag ),
        M_name( name )
        {}
    SelectorModelField1( SelectorModelField1 const& ) = default;
    SelectorModelField1( SelectorModelField1 && ) = default;

    tag_type const& tag() const { return M_tag; }
    std::string const& name() const { return M_name; }

private :
    tag_type M_tag;
    std::string M_name;
};

struct SelectorModelFieldFeelppTag {};
struct SelectorModelFieldsFeelppTag {};

template <typename FieldTagType>
class SelectorModelField : public std::vector<SelectorModelField1<FieldTagType>>
{
    using super_type = std::vector<SelectorModelField1<FieldTagType>>;
    using selector_model_field1_type = SelectorModelField1<FieldTagType>;
public :
    using feelpp_tag = SelectorModelFieldFeelppTag;
    using tag_type = typename selector_model_field1_type::tag_type;

    SelectorModelField() = default;
    SelectorModelField( tag_type const& tag, std::string const& name ) : super_type( 1, selector_model_field1_type(tag,name) ) {}
    SelectorModelField( SelectorModelField const& ) = default;
    SelectorModelField( SelectorModelField && ) = default;
};

template <typename FieldTagType>
SelectorModelField<FieldTagType>
selectorModelField( FieldTagType const& tag, std::string const& name )
{
    return SelectorModelField<FieldTagType>( tag, name );
}

template <typename TupleSelectorModelFieldType>
class SelectorModelFields
{
public :
    using feelpp_tag = SelectorModelFieldsFeelppTag;
    using tuple_type = TupleSelectorModelFieldType;

    SelectorModelFields() = default;
    SelectorModelFields( SelectorModelFields const& ) = default;
    SelectorModelFields( SelectorModelFields && ) = default;
    SelectorModelFields& operator=( SelectorModelFields&& ) = default;

    explicit SelectorModelFields( tuple_type const& tmf )
        :
        M_tuple( tmf )
        {}

    explicit SelectorModelFields( tuple_type && tmf )
        :
        M_tuple( tmf )
        {}

    tuple_type const& tuple() const { return M_tuple; }
private :
    tuple_type M_tuple;
};

template<typename... SelectorModelFieldsType>
struct SelectorModelFieldsTraits
{
    static constexpr auto callApply = [](const auto& ...args) { return Feel::detail::AdvancedConcatOfTupleContainerType<SelectorModelFieldsFeelppTag,SelectorModelFieldFeelppTag>::template apply( args... ); };
    using tuple_type = std::decay_t<decltype( hana::unpack( hana::tuple<SelectorModelFieldsType...>{},  callApply ) )>;
    using type = SelectorModelFields<tuple_type>;
};

template<typename... SelectorModelFieldsType>
using selector_model_fields_t = typename SelectorModelFieldsTraits<SelectorModelFieldsType...>::type;

template<typename... SelectorModelFieldsType>
selector_model_fields_t<SelectorModelFieldsType...>
selectorModelFields( const SelectorModelFieldsType&... smf )
{
    return selector_model_fields_t<SelectorModelFieldsType...>( Feel::detail::AdvancedConcatOfTupleContainerType<SelectorModelFieldsFeelppTag,SelectorModelFieldFeelppTag>::template apply( smf... ) );
}






struct ModelFieldsFeelppTag {};

template <typename TagType>
struct ModelFieldsFindTag
{
    template <typename T>
    struct is_same_tag : public std::integral_constant<bool, std::is_same_v<typename T::tag_type, TagType>> {};

    template <typename C>
    static auto /*const&*/ find( C const& c ) { return hana::find_if( c,  hana::trait< is_same_tag > ); }
};

template <typename TupleFieldType>
class ModelFields
{
public :
    using feelpp_tag = ModelFieldsFeelppTag;
    using tuple_type = TupleFieldType;

    ModelFields() = default;
    ModelFields( ModelFields const& ) = default;
    ModelFields( ModelFields && ) = default;
    ModelFields& operator=( ModelFields&& ) = default;

    explicit ModelFields( tuple_type const& tmf )
        :
        M_tuple( tmf )
        {}

    explicit ModelFields( tuple_type && tmf )
        :
        M_tuple( tmf )
        {}

    auto symbolsExpr() const
        {
            return this->symbolsExprImpl<0>( symbols_expression_empty_t{}, M_tuple );
        }

    template <typename TupleSelectorModelFieldType>
    auto trialSymbolsExpr( SelectorModelFields<TupleSelectorModelFieldType> const& smf ) const
        {
            return this->trialSymbolsExprImpl<0,0>( trial_symbols_expr_empty_t{}, smf.tuple() );
        }

    template <typename TagType>
    auto const&
    field( TagType const& thetag, std::string const& name ) const
        {
            // found the field type related to TagType
            using findFieldT_opt = std::decay_t<decltype( ModelFieldsFindTag<TagType>::find( hana::to_tuple( M_tuple ) ) )>;
            static_assert( !decltype(hana::is_nothing( findFieldT_opt{} ))::value, "tag not found" );
            using found_field_type = typename std::decay_t<decltype( findFieldT_opt{}.value() )>::type::field_type;

            // get the field associated to the name (warning can have several same tag, need to loop over all)
            const found_field_type * dummyRet = nullptr;
            return this->fieldImpl<found_field_type,TagType,0>( thetag,name,dummyRet );
        }

    tuple_type const& tuple() const { return M_tuple; }

private :
    tuple_type M_tuple;

private :

    static const int nModelField = std::decay_t<decltype(hana::size(tuple_type{}))>::value;

    template<int Index,typename ResType,typename TheTupleModelFieldType>
    static auto symbolsExprImpl( ResType && res, TheTupleModelFieldType const& t )
        {
            if constexpr ( Index < nModelField )
                return symbolsExprImpl<Index+1>( Feel::vf::symbolsExpr( std::forward<ResType>( res ), hana::at( t, hana::int_c<Index> ).symbolsExpr() ), t );
            else
                return std::move( res );
        }

    template<int Index,int Index2,typename ResType,typename TupleSelectorModelFieldType>
    auto trialSymbolsExprImpl( ResType && res, TupleSelectorModelFieldType const& t ) const
        {
            if constexpr ( Index < nModelField )
            {
                static const int nSelectorModelField = std::decay_t<decltype(hana::size( TupleSelectorModelFieldType/*tuple_type*/{} ))>::value;
                if constexpr ( Index2 < nSelectorModelField )
                {
                    using CurrentModelFieldType = typename std::decay_t<decltype( hana::at( M_tuple, hana::int_c<Index> ) )>;
                    using CurrentTagType = typename CurrentModelFieldType::tag_type;
                    using SelectorTagType = typename std::decay_t<decltype( hana::at( t, hana::int_c<Index2> ) )>::tag_type;
                    // found the field type related to TagType
                    //using findFieldT_opt = std::decay_t<decltype( ModelFieldsFindTag<TagType>::find( hana::to_tuple( M_tuple ) ) )>;
                    if constexpr ( std::is_same_v<SelectorTagType,CurrentTagType> )
                    {
                        return this->trialSymbolsExprImpl<Index,Index2+1>( Feel::FeelModels::trialSymbolExprs( std::forward<ResType>( res ),
                                                                                                               hana::at( M_tuple, hana::int_c<Index> ).trialSymbolsExpr( hana::at( t, hana::int_c<Index2> ) )
                                                                                                               ), t );
#if 0
                        using trial_symbols_expr_to_insert_type = std::decay_t<decltype(CurrentModelFieldType{}.trialSymbolsExpr(""))>;
                        //static_assert( !decltype(hana::is_nothing( findFieldT_opt{} ))::value, "tag not found" );
                        //using found_field_type = typename std::decay_t<decltype( findFieldT_opt{}.value() )>::type::field_type;
                        trial_symbols_expr_to_insert_type tse_to_insert;
                        for ( auto const& mfield : hana::at( M_tuple, hana::int_c<Index> ) )
                        {
                            for ( auto const& smf : hana::at( t, hana::int_c<Index2> ) )
                            {
                                if ( smf.name() != mfield.name() )
                                    continue;
                                if ( smf.tag() != mfield.tag() )
                                    continue;

                                tse_to_insert = Feel::FeelModels::trialSymbolExprs( tse_to_insert, mfield.trialSymbolsExpr(
                            }
                        }
                        return this->trialSymbolsExprImpl<Index,Index2+1>( Feel::FeelModels::trialSymbolExprs( std::forward<ResType>( res ), tse_to_insert ), t );
#endif
                    }
                    else
                        return this->trialSymbolsExprImpl<Index,Index2+1>( std::forward<ResType>( res ), t );
                }
                else
                {
                    return this->trialSymbolsExprImpl<Index+1,0>( std::forward<ResType>( res ), t );
                }
            }
            else
                return std::move( res );
        }


    template <typename FieldType,typename TagType,int Index>
    auto const&
    fieldImpl( TagType const& thetag, std::string const& name, const FieldType * dummyRet ) const
        {
            if constexpr ( Index < nModelField )
            {
                if constexpr (
                    std::is_same_v< TagType, typename std::decay_t<decltype(hana::at( M_tuple, hana::int_c<Index> ))>::tag_type > &&
                    std::is_same_v< FieldType, typename std::decay_t<decltype(hana::at( M_tuple, hana::int_c<Index> ))>::field_type > )
                    {
                        //std::cout << "size of mfield : " <<  hana::at( M_tuple, hana::int_c<Index> ).size() << std::endl;
                        for ( auto const& mfield : hana::at( M_tuple, hana::int_c<Index> ) )
                        {
                            if ( name != mfield.name() )
                                continue;
                            if ( thetag != mfield.tag() )
                                continue;
                            return mfield.field();
                        }
                    }
                    return fieldImpl<FieldType,TagType,Index+1>( thetag, name, dummyRet );
            }
            else
            {
                CHECK( false ) << "shouldn't go here, something is wrong";
                return *dummyRet;
            }
        }

};


template<typename... MFieldsType>
struct ModelFieldsTraits
{
    static constexpr auto callApply = [](const auto& ...mfields) { return Feel::detail::AdvancedConcatOfTupleContainerType<ModelFieldsFeelppTag,ModelFieldFeelppTag>::template apply( mfields... ); };
    using tuple_type = std::decay_t<decltype( hana::unpack( hana::tuple<MFieldsType...>{},  callApply ) )>;
    using type = ModelFields<tuple_type>;
};

template<typename... MFieldsType>
using model_fields_t = typename ModelFieldsTraits<MFieldsType...>::type;

using model_fields_empty_t = ModelFields<hana::tuple<>>;


template<typename... MFieldsType>
model_fields_t<MFieldsType...>
modelFields( const MFieldsType&... mfields )
{
    return model_fields_t<MFieldsType...>( Feel::detail::AdvancedConcatOfTupleContainerType<ModelFieldsFeelppTag,ModelFieldFeelppTag>::template apply( mfields... ) );
}




} // namespace FeelModels

} // namespace Feel

#endif
