/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

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

#include <feel/feelmodels/modelcore/trialssymbolsexpr.hpp>
#include <feel/feelmodels/modelcore/cachedmodelfield.hpp>

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
const size_type CURL         = ( 1<<4 );
const size_type CURL_MAGNITUDE = ( 1<<5 );
const size_type DIV          = ( 1<<6 );

const size_type FULL = FieldCtx::ID|FieldCtx::MAGNITUDE|
                       FieldCtx::GRAD|FieldCtx::GRAD_NORMAL|
                       FieldCtx::CURL|FieldCtx::CURL_MAGNITUDE|
                       FieldCtx::DIV;
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
    using update_function_type = std::function<field_type const&()>;

    using raw_field_type = raw_field_t<field_type>;
    static constexpr uint16_type nComponents1 = raw_field_type::nComponents1;
    static constexpr uint16_type nComponents2 = raw_field_type::nComponents2;
    static constexpr uint16_type nRealDim = raw_field_type::nRealDim;
    using functionspace_type = typename raw_field_type::functionspace_type;

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
    auto field() const 
    {
        if constexpr ( is_cached_model_field_v<field_type> ) // M_field is a CachedModelField
            return M_field.fieldPtr( false ); // we return the field ptr (not updated)
        else
            return M_field; // otherwise we return directly the field
    }
    std::string const& symbol() const { return M_symbol; }
    std::string const& prefixSymbol() const { return M_prefixSymbol; }
    SymbolExprUpdateFunction updateFunctionSymbolExpr() const
    {
        if constexpr ( is_cached_model_field_v<field_type> ) // M_field is a CachedModelField
        {
            return std::bind( &self_type::applyUpdateFunction, this );
        }
        else
        {
            if ( M_updateFunction )
                return std::bind( &self_type::applyUpdateFunction, this );
            else
                return SymbolExprUpdateFunction{};
        }
    }

    std::string nameWithPrefix() const { return prefixvm( M_prefix,M_name ); }

    void applyUpdateFunction() const
    {
        if constexpr ( is_cached_model_field_v<field_type> ) // M_field is a CachedModelField
        {
            const_cast<self_type*>(this)->M_field.update();
        }
        else
        {
            if ( M_updateFunction )
                const_cast<self_type*>(this)->M_field = M_updateFunction();
        }
    }


    template <size_type Ctx>
    auto symbolExpr1( std::enable_if_t< Ctx == FieldCtx::ID >* = nullptr ) const
        {
            SymbolExprComponentSuffix secs( nComponents1, nComponents2 );
            using _expr_type = std::decay_t<decltype( idv(this->field()) )>;
            return typename symbol_expression_t<_expr_type>::symbolexpr1_type( prefixvm( this->prefixSymbol(),this->symbol(),"_" ),
                                                                               idv(this->field()),
                                                                               secs, this->updateFunctionSymbolExpr() );
        }
    template <size_type Ctx>
    auto trialSymbolExpr1( size_type blockSpaceIndex, std::enable_if_t< Ctx == FieldCtx::ID >* = nullptr ) const
        {
            using _expr_type = std::decay_t<decltype( idt(this->field()) )>;
            using trial_symbol_expr1_type = typename trial_symbol_expr_t<functionspace_type,_expr_type>::trial_symbol_expr1_type;
            static_assert( trial_symbol_expr1_type::expr_shape_type::M == nComponents1 && trial_symbol_expr1_type::expr_shape_type::N == nComponents2, "invalid shape");
            std::string symbolNameBase = prefixvm( this->prefixSymbol(),this->symbol(),"_" );
            auto theSpace = unwrap_ptr( this->field() ).functionSpace();
            return trial_symbol_expr1_type( theSpace,symbolNameBase,idt( this->field() ), blockSpaceIndex );
        }
    template <size_type Ctx>
    auto symbolExpr1( std::enable_if_t< Ctx == FieldCtx::MAGNITUDE >* = nullptr ) const
        {
            SymbolExprComponentSuffix secs( 1, 1 );
            using _expr_type = std::decay_t<decltype( norm2(idv(this->field())) )>;
            return typename symbol_expression_t<_expr_type>::symbolexpr1_type( prefixvm( this->prefixSymbol(), this->symbol()+"_magnitude","_" ),
                                                                               norm2(idv(this->field())),
                                                                               secs, this->updateFunctionSymbolExpr() );
        }
    template <size_type Ctx>
    auto symbolExpr1( std::enable_if_t< Ctx == FieldCtx::GRAD >* = nullptr ) const
        {
            SymbolExprComponentSuffix secs( nComponents1, nRealDim );
            using _expr_type = std::decay_t<decltype( gradv(this->field()) )>;
            return typename symbol_expression_t<_expr_type>::symbolexpr1_type( prefixvm( this->prefixSymbol(), "grad_"+this->symbol(),"_" ),
                                                                               gradv(this->field()),
                                                                               secs, this->updateFunctionSymbolExpr() );
        }
    template <size_type Ctx>
    auto trialSymbolExpr1( size_type blockSpaceIndex, std::enable_if_t< Ctx == FieldCtx::GRAD >* = nullptr ) const
        {
            using _expr_type = std::decay_t<decltype( gradt(this->field()) )>;
            using trial_symbol_expr1_type = typename trial_symbol_expr_t<functionspace_type,_expr_type>::trial_symbol_expr1_type;
            static_assert( trial_symbol_expr1_type::expr_shape_type::M == nComponents1 && trial_symbol_expr1_type::expr_shape_type::N == nRealDim, "invalid shape");
            std::string symbolNameBase = prefixvm( this->prefixSymbol(), "grad_"+this->symbol(),"_" );
            auto theSpace = unwrap_ptr( this->field() ).functionSpace();
            return trial_symbol_expr1_type( theSpace,symbolNameBase,gradt( this->field() ), blockSpaceIndex );
        }

    template <size_type Ctx>
    auto symbolExpr1( std::enable_if_t< Ctx == FieldCtx::GRAD_NORMAL >* = nullptr ) const
        {
            static_assert( nComponents2 == 1, "not support tensor2 shape" );
            SymbolExprComponentSuffix secs( nComponents1, 1 );
            using _expr_type = std::decay_t<decltype( dnv(this->field()) )>;
            return typename symbol_expression_t<_expr_type>::symbolexpr1_type( prefixvm( this->prefixSymbol(), "dn_"+this->symbol(),"_" ),
                                                                               dnv(this->field()),
                                                                               secs, this->updateFunctionSymbolExpr() );
        }
    template <size_type Ctx>
    auto symbolExpr1( std::enable_if_t< Ctx == FieldCtx::CURL >* = nullptr ) const
        {
            static_assert( nComponents1 > 1 && nComponents2 == 1, "only for vectorial shape" );
            SymbolExprComponentSuffix secs( nRealDim==3 ? 3 : 1, 1 );
            using _expr_type = std::decay_t<decltype( curlv(this->field()) )>;
            return typename symbol_expression_t<_expr_type>::symbolexpr1_type( prefixvm( this->prefixSymbol(), "curl_"+this->symbol(),"_" ),
                                                                               curlv(this->field()),
                                                                               secs, this->updateFunctionSymbolExpr() );
        }
    template <size_type Ctx>
    auto trialSymbolExpr1( size_type blockSpaceIndex, std::enable_if_t< Ctx == FieldCtx::CURL >* = nullptr ) const
        {
            using _expr_type = std::decay_t<decltype( curlt(this->field()) )>;
            using trial_symbol_expr1_type = typename trial_symbol_expr_t<functionspace_type,_expr_type>::trial_symbol_expr1_type;
            static_assert( trial_symbol_expr1_type::expr_shape_type::M == (nRealDim==3 ? 3 : 1) && trial_symbol_expr1_type::expr_shape_type::N == 1, "invalid shape");
            std::string symbolNameBase = prefixvm( this->prefixSymbol(), "curl_"+this->symbol(),"_" );
            auto theSpace = unwrap_ptr( this->field() ).functionSpace();
            return trial_symbol_expr1_type( theSpace,symbolNameBase,curlt( this->field() ), blockSpaceIndex );
        }

    template <size_type Ctx>
    auto symbolExpr1( std::enable_if_t< Ctx == FieldCtx::CURL_MAGNITUDE >* = nullptr ) const
        {
            static_assert( nComponents1 > 1 && nComponents2 == 1, "only for vectorial shape" );
            SymbolExprComponentSuffix secs( 1, 1 );
            using _expr_type = std::decay_t<decltype( norm2(curlv(this->field())) )>;
            return typename symbol_expression_t<_expr_type>::symbolexpr1_type( prefixvm( this->prefixSymbol(), "curl_"+this->symbol()+"_magnitude","_" ),
                                                                               norm2(curlv(this->field())),
                                                                               secs, this->updateFunctionSymbolExpr() );
        }
    template <size_type Ctx>
    auto symbolExpr1( std::enable_if_t< Ctx == FieldCtx::DIV >* = nullptr ) const
        {
            static_assert( nComponents1 > 1 /*&& nComponents2 == 1*/, "only for vectorial/tensorial shape" );
            SymbolExprComponentSuffix secs( nComponents2==1 ? 1 : nComponents1, 1 );
            using _expr_type = std::decay_t<decltype( divv(this->field()) )>;
            return typename symbol_expression_t<_expr_type>::symbolexpr1_type( prefixvm( this->prefixSymbol(), "div_"+this->symbol(),"_" ),
                                                                               divv(this->field()),
                                                                               secs, this->updateFunctionSymbolExpr() );
        }
    template <size_type Ctx>
    auto trialSymbolExpr1( size_type blockSpaceIndex, std::enable_if_t< Ctx == FieldCtx::DIV >* = nullptr ) const
        {
            using _expr_type = std::decay_t<decltype( divt(this->field()) )>;
            using trial_symbol_expr1_type = typename trial_symbol_expr_t<functionspace_type,_expr_type>::trial_symbol_expr1_type;
            static_assert( trial_symbol_expr1_type::expr_shape_type::M == (nComponents2==1 ? 1 : nComponents1) && trial_symbol_expr1_type::expr_shape_type::N == 1, "invalid shape");
            std::string symbolNameBase = prefixvm( this->prefixSymbol(), "div_"+this->symbol(),"_" );
            auto theSpace = unwrap_ptr( this->field() ).functionSpace();
            return trial_symbol_expr1_type( theSpace,symbolNameBase,divt( this->field() ), blockSpaceIndex );
        }


    template <size_type Ctx>
    void updateInformationObject( nl::json & p ) const
        {
            p.emplace( "name", this->name() );
            if ( !M_prefix.empty() )
                p.emplace( "prefix", this->prefix() );
            p.emplace( "base symbol", this->symbol() );
            p.emplace( "prefix symbol", this->prefixSymbol() );
            if constexpr ( is_cached_model_field_v<field_type> ) // M_field is a CachedModelField
                p.emplace( "function space", "unknown (cached model field)" );
            else
                p.emplace( "function space", unwrap_ptr( M_field ).functionSpace()->journalSection().to_string() );

            nl::json::array_t jaSE;

            hana::for_each(hana::make_tuple( std::make_tuple( std::integral_constant<size_type,FieldCtx::ID>{}, "idv(.)",  "eval of "+ M_name ),
                                             std::make_tuple( std::integral_constant<size_type,FieldCtx::MAGNITUDE>{}, "norm2(.)",  "norm2 of "+ M_name ),
                                             std::make_tuple( std::integral_constant<size_type,FieldCtx::GRAD>{}, "gradv(.)",  "grad of "+ M_name ),
                                             std::make_tuple( std::integral_constant<size_type,FieldCtx::GRAD_NORMAL>{}, "dnv(.)",  "normal derivative of "+ M_name ),
                                             std::make_tuple( std::integral_constant<size_type,FieldCtx::CURL>{}, "curlv(.)",  "curl of "+ M_name ),
                                             std::make_tuple( std::integral_constant<size_type,FieldCtx::CURL_MAGNITUDE>{}, "norm2(curlv(.))",  "norm2 of curl of "+ M_name ),
                                             std::make_tuple( std::integral_constant<size_type,FieldCtx::DIV>{}, "divv(.)",  "div of "+ M_name )
                                             ), [this,&jaSE](auto const& x) {
                    constexpr size_type thectx = std::decay_t<decltype( std::get<0>( x ) )>::value;
                    if constexpr ( has_value_v<Ctx,thectx> )
                        {
                            auto se1 = this->symbolExpr1<thectx>();
                            jaSE.push_back( symbolExprInformations( se1.symbol(), std::get<1>(x), se1.componentSuffix(), std::get<2>(x) ) );
                        }
                });
            p.emplace( "SymbolsExpr", jaSE );
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
    using raw_field_type = raw_field_t<FieldType>;
    static constexpr uint16_type nComponents1 = raw_field_type::nComponents1;
    static constexpr uint16_type nComponents2 = raw_field_type::nComponents2;
    static constexpr uint16_type nRealDim = raw_field_type::nRealDim;
    using functionspace_type = typename raw_field_type::functionspace_type;

    using model_field1_type = ModelField1<ModelFieldTagType,FieldType>;

    static constexpr size_type ctx_clean_scalar = nComponents1 == 1? clear_values_v<Ctx,FieldCtx::CURL,FieldCtx::CURL_MAGNITUDE|FieldCtx::DIV> : Ctx;
    static constexpr size_type ctx_clean_vectorial = ctx_clean_scalar;
    static constexpr size_type ctx_clean_tensor2 = nComponents2 > 1? clear_values_v<ctx_clean_vectorial,FieldCtx::CURL,FieldCtx::CURL_MAGNITUDE|FieldCtx::GRAD_NORMAL> : ctx_clean_vectorial;
    static constexpr size_type ctx_clean = ctx_clean_tensor2;

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
        return Feel::vf::symbolsExpr( this->symbolsExprImpl<FieldCtx::ID>(),
                                      this->symbolsExprImpl<FieldCtx::MAGNITUDE>(),
                                      this->symbolsExprImpl<FieldCtx::GRAD>(),
                                      this->symbolsExprImpl<FieldCtx::GRAD_NORMAL>(),
                                      this->symbolsExprImpl<FieldCtx::CURL>(),
                                      this->symbolsExprImpl<FieldCtx::CURL_MAGNITUDE>(),
                                      this->symbolsExprImpl<FieldCtx::DIV>()
                                      );
    }

    template <typename SelectorModelFieldType>
    auto trialSymbolsExpr( SelectorModelFieldType const& smf ) const
        {
            return Feel::FeelModels::trialSymbolsExpr<functionspace_type>( this->trialSymbolsExprImpl<FieldCtx::ID>( smf ),
                                                                           this->trialSymbolsExprImpl<FieldCtx::GRAD>( smf ),
                                                                           this->trialSymbolsExprImpl<FieldCtx::CURL>( smf ),
                                                                           this->trialSymbolsExprImpl<FieldCtx::DIV>( smf )
                                                                           );
        }

    void updateInformationObject( nl::json & p ) const
        {
            for ( auto const& mfield : *this )
            {
                mfield.template updateInformationObject<ctx_clean>( p[mfield.nameWithPrefix()] );
            }
        }

  private :
    template <size_type TheContext>
    auto symbolsExprImpl() const
    {
        if constexpr ( has_value_v<ctx_clean,TheContext> )
        {
            using _expr_type = typename std::decay_t<decltype( this->front().template symbolExpr1<TheContext>() ) >::expr_type;
            symbol_expression_t<_expr_type> se;
            for ( auto const& mfield : *this )
                se.push_back( mfield.template symbolExpr1<TheContext>() );
            return Feel::vf::symbolsExpr( se );
        }
        else
            return symbols_expression_empty_t{};
    }

    template <size_type TheContext,typename SelectorModelFieldType>
    auto trialSymbolsExprImpl( SelectorModelFieldType const& smfs ) const
        {
            if constexpr ( has_value_v<ctx_clean,TheContext> )
            {
                using _expr_type = typename std::decay_t<decltype( this->front().template trialSymbolExpr1<TheContext>(0) ) >::expr_type;
                trial_symbol_expr_t<functionspace_type,_expr_type> tse;
                for ( auto const& mfield : *this )
                {
                    for ( auto const& smf : smfs )
                    {
                        if ( smf.name() != mfield.name() )
                            continue;
                        if ( smf.tag() != mfield.tag() )
                            continue;

                        tse.push_back( mfield.template trialSymbolExpr1<TheContext>( smf.blockSpaceIndex() ) );
                    }
                }
                return Feel::FeelModels::trialSymbolsExpr<functionspace_type>( tse );
            }
            else
                return trial_symbols_expr_empty_t<functionspace_type>{};
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
    SelectorModelField1( tag_type const& tag, std::string const& name, size_type blockSpaceIndex )
        :
        M_tag( tag ),
        M_name( name ),
        M_blockSpaceIndex( blockSpaceIndex )
        {}
    SelectorModelField1( SelectorModelField1 const& ) = default;
    SelectorModelField1( SelectorModelField1 && ) = default;

    tag_type const& tag() const { return M_tag; }
    std::string const& name() const { return M_name; }
    size_type blockSpaceIndex() const { return M_blockSpaceIndex; }

private :
    tag_type M_tag;
    std::string M_name;
    size_type M_blockSpaceIndex;
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
    SelectorModelField( tag_type const& tag, std::string const& name, size_type blockSpaceIndex ) : super_type( 1, selector_model_field1_type(tag,name,blockSpaceIndex) ) {}
    SelectorModelField( SelectorModelField const& ) = default;
    SelectorModelField( SelectorModelField && ) = default;
    SelectorModelField& operator=( SelectorModelField&& ) = default;
};

template <typename FieldTagType>
SelectorModelField<FieldTagType>
selectorModelField( FieldTagType const& tag, std::string const& name, size_type blockSpaceIndex = 0 )
{
    return SelectorModelField<FieldTagType>( tag, name, blockSpaceIndex );
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
    using self_type = ModelFields<TupleFieldType>;
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
            return this->trialSymbolsExprImpl<0,0>( trials_symbols_expr_empty_t{}, smf.tuple() );
        }

    template <typename TagType>
    auto
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

    //! create a new ModelFields object by copying the current object but exclude fields defined in mfieldsToExclude (with same tag and name)
    template <typename AnOtherTupleFieldType>
    self_type exclude( ModelFields<AnOtherTupleFieldType> const& mfieldsToExclude ) const
        {
            auto newTupleField = hana::transform( this->tuple(), [&mfieldsToExclude](auto const& ef)
                                                  {
                                                      using mfield_type = std::decay_t<decltype(ef)>;
                                                      mfield_type mfieldNew;
                                                      for ( auto const& mfield : ef )
                                                      {
                                                          bool hasFoundThisField = false;
                                                          hana::for_each( mfieldsToExclude.tuple(), [&mfield,&hasFoundThisField](auto const& ef2)
                                                                          {
                                                                              using mfield2_type = std::decay_t<decltype(ef2)>;
                                                                              if constexpr ( std::is_same_v< typename mfield_type::tag_type, typename mfield2_type::tag_type > )
                                                                                  {
                                                                                      for ( auto const& mfield2 : ef2 )
                                                                                      {
                                                                                          if ( mfield.name() != mfield2.name() || mfield.tag() != mfield2.tag() )
                                                                                              continue;
                                                                                          hasFoundThisField = true;
                                                                                      }
                                                                                  }
                                                                          });
                                                          if ( hasFoundThisField )
                                                              continue;
                                                          mfieldNew.push_back( mfield );
                                                      }
                                                      return mfieldNew;
                                                  });

            return self_type{ std::move(newTupleField) };
        }

    tuple_type const& tuple() const { return M_tuple; }
    tuple_type & tuple() { return M_tuple; }

    void updateInformationObject( nl::json & p ) const
        {
            hana::for_each( this->tuple(), [&p]( auto const& e ) { e.updateInformationObject( p ); });
        }

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
                    if constexpr ( std::is_same_v<SelectorTagType,CurrentTagType> )
                    {
                        return this->trialSymbolsExprImpl<Index,Index2+1>( Feel::FeelModels::trialsSymbolsExpr( std::forward<ResType>( res ),
                                                                                                                hana::at( M_tuple, hana::int_c<Index> ).trialSymbolsExpr( hana::at( t, hana::int_c<Index2> ) )
                                                                                                                ), t );
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
    auto
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
