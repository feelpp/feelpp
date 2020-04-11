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
#if 0
enum class ToolboxTag
{
    heat = 0,
    electric,
    fluid,
    solid
};
#endif
namespace FieldCtx
{

const size_type ID           = ( 1<<0 );
const size_type MAGNITUDE    = ( 1<<1 );
const size_type GRAD         = ( 1<<2 );
const size_type GRAD_NORMAL  = ( 1<<3 );

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
  private :

    auto symbolsExpr_ID() const
    {
        if constexpr ( has_value_v<Ctx,FieldCtx::ID> )
        {
            SymbolExprComponentSuffix secs( nComponents1, nComponents2, true );
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
                SymbolExprComponentSuffix secs( 1, 1, true );
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
            SymbolExprComponentSuffix secs( nComponents1, nRealDim, true );
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
            SymbolExprComponentSuffix secs( nComponents1, 1, true );
            using _expr_type = std::decay_t<decltype(dnv(this->front().field()) )>;
            symbol_expression_t<_expr_type> se;
            for ( auto const& mfield : *this )
                se.add( prefixvm( mfield.prefixSymbol(), "dn_"+mfield.symbol(),"_" ), dnv(mfield.field()), secs, mfield.updateFunctionSymbolExpr() );
            return Feel::vf::symbolsExpr( se );
        }
        else
            return symbols_expression_empty_t{};
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


struct ModelFieldsFeelppTag {};

//! defined type from input args (variadic expression)
struct ModelFieldsTraits
{
    template <typename Tag, typename T>
    struct is_a_t :
        hana::integral_constant<bool, std::is_same<Tag, typename T::feelpp_tag >::value >
    {};

    template<typename... MFieldsType>
    static constexpr auto apply( const MFieldsType&... mfields )
        {
            return applyImpl( hana::tuple<>{}, mfields... );
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
    template<typename ResType, typename T1, typename... MFieldsType >
    static constexpr auto applyImpl2( ResType && res,  hana::tuple<T1,MFieldsType...> const& t )
        {
            return applyImpl2( applyImpl( std::forward<ResType>( res ), hana::at( t, 0_c ) ),  hana::remove_at( t, 0_c ) );
        }

    template < typename T1, typename... MFieldsType >
    static constexpr auto applyModelField( T1 const& t1, hana::tuple<MFieldsType...> && res )
        {
            if constexpr ( hana::find( hana::to_tuple(hana::tuple_t<MFieldsType...> ), hana::type_c<T1>) == hana::nothing )
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
    template<typename ResType, typename T1, typename... MFieldsType2 >
    static constexpr auto applyImpl( ResType && res, T1 const& t1, const MFieldsType2&... mfields )
        {
            if constexpr ( is_a_t<ModelFieldFeelppTag, T1 >::value )
                {
                    return applyImpl( applyModelField(t1,std::forward<ResType>(res) ), mfields... );
                }
            else if constexpr ( is_a_t<ModelFieldsFeelppTag, T1 >::value )
                {
                    if constexpr ( std::decay_t<decltype(hana::size(t1.tupleModelField))>::value == 0 )
                                     return applyImpl( std::forward<ResType>( res ), mfields... );
                        else
                            return applyImpl(  applyImpl2( std::forward<ResType>( res ), t1.tupleModelField ), mfields... );
                }
            else
                return std::move( res );
        }

};


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
        tupleModelField( tmf )
        {}

    explicit ModelFields( tuple_type && tmf )
        :
        tupleModelField( tmf )
        {}

    auto symbolsExpr() const
        {
            return this->symbolsExprImpl<0>( symbols_expression_empty_t{}, tupleModelField );
        }

    template <typename TagType>
    auto const&
    field( TagType const& thetag, std::string const& name ) const
        {
            // found the field type related to TagType
            using findFieldT_opt = std::decay_t<decltype( ModelFieldsFindTag<TagType>::find( hana::to_tuple( tupleModelField ) ) )>;
            static_assert( !decltype(hana::is_nothing( findFieldT_opt{} ))::value, "tag not found" );
            using found_field_type = typename std::decay_t<decltype( findFieldT_opt{}.value() )>::type::field_type;

            // get the field associated to the name (warning can have several same tag, need to loop over all)
            const found_field_type * dummyRet = nullptr;
            return this->fieldImpl<found_field_type,TagType,0>( thetag,name,dummyRet );
        }

    tuple_type tupleModelField;

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

    template <typename FieldType,typename TagType,int Index>
    auto const&
    fieldImpl( TagType const& thetag, std::string const& name, const FieldType * dummyRet ) const
        {
            if constexpr ( Index < nModelField )
            {
                if constexpr (
                    std::is_same_v< TagType, typename std::decay_t<decltype(hana::at( tupleModelField, hana::int_c<Index> ))>::tag_type > &&
                    std::is_same_v< FieldType, typename std::decay_t<decltype(hana::at( tupleModelField, hana::int_c<Index> ))>::field_type > )
                    {
                        //std::cout << "size of mfield : " <<  hana::at( tupleModelField, hana::int_c<Index> ).size() << std::endl;
                        for ( auto const& mfield : hana::at( tupleModelField, hana::int_c<Index> ) )
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
struct ModelFieldsTraits2
{
    static constexpr auto callApply = [](const auto& ...mfields) { return ModelFieldsTraits::template apply( mfields... ); };
    using tuple_type = std::decay_t<decltype( hana::unpack( hana::tuple<MFieldsType...>{},  callApply ) )>;
    using type = ModelFields<tuple_type>;
};

template<typename... MFieldsType>
using model_fields_t = typename ModelFieldsTraits2<MFieldsType...>::type;

using model_fields_empty_t = ModelFields<hana::tuple<>>;


template<typename... MFieldsType>
model_fields_t<MFieldsType...>
modelFields( const MFieldsType&... mfields )
{
    return model_fields_t<MFieldsType...>( ModelFieldsTraits::template apply( mfields... ) );
}

} // namespace FeelModels

} // namespace Feel

#endif
