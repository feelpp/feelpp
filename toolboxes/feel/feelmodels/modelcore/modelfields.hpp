
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

namespace FieldTag
{
//static uint16_type next_free_tag = 0;
const uint16_type heat_temperature=0;//next_free_tag++;
const uint16_type electric_potential=1;//next_free_tag++;
const uint16_type fluid_velocity=2;//next_free_tag++;
const uint16_type fluid_pressure=3;//next_free_tag++;
}
namespace FieldCtx
{

const size_type ID           = ( 1<<0 );
const size_type MAGNITUDE    = ( 1<<1 );
const size_type GRAD         = ( 1<<2 );
const size_type GRAD_NORMAL  = ( 1<<3 );

}

struct ModelFieldTag {};

template <uint16_type Tag,size_type Ctx,typename FieldType>
class ModelField : public std::vector<std::tuple<std::string,FieldType,std::string,std::string>> // name, field, symbol, prefix in symbol expr
{
    static constexpr uint16_type nComponents1 = Feel::remove_shared_ptr_type<FieldType>::nComponents1;
    static constexpr uint16_type nComponents2 = Feel::remove_shared_ptr_type<FieldType>::nComponents2;
    static constexpr uint16_type nRealDim = Feel::remove_shared_ptr_type<FieldType>::nRealDim;

  public :
    static constexpr uint16_type tag = Tag;
    using type = ModelField<Tag,Ctx,FieldType>; // used by boost::hana and find_if
    using field_type = FieldType;
    using feelpp_tag = ModelFieldTag;

    ModelField() = default;
    ModelField( ModelField const& ) = default;
    ModelField( ModelField && ) = default;
    ModelField( std::string const& name, field_type const& u, std::string const& symbol = "", std::string const& prefix_symbol = "" )
    {
        if constexpr ( is_shared_ptr_v<field_type> )
            {
                if ( u )
                    this->push_back( std::make_tuple(name,u,symbol.empty()? name : symbol, prefix_symbol) );
            }
        else
            this->push_back( std::make_tuple(name,u,symbol.empty()? name : symbol, prefix_symbol) );
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
            using _expr_type = std::decay_t<decltype(idv( std::get<1>( this->front() ) ) )>;
            std::vector<std::tuple<std::string,_expr_type,SymbolExprComponentSuffix>> _symbsExpr;
            for ( auto const& [name,u,symbol,prefix] : *this )
                _symbsExpr.push_back( std::make_tuple( prefixvm( prefix,symbol,"_" ), idv(u), SymbolExprComponentSuffix( nComponents1, nComponents2, true ) ) );
            return Feel::vf::symbolsExpr( symbolExpr( _symbsExpr ) );
        }
        else
            return symbols_expression_empty_t{};
    }

    auto symbolsExpr_MAGNITUDE() const
    {
        if constexpr ( has_value_v<Ctx,FieldCtx::MAGNITUDE> )
            {
                SymbolExprComponentSuffix secs( 1, 1, true );
                using _expr_type = std::decay_t<decltype(  inner(idv( std::get<1>( this->front() ) ), mpl::int_<InnerProperties::SQRT>() ) )>;
                std::vector<std::tuple<std::string,_expr_type,SymbolExprComponentSuffix>> _symbsExpr;
                for ( auto const& [name,u,symbol,prefix] : *this )
                    _symbsExpr.push_back( std::make_tuple( prefixvm( prefix, symbol+"_magnitude","_" ), inner(idv(u),mpl::int_<InnerProperties::SQRT>()) , secs ) );
                return Feel::vf::symbolsExpr( symbolExpr( _symbsExpr ) );
            }
        else
            return symbols_expression_empty_t{};
    }

    auto symbolsExpr_GRAD() const
    {
        if constexpr ( has_value_v<Ctx,FieldCtx::GRAD> && nComponents2 == 1 )
        {
            SymbolExprComponentSuffix secs( nComponents1, nRealDim, true );
            using _expr_type = std::decay_t<decltype(gradv( std::get<1>( this->front() ) ) )>;
            std::vector<std::tuple<std::string,_expr_type,SymbolExprComponentSuffix>> _symbsExpr;
            for ( auto const& [name,u,symbol,prefix] : *this )
                _symbsExpr.push_back( std::make_tuple( prefixvm( prefix, "grad_"+symbol,"_" ), gradv(u), secs ) );
            return Feel::vf::symbolsExpr( symbolExpr( _symbsExpr ) );
        }
        else
            return symbols_expression_empty_t{};
    }

    auto symbolsExpr_GRAD_NORMAL() const
    {
        if constexpr ( has_value_v<Ctx,FieldCtx::GRAD_NORMAL> && nComponents2 == 1 )
        {
            SymbolExprComponentSuffix secs( nComponents1, 1, true );
            using _expr_type = std::decay_t<decltype(dnv( std::get<1>( this->front() ) ) )>;
            std::vector<std::tuple<std::string,_expr_type,SymbolExprComponentSuffix>> _symbsExpr;
            for ( auto const& [name,u,symbol,prefix] : *this )
                _symbsExpr.push_back( std::make_tuple( prefixvm( prefix, "dn_"+symbol,"_" ), dnv(u), secs ) );
            return Feel::vf::symbolsExpr( symbolExpr( _symbsExpr ) );
        }
        else
            return symbols_expression_empty_t{};
    }

};

template <uint16_type Tag,size_type Ctx,typename FieldType>
auto modelField( std::string const& name, FieldType const& u, std::string const& symbol = "", std::string const& prefix_symbol = "" )
{
    return ModelField<Tag,Ctx,FieldType>( name,u,symbol,prefix_symbol );
}


struct ModelFieldsTag {};

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
            if constexpr ( is_a_t<ModelFieldTag, T1 >::value )
                {
                    return applyImpl( applyModelField(t1,std::forward<ResType>(res) ), mfields... );
                }
            else if constexpr ( is_a_t<ModelFieldsTag, T1 >::value )
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


template <uint16_type Tag>
struct ModelFieldsFindTag
{
    template <typename T>
    struct is_same_tag : public std::integral_constant<bool, T::tag == Tag> {};

    template <typename C>
    static auto find(C const& c ) { return hana::find_if( c,  hana::trait< is_same_tag > ); }
};

template <typename TupleFieldType>
class ModelFields
{
public :
    using feelpp_tag = ModelFieldsTag;
    using tuple_type = TupleFieldType;

    ModelFields() = default;

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
            return this->symbolsExprImpl( symbols_expression_empty_t{}, tupleModelField );
        }

    template <uint16_type TheTAG>
    auto
    field( std::string const& name ) const
        {
            auto findFieldT = ModelFieldsFindTag<TheTAG>::find( tupleModelField );
            static_assert( findFieldT != hana::nothing, "tag not found" );
            using found_field_type = typename std::decay_t<decltype( findFieldT.value() )>::field_type;
            for ( auto const& [name2,field,symbol,prefix_symbol] : findFieldT.value() )
                if ( name == name2 )
                    return field;// std::optional<found_field_type>( field );
            CHECK( false ) << "field name not found";
            return std::get<1>(findFieldT.value().front());
            //return std::optional<found_field_type>{};
        }

#if 0
    template <typename FieldType>
    std::optional<FieldType>
    field( std::string const& name ) const
        {
            std::optional<FieldType> res;
            this->fieldImpl<FieldType,0>( name,res );
            return res;
        }
#endif
    tuple_type tupleModelField;

private :
    template<typename ResType >
    static auto symbolsExprImpl( ResType && res, hana::tuple<> const& t )
        {
            return std::move( res );
        }
    template<typename ResType, typename T1, typename... MFieldsType >
    static auto symbolsExprImpl( ResType && res, hana::tuple<T1,MFieldsType...> const& t )
        {
            return symbolsExprImpl( Feel::vf::symbolsExpr( std::forward<ResType>( res ),  hana::at( t, 0_c ).symbolsExpr() ), hana::remove_at( t, 0_c ) );
        }

#if 0
    static const int nModelField = std::decay_t<decltype(hana::size(tuple_type{}))>::value;

    template <typename FieldType,int Index>
    void
    fieldImpl( std::string const& name, std::optional<FieldType> & res ) const
        {
            if ( res )
                return;
            if constexpr ( Index < nModelField )
            {
                if constexpr ( std::is_same_v< FieldType, typename std::decay_t<decltype(hana::at( tupleModelField, hana::int_c<Index> ))>::field_type > )
                    {
                        for ( auto const& [name2,field,symbol,prefix_symbol] : hana::at( tupleModelField, hana::int_c<Index> ) )
                            if ( name == name2 )
                                res = std::optional<FieldType>( field );
                    }
                fieldImpl<FieldType,Index+1>( name, res );
            }
        }
#endif
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
