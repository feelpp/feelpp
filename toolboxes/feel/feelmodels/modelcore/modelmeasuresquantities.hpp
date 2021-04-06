/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_MODELCORE_MODELMEASURESQUANTITIES_H
#define FEELPP_TOOLBOXES_MODELCORE_MODELMEASURESQUANTITIES_H 1

#include <string>
#include <vector>
#include <functional>
#include <type_traits>

#include <feel/feelcore/traits.hpp>

namespace Feel {
namespace FeelModels {


template< typename QuantityType >
class ModelMeasuresQuantity1
{
public:
    using self_type = ModelMeasuresQuantity1< QuantityType >;
    using quantity_type = QuantityType;
    using value_type = typename std::conditional<
        std::is_invocable_v<quantity_type>,
        std::invoke_result<quantity_type>,
        Feel::type_identity<quantity_type>
            >::type::type;
    //static_assert( std::is_convertible_v<value_type, double>
            //|| Feel::is_iterable_of_v<value_type, double>,
            //"Measure quantity must provide a double-convertible value or an iterable of double-convertible values" );
    // For the moment, we restrict iterables to vectors
    static_assert( std::is_convertible_v<value_type, double> ||
                   std::is_convertible_v<value_type, std::vector<double>> ||
                   is_eigen_matrix_v<value_type> ,
                   "Measure quantity must provide a double-convertible value or a vector of doubles or Eigen::Matrix" );

public:
    ModelMeasuresQuantity1() = default;
    ModelMeasuresQuantity1( ModelMeasuresQuantity1 const& ) = default;
    ModelMeasuresQuantity1( ModelMeasuresQuantity1 && ) = default;
    ModelMeasuresQuantity1( std::string const& prefix, std::string const& name, quantity_type const& q ) :
        M_prefix( prefix ),
        M_name( name ),
        M_quantity( q )
    {}

    std::string const& name() const { return M_name; }
    std::string const& prefix() const { return M_prefix; }
    quantity_type const& quantity() const { return M_quantity; }

    std::string nameWithPrefix() const { return prefixvm( M_prefix, M_name ); }

    value_type value() const
    {
        if constexpr ( std::is_invocable_v<quantity_type> )
        {
            return M_quantity();
        }
        else
        {
            return M_quantity;
        }
    }

private:
    std::string M_prefix;
    std::string M_name;
    quantity_type M_quantity;

};

struct ModelMeasuresQuantityFeelppTag {};

template< typename QuantityType >
class ModelMeasuresQuantity : public std::vector<ModelMeasuresQuantity1<QuantityType>>
{
    using super_type = std::vector<ModelMeasuresQuantity1<QuantityType>>;
    using model_measure_quantity1_type = ModelMeasuresQuantity1<QuantityType>;
public :
    using feelpp_tag = ModelMeasuresQuantityFeelppTag;
    using quantity_type = typename model_measure_quantity1_type::quantity_type;

    ModelMeasuresQuantity() = default;
    ModelMeasuresQuantity( std::string const& prefix, std::string const& name, quantity_type const& q ) : super_type( 1, model_measure_quantity1_type(prefix,name,q) ) {}
    ModelMeasuresQuantity( ModelMeasuresQuantity const& ) = default;
    ModelMeasuresQuantity( ModelMeasuresQuantity && ) = default;
    ModelMeasuresQuantity& operator=( ModelMeasuresQuantity&& ) = default;
};

template <typename QuantityType>
ModelMeasuresQuantity<QuantityType>
modelMeasuresQuantity( std::string const& prefix, std::string const& name, QuantityType const& q )
{
    return ModelMeasuresQuantity<QuantityType>( prefix, name, q );
}


struct ModelMeasuresQuantitiesFeelppTag {};

template <typename TupleQuantityType>
class ModelMeasuresQuantities
{
public :
    using self_type = ModelMeasuresQuantities<TupleQuantityType>;
    using feelpp_tag = ModelMeasuresQuantitiesFeelppTag;
    using tuple_type = TupleQuantityType;

    ModelMeasuresQuantities() = default;
    ModelMeasuresQuantities( ModelMeasuresQuantities const& ) = default;
    ModelMeasuresQuantities( ModelMeasuresQuantities && ) = default;
    ModelMeasuresQuantities& operator=( ModelMeasuresQuantities&& ) = default;

    explicit ModelMeasuresQuantities( tuple_type const& tmf )
        :
        M_tuple( tmf )
        {}

    explicit ModelMeasuresQuantities( tuple_type && tmf )
        :
        M_tuple( tmf )
        {}

    tuple_type const& tuple() const { return M_tuple; }
    tuple_type & tuple() { return M_tuple; }

private :
    tuple_type M_tuple;
};


template<typename... ModelMeasuresQuantitiesType>
struct ModelMeasuresQuantitiesTraits
{
    static constexpr auto callApply = [](const auto& ...args) { return Feel::detail::AdvancedConcatOfTupleContainerType<ModelMeasuresQuantitiesFeelppTag,ModelMeasuresQuantityFeelppTag>::template apply( args... ); };
    using tuple_type = std::decay_t<decltype( hana::unpack( hana::tuple<ModelMeasuresQuantitiesType...>{},  callApply ) )>;
    using type = ModelMeasuresQuantities<tuple_type>;
};

template<typename... ModelMeasuresQuantitiesType>
using model_measures_quantities_t = typename ModelMeasuresQuantitiesTraits<ModelMeasuresQuantitiesType...>::type;

using model_measures_quantities_empty_t = ModelMeasuresQuantities<hana::tuple<>>;

template<typename... ModelMeasuresQuantitiesType>
model_measures_quantities_t<ModelMeasuresQuantitiesType...>
modelMeasuresQuantities( const ModelMeasuresQuantitiesType&... smf )
{
    return model_measures_quantities_t<ModelMeasuresQuantitiesType...>( Feel::detail::AdvancedConcatOfTupleContainerType<ModelMeasuresQuantitiesFeelppTag,ModelMeasuresQuantityFeelppTag>::template apply( smf... ) );
}


} // namespace FeelModels
} // namespace Feel

#endif
