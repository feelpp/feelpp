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
class ModelMeasuresQuantity
{
public:
    using self_type = ModelMeasuresQuantity< QuantityType >;
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
    static_assert( std::is_convertible_v<value_type, double> 
            || std::is_convertible_v<value_type, std::vector<double>>,
            "Measure quantity must provide a double-convertible value or a vector of doubles" );

public:
    ModelMeasuresQuantity() = default;
    ModelMeasuresQuantity( ModelMeasuresQuantity const& ) = default;
    ModelMeasuresQuantity( ModelMeasuresQuantity && ) = default;
    ModelMeasuresQuantity( std::string const& prefix, std::string const& name, quantity_type const& q ) :
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

} // namespace FeelModels
} // namespace Feel

#endif
