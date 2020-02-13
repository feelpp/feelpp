/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_MODELCORE_MODELPDE_H
#define FEELPP_TOOLBOXES_MODELCORE_MODELPDE_H 1

#include <type_traits> // TO ADD
#include <feel/feelcore/feeltypes.hpp>
#include <string>
#include <set>
#include <map>
#include <vector>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
class ModelPDE
{
protected:
public :
    using material_property_shape_dim_type = std::pair<uint16_type,uint16_type>;
    using material_property_description_type = std::tuple<std::string,std::vector<material_property_shape_dim_type>>;

    ModelPDE() = default;
    ModelPDE( ModelPDE const& ) = default;
    ModelPDE( ModelPDE && ) = default;

    std::map<std::string,std::tuple<std::string,std::vector<material_property_shape_dim_type>>> const& materialPropertyDescription() const { return M_materialPropertyDescription; }

protected :
    void addMaterialPropertyDescription( std::string const& propName, std::string const& symbol, std::initializer_list<material_property_shape_dim_type> const& shapes )
        {
            M_materialPropertyDescription[propName] = std::make_tuple( symbol, std::vector< material_property_shape_dim_type>( shapes ) );
        }

private :
    std::map<std::string,material_property_description_type> M_materialPropertyDescription; // name -> (symbol, shapes.. )
};

} // namespace FeelModels
} // namespace Feel

#endif
