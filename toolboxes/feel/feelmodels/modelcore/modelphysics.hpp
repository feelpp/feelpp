/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_MODELCORE_MODELPHYSICS_H
#define FEELPP_TOOLBOXES_MODELCORE_MODELPHYSICS_H 1

#include <type_traits> // TO ADD
#include <feel/feelcore/feel.hpp>
//#include <feel/feelcore/feeltypes.hpp>
#include <string>
#include <set>
#include <map>
#include <vector>

//#include <feel/feelmodels/modelcore/modelpde.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
class ModelPhysics // : public ModelPDE<Dim>
{
public :
    using material_property_shape_dim_type = std::pair<uint16_type,uint16_type>;
    using material_property_description_type = std::tuple<std::string,std::vector<material_property_shape_dim_type>>;

    //using super_type = ModelPDE<Dim>;
    //using material_property_shape_dim_type = typename super_type::material_property_shape_dim_type;
    static const uint16_type nDim = Dim;
    //public :
    ModelPhysics() = default;
    explicit ModelPhysics( std::string const& physic );
    ModelPhysics( ModelPhysics const& ) = default;
    ModelPhysics( ModelPhysics && ) = default;

    //! return name of physic assiciated to system of pde
    std::string const& physic() const { return M_physic; }
    //! return name of all physic present in system of pde
    std::set<std::string> const& physics() const { return M_physics; }
    //! map between physic to subphysics
    std::map<std::string,std::set<std::string>> const& mapPhysicsToSubphysics() const { return M_mapPhysicsToSubphysics; }
    //! return the material properties description (.i.e. coefficients of pdes)
    std::map<std::string,std::tuple<std::string,std::vector<material_property_shape_dim_type>>> const& materialPropertyDescription() const { return M_materialPropertyDescription; }

protected :
    void addMaterialPropertyDescription( std::string const& propName, std::string const& symbol, std::initializer_list<material_property_shape_dim_type> const& shapes )
        {
            M_materialPropertyDescription[propName] = std::make_tuple( symbol, std::vector< material_property_shape_dim_type>( shapes ) );
        }

    //private :
    std::string M_physic; // the name of the physic
    std::set<std::string> M_physics; // all physics (example thermo-electric include also heat and electric)
    std::map<std::string,std::set<std::string>> M_mapPhysicsToSubphysics;

    std::map<std::string,material_property_description_type> M_materialPropertyDescription; // name -> (symbol, shapes.. )

};

} // namespace FeelModels
} // namespace Feel

#endif
