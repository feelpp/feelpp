/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/modelcore/modelphysics.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
ModelPhysics<Dim>::ModelPhysics( std::string const& physic )
    :
    M_physic( physic )
{
    if ( M_physic == "heat" )
    {
        M_physics = { "heat" };
    }
    else if ( M_physic == "electric" )
    {
        M_physics = { "electric" };
    }
    else if ( M_physic == "fluid" )
    {
        M_physics = { "fluid" };
    }
    else if ( M_physic == "body" )
    {
        M_physics = { "body" };
    }
    else if ( M_physic == "thermo-electric" )
    {
        M_physics = { "heat", "electric", "thermo-electric" };
        M_mapPhysicsToSubphysics["thermo-electric"] = { "heat", "electric" };
    }
    else if ( M_physic == "aerothermal" )
    {
        M_physics = { "heat", "fluid", "aerothermal" };
        M_mapPhysicsToSubphysics["aerothermal"] = { "heat", "fluid" };
    }


    material_property_shape_dim_type scalarShape = std::make_pair(1,1);
    material_property_shape_dim_type matrixShape = std::make_pair(nDim,nDim);

    this->addMatertialPropertyDescription( "density", "rho", { scalarShape } );
    if ( M_physic == "heat" ||  M_physic == "thermo-electric" || M_physic == "aerothermal" )
    {
        this->addMatertialPropertyDescription( "specific-heat-capacity", "Cp", { scalarShape } );
        this->addMatertialPropertyDescription( "thermal-expansion", "beta", { scalarShape } );
        this->addMatertialPropertyDescription( "thermal-conductivity", "k", { scalarShape,matrixShape } );
    }
    if ( M_physic == "electric" || M_physic == "thermo-electric" )
    {
        this->addMatertialPropertyDescription( "electric-conductivity", "sigma", { scalarShape } );
    }

}

template class ModelPhysics<2>;
template class ModelPhysics<3>;

} // namespace FeelModels
} // namespace Feel
