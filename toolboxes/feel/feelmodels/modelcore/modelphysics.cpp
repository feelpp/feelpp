/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/modelcore/modelphysics.hpp>

namespace Feel
{
namespace FeelModels
{

ModelPhysics::ModelPhysics( std::string const& physic )
    :
    M_physic( physic )
{
    if ( M_physic == "heat" )
    {
        M_physics = { "heat" };
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

}

} // namespace FeelModels
} // namespace Feel
