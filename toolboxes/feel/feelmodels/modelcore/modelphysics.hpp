/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_MODELCORE_MODELPHYSICS_H
#define FEELPP_TOOLBOXES_MODELCORE_MODELPHYSICS_H 1

#include <string>
#include <set>
#include <map>

namespace Feel
{
namespace FeelModels
{

class ModelPhysics
{
public :
    ModelPhysics( std::string const& physic );
    ModelPhysics( ModelPhysics const& ) = default;
    ModelPhysics( ModelPhysics && ) = default;

    std::string const& physic() const { return M_physic; }
    std::set<std::string> const& physics() const { return M_physics; }
    std::map<std::string,std::set<std::string>> const& mapPhysicsToSubphysics() const { return M_mapPhysicsToSubphysics; }

private :
    std::string M_physic; // the name of the physica
    std::set<std::string> M_physics; // all physics (example thermo-electric include also heat and electric)
    std::map<std::string,std::set<std::string>> M_mapPhysicsToSubphysics;
};

} // namespace FeelModels
} // namespace Feel

#endif
