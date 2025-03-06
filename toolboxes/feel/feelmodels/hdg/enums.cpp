#include <map>
#include <string>
#include <feel/feelmodels/hdg/enums.hpp>

namespace Feel {

namespace FeelModels {

std::map<MixedPoissonPhysics, std::map<std::string, std::string>> MixedPoissonPhysicsMap = {
    { MixedPoissonPhysics::None, {
        {"potentialK",       "potential"},
        {"fluxK",            "flux"},
        {"keyword",          "poisson"},
        {"potentialSymbol",  "P"},
        {"fluxSymbol",       "F"}
    }},
    { MixedPoissonPhysics::Electric, {
        {"potentialK",       "electric-potential"},
        {"fluxK",            "current-density"},
        {"keyword",          "electric"},
        {"potentialSymbol",  "P"},
        {"fluxSymbol",       "C"}
    }},
    { MixedPoissonPhysics::Heat, {
        {"potentialK",       "temperature"},
        {"fluxK",            "heat-flux"},
        {"keyword",          "heat"},
        {"potentialSymbol",  "T"},
        {"fluxSymbol",       "F"}
    }},
    { MixedPoissonPhysics::Elasticity, {
        {"potentialK",       "displacement"},
        {"fluxK",            "stress"},
        {"keyword",          "elasticity"},
        {"potentialSymbol",  "d"},
        {"fluxSymbol",       "s"}
    }},
    { MixedPoissonPhysics::Concentration, {
        {"potentialK",       "concentration"},
        {"fluxK",            "concentration-flux"},
        {"keyword",          "concentration"},
        {"potentialSymbol",  "c"}, // or "C"
        {"fluxSymbol",       "j"}  // or "F", "J", etc.
    }},
    { MixedPoissonPhysics::Chemoattractant, {
        {"potentialK",       "chemoattractant"},
        {"fluxK",            "chemoattractant-flux"},
        {"keyword",          "chemoattractant"},
        {"potentialSymbol",  "phi"}, // or "Ch"
        {"fluxSymbol",       "j"}  // or "F", etc.
    }}
};

}
}
