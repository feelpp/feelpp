#include <map>
#include <string>
#include <feel/feelmodels/hdg/enums.hpp>

namespace Feel {

namespace FeelModels {

std::map<MixedPoissonPhysics,std::map<std::string,std::string> > MixedPoissonPhysicsMap = {
    { MixedPoissonPhysics::None, {{"potentialK", "potential"},{"fluxK","flux"},{"keyword","poisson"},{"potentialSymbol","P"},{"fluxSymbol","F"}} },
    { MixedPoissonPhysics::Electric, {{"potentialK", "electric-potential"},{"fluxK","current-density"},{"keyword","electric"},{"potentialSymbol","P"},{"fluxSymbol","C"}} },
    { MixedPoissonPhysics::Heat, {{"potentialK", "temperature"},{"fluxK","heat-flux"},{"keyword","heat"},{"potentialSymbol","T"},{"fluxSymbol","F"}} }
};

}
}
