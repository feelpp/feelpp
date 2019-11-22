#include <feel/feelmodels/hdg/enums.hpp>

namespace Feel {

namespace FeelModels {

std::map<MixedPoissonPhysics,std::map<std::string,std::string> > MixedPoissonPhysicsMap = {
    { MixedPoissonPhysics::None, {{"potentialK", "potential"},{"fluxK","flux"}} },
    { MixedPoissonPhysics::Electric, {{"potentialK", "electric-potential"},{"fluxK","current-density"}} },
    { MixedPoissonPhysics::Heat, {{"potentialK", "temperature"},{"fluxK","heat-flux"}} }
};

}
}
