#ifndef FEELPP_MIXEDPOISSON_ENUM_HPP
#define FEELPP_MIXEDPOISSON_ENUM_HPP

#include <map>

namespace Feel {

namespace FeelModels {

enum MixedPoissonPhysics{
    None = 0,
    Electric,
    Heat,
    Elasticity
};

extern std::map<MixedPoissonPhysics,std::map<std::string,std::string> > MixedPoissonPhysicsMap;

}
}

#endif
