
#include "multifluidconfig.h"

#include <feel/feelmodels/multifluid/multifluid.cpp>

namespace Feel {
namespace FeelModels {

template class MultiFluid< FLUIDMECHANICS_CLASS_TYPE, LEVELSET_CLASS_TYPE >;

}
}

