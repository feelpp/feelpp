
#include "fsiconfig.h"

#include <feel/feelmodels2/fsi/fsi.cpp>
#include <feel/feelmodels2/fsi/interpolationfsi.cpp>
#include <feel/feelmodels2/fsi/aitkenrelaxationfsi.cpp>

namespace Feel {
namespace FeelModels {

template class FSI< FLUIDMECHANICS_CLASS_TYPE,SOLIDMECHANICS_CLASS_TYPE >;
template class InterpolationFSI< FLUIDMECHANICS_CLASS_TYPE,SOLIDMECHANICS_CLASS_TYPE >;
template class AitkenRelaxationFSI< SOLIDMECHANICS_CLASS_TYPE >;

}
}

