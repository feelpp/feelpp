
#include "fsiconfig.h"

#include <feel/feelmodels/fsi/fsi.cpp>
#include <feel/feelmodels/fsi/fsiinterpolation.cpp>
#include <feel/feelmodels/fsi/aitkenrelaxationfsi.cpp>

namespace Feel {
namespace FeelModels {

template class FSI< FLUIDMECHANICS_CLASS_TYPE,SOLIDMECHANICS_CLASS_TYPE >;
template class AitkenRelaxationFSI< SOLIDMECHANICS_CLASS_TYPE >;

}
}

