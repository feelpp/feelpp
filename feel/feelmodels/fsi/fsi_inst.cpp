
#include "fsiconfig.h"

#include <feel/feelmodels/fsi/aitkenrelaxationfsi.cpp>
#include <feel/feelmodels/fsi/fsi.cpp>
#include <feel/feelmodels/fsi/interpolationfsi.cpp>

namespace Feel
{
namespace FeelModels
{

template class FSI<FLUIDMECHANICS_CLASS_TYPE, SOLIDMECHANICS_CLASS_TYPE>;
template class InterpolationFSI<FLUIDMECHANICS_CLASS_TYPE, SOLIDMECHANICS_CLASS_TYPE>;
template class AitkenRelaxationFSI<SOLIDMECHANICS_CLASS_TYPE>;
}
}
