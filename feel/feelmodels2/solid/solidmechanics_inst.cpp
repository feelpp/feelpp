

#include <feel/feelmodels2/solid/solidmechanics.cpp>

namespace Feel {
namespace FeelModels {

#include "solidmecconfig.h"
template class SolidMechanics< Simplex<SOLIDMECHANICS_DIM,SOLIDMECHANICS_ORDERGEO,SOLIDMECHANICS_DIM>, SOLIDMECHANICS_ORDER_DISPLACEMENT,SOLIDMECHANICS_USE_CST_DENSITY_COEFFLAME >;

}
}
