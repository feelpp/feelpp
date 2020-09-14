
#include "levelsetbaseconfig.h"
#include <feel/feelmodels/levelset/levelsetspacemanager.cpp>

namespace Feel {
namespace FeelModels {

template class LevelSetSpaceManager<
    LEVELSETBASE_CONVEX_TYPE,
    LEVELSETBASE_BASIS_TYPE,
    LEVELSETBASE_PERIODICITY_TYPE,
    LEVELSETBASE_BASISPN_TYPE
        >;

}
}
