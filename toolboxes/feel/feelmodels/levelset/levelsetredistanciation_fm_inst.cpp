
#include "levelsetbaseconfig.h"
#include <feel/feelmodels/levelset/levelsetspacemanager.hpp>
#include <feel/feelmodels/levelset/levelsetredistanciation_fm.cpp>

#define LEVELSETSPACEMANAGER_CLASS_TYPE \
    LevelSetSpaceManager< \
        LEVELSETBASE_CONVEX_TYPE, \
        LEVELSETBASE_BASIS_TYPE, \
        LEVELSETBASE_PERIODICITY_TYPE, \
        LEVELSETBASE_BASISPN_TYPE \
        >                    \
    /**/

namespace Feel {

extern template class FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE;

template class LevelSetRedistanciationFM< 
    typename FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type
    >;

}
