
#include "levelsetconfig.h"
#include <feel/feelmodels/levelset/levelset.cpp>

#define LEVELSETSPACEMANAGER_CLASS_TYPE \
    LevelSetSpaceManager< \
        LEVELSET_CONVEX_TYPE, \
        LEVELSET_BASIS_TYPE, \
        LEVELSET_PERIODICITY_TYPE, \
        LEVELSET_BASISPN_TYPE \
        >                    \
    /**/

namespace Feel {

// Redistanciation FM
extern template class LevelSetRedistanciationFM< 
    typename Feel::FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type
    >;

// Redistanciation HJ
extern template class LevelSetRedistanciationHJ< 
    typename Feel::FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type
    >;
//extern template class Feel::FeelModels::AdvDiffReac<
    //typename Feel::FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type
    //>;

namespace FeelModels {

// LevelSetBase
extern template class LEVELSETSPACEMANAGER_CLASS_TYPE;
extern template class LevelSetBase<
    LEVELSET_CONVEX_TYPE,
    LEVELSET_BASIS_TYPE,
    LEVELSET_PERIODICITY_TYPE,
    LEVELSET_BASISPN_TYPE
        >;

// Explicit instantiation
LEVELSET_CLASS_INSTANTIATION

}
}
