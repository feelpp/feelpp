
#include "levelsetbaseconfig.h"
#include <feel/feelmodels/levelset/levelsetbase.cpp>

#define LEVELSETSPACEMANAGER_CLASS_TYPE \
    LevelSetSpaceManager< \
        LEVELSETBASE_CONVEX_TYPE, \
        LEVELSETBASE_BASIS_TYPE, \
        LEVELSETBASE_PERIODICITY_TYPE, \
        LEVELSETBASE_BASISPN_TYPE \
        >                    \
    /**/

namespace Feel {

// Redistanciation FM
extern template class LevelSetRedistanciationFM< 
    typename FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type
    >;

// Redistanciation HJ
extern template class LevelSetRedistanciationHJ< 
    typename FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type
    >;
#if 0
extern template class FeelModels::AdvDiffReac<
    typename FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type
    >;
#endif

namespace FeelModels {

extern template class LEVELSETSPACEMANAGER_CLASS_TYPE;

LEVELSETBASE_CLASS_INSTANTIATION

}
}
