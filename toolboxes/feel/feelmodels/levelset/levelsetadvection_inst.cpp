
#include "levelsetconfig.h"
#include <feel/feelmodels/levelset/levelsetspacemanager.hpp>

#define ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS \
    template< typename FunctionSpaceType, typename FunctionSpaceAdvectionVelocityType, typename BasisDiffusionCoeffType, typename BasisReactionCoeffType> \
        /**/
#define ADVDIFFREAC_CLASS_TEMPLATE_TYPE \
    AdvDiffReac<FunctionSpaceType, FunctionSpaceAdvectionVelocityType, BasisDiffusionCoeffType, BasisReactionCoeffType> \
        /**/
#include <feel/feelmodels/advection/advection.cpp>

#define LEVELSETSPACEMANAGER_CLASS_TYPE \
    LevelSetSpaceManager< \
        LEVELSET_CONVEX_TYPE, \
        LEVELSET_BASIS_TYPE, \
        LEVELSET_PERIODICITY_TYPE, \
        LEVELSET_BASISPN_TYPE \
        >                    \
    /**/

namespace Feel {
namespace FeelModels {

extern template class LEVELSETSPACEMANAGER_CLASS_TYPE;

// Scalar advection
template class AdvDiffReac<
    typename FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type,
    LEVELSET_FUNCTIONSPACEADVECTIONVELOCITY_TYPE
    >;
// Vectorial advection
template class AdvDiffReac<
    typename FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_vectorial_type,
    LEVELSET_FUNCTIONSPACEADVECTIONVELOCITY_TYPE
    >;

//#if LEVELSET_ORDERPOLY != LEVELSET_VELOCITY_ORDER
//// Scalar iso advection (for HJ redistanciation)
//template class AdvDiffReac<
    //typename FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type,
    //FunctionSpace<
        //Mesh<LEVELSET_CONVEX_TYPE>,
        //bases<Lagrange<LEVELSET_ORDERPOLY, Vectorial, Continuous, LEVELSET_INTERPOLATION_POINTS>>
            //>
    //>;
//#endif
}
}
