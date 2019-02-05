
#include "levelsetbaseconfig.h"
#include <feel/feelmodels/levelset/levelsetbase.cpp>

// Scalar advection required for ReinitializerHJ
#define ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS \
    template< typename FunctionSpaceType, typename FunctionSpaceAdvectionVelocityType, typename BasisDiffusionCoeffType, typename BasisReactionCoeffType> \
        /**/
#define ADVDIFFREAC_CLASS_TEMPLATE_TYPE \
    AdvDiffReac<FunctionSpaceType, FunctionSpaceAdvectionVelocityType, BasisDiffusionCoeffType, BasisReactionCoeffType> \
        /**/
#include <feel/feelmodels/advection/advection.cpp>

namespace Feel {
namespace FeelModels {

LEVELSETBASE_CLASS_INSTANTIATION

// Scalar iso advection (for reinitializerHJ)
template class AdvDiffReac<
    typename LevelSetSpaceManager<
        Simplex<LEVELSETBASE_DIM,LEVELSETBASE_ORDERGEO,LEVELSETBASE_DIM>,
        Lagrange<LEVELSETBASE_ORDERPOLY, Scalar, Continuous, LEVELSETBASE_INTERPOLATION_POINTS>,
        LEVELSETBASE_PERIODICITY,
        Lagrange<LEVELSETBASE_PN_ORDERPOLY, Scalar, Continuous, LEVELSETBASE_INTERPOLATION_POINTS>
            >::space_scalar_type,
    FunctionSpace<
        Mesh<Simplex<LEVELSETBASE_DIM,LEVELSETBASE_ORDERGEO,LEVELSETBASE_DIM>>,
        bases<Lagrange<LEVELSETBASE_ORDERPOLY, Vectorial, Continuous, LEVELSETBASE_INTERPOLATION_POINTS>>,
        Periodicity<LEVELSETBASE_PERIODICITY>
            >
    >;

}
}
