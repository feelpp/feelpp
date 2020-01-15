
#include "levelsetbaseconfig.h"
#include <feel/feelmodels/levelset/levelsetbase.hpp>
#include <feel/feelmodels/levelset/levelsetredistanciation_hj.cpp>

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

// Scalar iso advection (for reinitializerHJ)
//template class AdvDiffReac<
    //typename LevelSetSpaceManager<
        //Simplex<LEVELSETBASE_DIM,LEVELSETBASE_ORDERGEO,LEVELSETBASE_DIM>,
        //Lagrange<LEVELSETBASE_ORDERPOLY, Scalar, Continuous, LEVELSETBASE_INTERPOLATION_POINTS>,
        //LEVELSETBASE_PERIODICITY,
        //Lagrange<LEVELSETBASE_PN_ORDERPOLY, Scalar, Continuous, LEVELSETBASE_INTERPOLATION_POINTS>
            //>::space_scalar_type,
    //FunctionSpace<
        //Mesh<Simplex<LEVELSETBASE_DIM,LEVELSETBASE_ORDERGEO,LEVELSETBASE_DIM>>,
        //bases<Lagrange<LEVELSETBASE_ORDERPOLY, Vectorial, Continuous, LEVELSETBASE_INTERPOLATION_POINTS>>,
        //Periodicity<LEVELSETBASE_PERIODICITY>
            //>
    //>;
template class AdvDiffReac<
    typename LEVELSETBASE_CLASS_TYPE::space_levelset_type
    >;

}

template class LevelSetRedistanciationHJ< 
    typename FeelModels::LEVELSETBASE_CLASS_TYPE::space_levelset_type
    >;
}
