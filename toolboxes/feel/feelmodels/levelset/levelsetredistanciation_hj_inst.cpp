
#include "levelsetbaseconfig.h"
#include <feel/feelmodels/levelset/levelsetspacemanager.hpp>
#include <feel/feelmodels/levelset/levelsetredistanciation_hj.cpp>

// Scalar advection required for ReinitializerHJ
#include <feel/feelmodels/advection/advection.cpp>

#define LEVELSETSPACEMANAGER_CLASS_TYPE \
    LevelSetSpaceManager< \
        LEVELSETBASE_CONVEX_TYPE, \
        LEVELSETBASE_BASIS_TYPE, \
        LEVELSETBASE_PERIODICITY_TYPE, \
        LEVELSETBASE_BASISPN_TYPE \
        >                    \
    /**/

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
    typename LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type
    >;

}

template class LevelSetRedistanciationHJ< 
    typename FeelModels::LEVELSETSPACEMANAGER_CLASS_TYPE::space_scalar_type
    >;
}
