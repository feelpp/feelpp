
#include "levelsetconfig.h"
#include <feel/feelmodels/levelset/levelsetspacemanager.hpp>

#define ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS \
    template< typename FunctionSpaceType, typename BasisDiffusionCoeffType, typename BasisReactionCoeffType> \
        /**/
#define ADVECTIONBASE_CLASS_TEMPLATE_TYPE \
    AdvectionBase<FunctionSpaceType, BasisDiffusionCoeffType, BasisReactionCoeffType> \
        /**/
#include <feel/feelmodels/advection/advectionbase.cpp>
#define ADVECTION_CLASS_TEMPLATE_DECLARATIONS \
    template< typename FunctionSpaceType, typename BasisDiffusionCoeffType, typename BasisReactionCoeffType> \
        /**/
#define ADVECTION_CLASS_TEMPLATE_TYPE \
    Advection<FunctionSpaceType, BasisDiffusionCoeffType, BasisReactionCoeffType> \
        /**/
#include <feel/feelmodels/advection/advection.cpp>

namespace Feel {
namespace FeelModels {

// Scalar advection
template class AdvectionBase<
    typename LevelSetSpaceManager<
        Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>,
        Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>,
        LEVELSET_PERIODICITY,
        Lagrange<LEVELSET_PN_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>
            >::space_scalar_type
    //FunctionSpace< 
        //Mesh<Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>>,
        //bases<Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>>,
        //Periodicity<LEVELSET_PERIODICITY>
        //>
    >;
template class Advection<
    typename LevelSetSpaceManager<
        Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>,
        Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>,
        LEVELSET_PERIODICITY,
        Lagrange<LEVELSET_PN_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>
            >::space_scalar_type
    //FunctionSpace< 
        //Mesh<Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>>,
        //bases<Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>>,
        //Periodicity<LEVELSET_PERIODICITY>
        //>
    >;
// Vectorial advection
template class AdvectionBase<
    typename LevelSetSpaceManager<
        Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>,
        Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>,
        LEVELSET_PERIODICITY,
        Lagrange<LEVELSET_PN_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>
            >::space_vectorial_type
    >;
template class Advection<
    typename LevelSetSpaceManager<
        Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>,
        Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>,
        LEVELSET_PERIODICITY,
        Lagrange<LEVELSET_PN_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>
            >::space_vectorial_type
    >;


}
}
