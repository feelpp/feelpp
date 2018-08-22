
#include "levelsetconfig.h"
#include <feel/feelmodels/levelset/levelsetspacemanager.hpp>

#define ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS \
    template< typename FunctionSpaceType, typename FunctionSpaceAdvectionVelocityType, typename BasisDiffusionCoeffType, typename BasisReactionCoeffType> \
        /**/
#define ADVECTIONBASE_CLASS_TEMPLATE_TYPE \
    AdvectionBase<FunctionSpaceType, FunctionSpaceAdvectionVelocityType, BasisDiffusionCoeffType, BasisReactionCoeffType> \
        /**/
#include <feel/feelmodels/advection/advectionbase.cpp>
#define ADVECTION_CLASS_TEMPLATE_DECLARATIONS \
    template< typename FunctionSpaceType, typename FunctionSpaceAdvectionVelocityType, typename BasisDiffusionCoeffType, typename BasisReactionCoeffType> \
        /**/
#define ADVECTION_CLASS_TEMPLATE_TYPE \
    Advection<FunctionSpaceType, FunctionSpaceAdvectionVelocityType, BasisDiffusionCoeffType, BasisReactionCoeffType> \
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
            >::space_scalar_type,
    FunctionSpace<
        Mesh<Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>>,
        Feel::detail::bases<Lagrange<LEVELSET_VELOCITY_ORDER, Vectorial, Continuous, LEVELSET_INTERPOLATION_POINTS>>,
        double,
        Periodicity<NoPeriodicity>,
        mortars<NoMortar>
            >
    >;
template class Advection<
    typename LevelSetSpaceManager<
        Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>,
        Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>,
        LEVELSET_PERIODICITY,
        Lagrange<LEVELSET_PN_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>
            >::space_scalar_type,
    FunctionSpace<
        Mesh<Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>>,
        Feel::detail::bases<Lagrange<LEVELSET_VELOCITY_ORDER, Vectorial, Continuous, LEVELSET_INTERPOLATION_POINTS>>,
        double,
        Periodicity<LEVELSET_PERIODICITY>,
        mortars<NoMortar>
            >
    >;
// Vectorial advection
template class AdvectionBase<
    typename LevelSetSpaceManager<
        Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>,
        Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>,
        LEVELSET_PERIODICITY,
        Lagrange<LEVELSET_PN_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>
            >::space_vectorial_type,
    FunctionSpace<
        Mesh<Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>>,
        Feel::detail::bases<Lagrange<LEVELSET_VELOCITY_ORDER, Vectorial, Continuous, LEVELSET_INTERPOLATION_POINTS>>,
        double,
        Periodicity<LEVELSET_PERIODICITY>,
        mortars<NoMortar>
            >
    >;
template class Advection<
    typename LevelSetSpaceManager<
        Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>,
        Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>,
        LEVELSET_PERIODICITY,
        Lagrange<LEVELSET_PN_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>
            >::space_vectorial_type,
    FunctionSpace<
        Mesh<Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>>,
        Feel::detail::bases<Lagrange<LEVELSET_VELOCITY_ORDER, Vectorial, Continuous, PointSetFekete>>,
        double,
        Periodicity<LEVELSET_PERIODICITY>,
        mortars<NoMortar>
            >
    >;

// Scalar iso advection (for reinitializerHJ)
template class AdvectionBase<
    typename LevelSetSpaceManager<
        Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>,
        Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>,
        LEVELSET_PERIODICITY,
        Lagrange<LEVELSET_PN_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>
            >::space_scalar_type,
    FunctionSpace<
        Mesh<Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>>,
        bases<Lagrange<LEVELSET_ORDERPOLY, Vectorial, Continuous, LEVELSET_INTERPOLATION_POINTS>>,
        Periodicity<NoPeriodicity>
            >
    >;
template class Advection<
    typename LevelSetSpaceManager<
        Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>,
        Lagrange<LEVELSET_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>,
        LEVELSET_PERIODICITY,
        Lagrange<LEVELSET_PN_ORDERPOLY, Scalar, Continuous, LEVELSET_INTERPOLATION_POINTS>
            >::space_scalar_type,
    FunctionSpace<
        Mesh<Simplex<LEVELSET_DIM,LEVELSET_ORDERGEO,LEVELSET_DIM>>,
        bases<Lagrange<LEVELSET_ORDERPOLY, Vectorial, Continuous, LEVELSET_INTERPOLATION_POINTS>>,
        Periodicity<LEVELSET_PERIODICITY>
            >
    >;

}
}
