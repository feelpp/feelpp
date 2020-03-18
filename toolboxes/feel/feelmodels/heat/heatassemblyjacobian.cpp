/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels/heat/heat.hpp>

namespace Feel {
namespace FeelModels {

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    if ( M_bcDirichlet.empty() ) return;

    this->log("Heat","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto mesh = this->mesh();
    size_type startBlockIndexTemperature = this->startSubBlockSpaceIndex( "temperature" );
    auto u = this->spaceTemperature()->element( U, this->rowStartInVector()+startBlockIndexTemperature );
    auto se = this->symbolsExpr();

    for( auto const& d : M_bcDirichlet )
    {
        auto theExpr = expression(d,se);
        u.on(_range=markedfaces(mesh, M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) ) ),
             _expr=theExpr );
    }

    // update info for synchronization
    this->updateDofEliminationIds( "temperature", data );

    this->log("Heat","updateNewtonInitialGuess","finish" );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateJacobian( data, this->modelContext( XVec, this->rowStartInVector() ) );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    if ( !M_bcDirichletMarkerManagement.hasMarkerDirichletBCelimination() ) return;

    this->log("Heat","updateJacobianDofElimination","start" );

    this->updateDofEliminationIds( "temperature", data );

    this->log("Heat","updateJacobianDofElimination","finish" );
}


} // namespace FeelModels
} // namespace Feel
