/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels/heat/heat.hpp>

namespace Feel {
namespace FeelModels {

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateResidual( data, this->modelContext( XVec, this->rowStartInVector() ) );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( !M_bcDirichletMarkerManagement.hasMarkerDirichletBCelimination() ) return;

    this->log("Heat","updateResidualDofElimination","start" );

    this->updateDofEliminationIds( "temperature", data );

    this->log("Heat","updateResidualDofElimination","finish" );
}

} // namespace FeelModels
} // namespace Feel
