/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/electric/electric.hpp>

namespace Feel {
namespace FeelModels {

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateResidual( data, this->modelContext( XVec, this->rowStartInVector() ) );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( !M_boundaryConditions.hasTypeDofElimination() )
        return;

    this->log("Electric","updateResidualDofElimination","start" );

    this->updateDofEliminationIds( "potential-electric", data );

    this->log("Electric","updateResidualDofElimination","finish" );
}

} // namespace FeelModels
} // namespace Feel
