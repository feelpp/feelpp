/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/electric/electric.hpp>

namespace Feel {
namespace FeelModels {

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    vector_ptrtype& U = data.initialGuess();
    this->updateNewtonInitialGuess( data, this->modelContext( U, this->rowStartInVector() ) );
}
ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateJacobian( data, this->modelContext( XVec, this->rowStartInVector() ) );
}

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    if ( !M_boundaryConditions.hasTypeDofElimination() )
        return;

    this->log("Electric","updateJacobianDofElimination","start" );

    this->updateDofEliminationIds( "electric-potential", data );

    this->log("Electric","updateJacobianDofElimination","finish" );
}

} // namespace FeelModels
} // namespace Feel
