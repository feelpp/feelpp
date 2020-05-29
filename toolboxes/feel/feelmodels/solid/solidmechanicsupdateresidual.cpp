/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/solid/solidmechanics.hpp>

namespace Feel
{
namespace FeelModels
{


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateResidual( data, this->modelContext( XVec, this->rowStartInVector() ) );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    vector_ptrtype& U = data.initialGuess();
    this->updateNewtonInitialGuess( data, this->modelContext( U, this->rowStartInVector() ) );
}

//--------------------------------------------------------------------------------------------------//


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( !M_bcDirichletMarkerManagement.hasMarkerDirichletBCelimination() )
        return;

    this->log("SolidMechanics","updateResidualDofElimination","start" );

    this->updateDofEliminationIds( "displacement", data );

    this->log("SolidMechanics","updateResidualDofElimination","finish" );
}


} // FeelModels
} // Feel



