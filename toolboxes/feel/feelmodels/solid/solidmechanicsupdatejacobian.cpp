/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/solid/solidmechanics.hpp>

namespace Feel
{
namespace FeelModels
{

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateJacobian( data, this->modelContext( XVec, this->rowStartInVector() ) );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    if ( !M_boundaryConditions.hasTypeDofElimination() )
        return;

    this->log("SolidMechanics","updateJacobianDofElimination","start" );

    this->updateDofEliminationIds( "displacement", data );

    this->log("SolidMechanics","updateJacobianDofElimination","finish" );
}

} // FeelModels

} // Feel



