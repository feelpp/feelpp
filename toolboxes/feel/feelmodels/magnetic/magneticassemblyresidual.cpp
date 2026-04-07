/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/magnetic/magnetic.hpp>

namespace Feel {
namespace FeelModels {

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateResidual( data, this->modelContext( XVec, this->rowStartInVector() ) );
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;

    this->log("Magnetic","updateResidualDofElimination","start" );

    this->updateDofEliminationIds( "vector_potential", data );

    this->log("Magnetic","updateResidualDofElimination","finish" );
}

} // namespace FeelModels
} // namespace Feel
