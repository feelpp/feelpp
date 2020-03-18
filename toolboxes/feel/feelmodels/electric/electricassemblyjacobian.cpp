/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels/electric/electric.hpp>

namespace Feel {
namespace FeelModels {

ELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
ELECTRIC_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    if ( M_bcDirichlet.empty() ) return;

    this->log("Electric","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto mesh = this->mesh();
    size_type startBlockIndexElectricPotential = this->startSubBlockSpaceIndex( "potential-electric" );
    auto v = this->spaceElectricPotential()->element( U, this->rowStartInVector()+startBlockIndexElectricPotential );
    auto se = this->symbolsExpr();

    for( auto const& d : M_bcDirichlet )
    {
        v.on(_range=markedfaces(mesh, M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) ) ),
             _expr=expression(d,se) );
    }

    // update info for synchronization
    this->updateDofEliminationIds( "potential-electric", data );

    this->log("Electric","updateNewtonInitialGuess","finish" );
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
    if ( this->M_bcDirichlet.empty() ) return;

    this->log("Electric","updateJacobianDofElimination","start" );

    this->updateDofEliminationIds( "potential-electric", data );

    this->log("Electric","updateJacobianDofElimination","finish" );
}


} // namespace FeelModels
} // namespace Feel
