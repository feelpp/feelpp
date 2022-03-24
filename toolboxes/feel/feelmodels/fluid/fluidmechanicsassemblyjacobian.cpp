/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    vector_ptrtype& U = data.initialGuess();
    this->updateNewtonInitialGuess( data, this->modelContext( U, this->rowStartInVector() ) );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateJacobian( data, this->modelContext( XVec, this->rowStartInVector() ) );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;

    this->log("FluidMechanics","updateJacobianDofElimination", "start" );

    this->timerTool("Solve").start();

    this->updateDofEliminationIds( "velocity", data );

    if ( !M_boundaryConditions->pressureImposed().empty() )
    {
        this->updateDofEliminationIds( "pressurelm1", this->dofEliminationIds( "pressurebc-lm" ), data );
        if ( nDim == 3 )
            this->updateDofEliminationIds( "pressurelm2", this->dofEliminationIds( "pressurebc-lm" ), data );
    }

    for ( auto const& [bpname,bpbc] : M_bodySetBC )
    {
        if ( bpbc.hasTranslationalVelocityExpr() )
        {
            std::string spaceName = "body-bc."+bpbc.name()+".translational-velocity";
            //this->updateDofEliminationIds( spaceName, this->dofEliminationIds( "body-bc.translational-velocity" ), data );
            this->updateDofEliminationIds( spaceName, data );
        }
        if ( bpbc.hasAngularVelocityExpr() )
        {
            std::string spaceName = "body-bc."+bpbc.name()+".angular-velocity";
            //this->updateDofEliminationIds( spaceName, this->dofEliminationIds( "body-bc.angular-velocity" ), data );
            this->updateDofEliminationIds( spaceName, data );
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("FluidMechanics","updateJacobianDofElimination","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

} // namespace FeelModels
} // namespace Feel


