/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    this->updateResidual( data, this->modelContext( XVec, this->rowStartInVector() ) );
} // updateResidual


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( !this->hasStrongDirichletBC() ) return;

    this->updateDofEliminationIds( "velocity", data );

    if ( this->hasMarkerPressureBC() )
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
            this->updateDofEliminationIds( spaceName, this->dofEliminationIds( "body-bc.translational-velocity" ), data );
        }
        if ( bpbc.hasAngularVelocityExpr() )
        {
            std::string spaceName = "body-bc."+bpbc.name()+".angular-velocity";
            this->updateDofEliminationIds( spaceName, this->dofEliminationIds( "body-bc.angular-velocity" ), data );
        }
    }
}




} // namespace FeelModels
} // namespace Feel


