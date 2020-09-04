/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto mctx = this->modelContext( XVec, this->rowStartInVector() );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSol = data.vectorInfo( "time-stepping.previous-solution");
        auto mctxPrevious = this->modelContext/*NoTrialSymbolsExpr*/( {
                { "solution", std::make_tuple( previousSol, this->rowStartInVector()) },
                { "velocity_extrapolated", std::make_tuple( M_vectorPreviousVelocityExtrapolated, 0 ) }
            } );
        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }
    this->updateLinearPDE( data, mctx );
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    this->updateLinearPDEDofElimination( data, this->modelContext() );
}


} // end namespace FeelModels
} // end namespace Feel


