/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/heat/heat.hpp>

namespace Feel {
namespace FeelModels {

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto mctx = this->modelContext( XVec, this->rowStartInVector() );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSol = data.vectorInfo( "time-stepping.previous-solution");
        auto mctxPrevious = this->modelContextNoTrialSymbolsExpr( previousSol,this->rowStartInVector() );
        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }

    this->updateLinearPDE( data, mctx );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    this->updateLinearPDEDofElimination( data, this->modelContext() );
}

} // namespace FeelModels
} // namespace Feel
