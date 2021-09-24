/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef _MULTIFLUID_ASSEMBLYLINEAR_HPP
#define _MULTIFLUID_ASSEMBLYLINEAR_HPP 1

#include <feel/feelmodels/multifluid/multifluid.hpp>

namespace Feel {
namespace FeelModels {

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc = (buildCstPart)?" (cst)":" (non cst)";
    this->log("MultiFluid", "updateLinearPDE", "start"+sc);
    this->timerTool("Solve").start();

    // TODO: build startBlockSpaceIndexLevelsets only once
    std::vector<size_type> startBlockSpaceIndexLevelsets;
    std::transform( M_levelsetModels.begin(), M_levelsetModels.end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } 
            );

    auto mctx = this->modelContext( XVec, M_fluidModel->startBlockSpaceIndexVector(), startBlockSpaceIndexLevelsets );

    M_fluidModel->updateLinearPDE( data, mctx );
    for( levelset_model_ptrtype const& lsModel: M_levelsetModels )
        lsModel->updateLinearPDE( data, mctx );

    //// Update interface forces
    //this->updateLinearPDEInterfaceForces( data );

    //// Update inextensibility
    //this->updateLinearPDEInextensibility( data );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","updateLinearPDE","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    M_fluidModel->updateLinearPDEDofElimination( data );
    for( levelset_model_ptrtype const& lsModel: M_levelsetModels )
        lsModel->updateLinearPDEDofElimination( data );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinear_Fluid( DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("MultiFluid","updateLinear_Fluid", "start"+sc);

    // TODO: build startBlockSpaceIndexLevelsets only once
    std::vector<size_type> startBlockSpaceIndexLevelsets;
    std::transform( M_levelsetModels.begin(), M_levelsetModels.end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } 
            );

    auto mctx = this->modelContext( vecCurrentSolution, M_fluidModel->startBlockSpaceIndexVector(),
                                    M_algebraicBlockVectorSolutionLevelsets, startBlockSpaceIndexLevelsets );
    //if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    //{
        //auto previousSolFluid = data.vectorInfo( "time-stepping.previous-solution");
        //auto mctxPrevious = this->modelContextNoTrialSymbolsExpr(
            //{ { "solution", std::make_tuple( M_timeStepThetaSchemePreviousSolution, this->startSubBlockSpaceIndex("heat") ) } },
            //{
                //{ "solution", std::make_tuple( previousSolFluid, this->fluidModel()->startBlockSpaceIndexVector()) },
                //{ "velocity_extrapolated", std::make_tuple( M_fluidModel->vectorPreviousVelocityExtrapolated(), 0 ) }
            //} );
        //mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    //}

    M_fluidModel->updateLinearPDE( data, mctx );

    this->log("MultiFluid","updateLinear_Fluid", "finish"+sc);
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinear_Levelset( size_type lsModelIndex, DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("MultiFluid","updateLinear_Levelset", "start "+std::to_string(lsModelIndex)+sc);

    std::vector<vector_ptrtype> vecSolsLevelsets = M_algebraicBlockVectorSolutionLevelsets;
    // TODO: build startBlockSpaceIndexLevelsets only once
    std::vector<size_type> startBlockSpaceIndexLevelsets;
    std::transform( M_levelsetModels.begin(), M_levelsetModels.end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } 
            );
    // Set current levelsetModel vector
    vecSolsLevelsets[lsModelIndex] = vecCurrentSolution;

    auto mctx = this->modelContext( 
            M_fluidModel->algebraicBlockVectorSolution()->vectorMonolithic(), M_fluidModel->startBlockSpaceIndexVector(),
            vecSolsLevelsets, startBlockSpaceIndexLevelsets 
            );
    //if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    //{
        //auto previousSolFluid = data.vectorInfo( "time-stepping.previous-solution");
        //auto mctxPrevious = this->modelContextNoTrialSymbolsExpr(
            //{ { "solution", std::make_tuple( M_timeStepThetaSchemePreviousSolution, this->startSubBlockSpaceIndex("heat") ) } },
            //{
                //{ "solution", std::make_tuple( previousSolFluid, this->fluidModel()->startBlockSpaceIndexVector()) },
                //{ "velocity_extrapolated", std::make_tuple( M_fluidModel->vectorPreviousVelocityExtrapolated(), 0 ) }
            //} );
        //mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    //}

    M_levelsetModel[lsModelIndex]->updateLinearPDE( data, mctx );

    this->log("MultiFluid","updateLinear_Levelset", "finish "+std::to_string(lsModelIndex)+sc);
}

} // namespace FeelModels
} // namespace Feel

#endif
