/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef _MULTIFLUID_ASSEMBLYLINEAR_HPP
#define _MULTIFLUID_ASSEMBLYLINEAR_HPP 1

#include <feel/feelmodels/multifluid/multifluid.hpp>

namespace Feel {
namespace FeelModels {

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

    // TODO: build vecSolsLevelsets and startBlockSpaceIndexLevelsets only once
    std::vector<vector_ptrtype> vecSolsLevelsets;
    std::transform( this->levelsetModels().begin(), this->levelsetModels().end(), std::back_inserter( vecSolsLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->algebraicBlockVectorSolution()->vectorMonolithic(); } 
            );
    std::vector<size_type> startBlockSpaceIndexLevelsets;
    std::transform( this->levelsetModels().begin(), this->levelsetModels().end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } 
            );

    auto mctx = this->modelContext( vecCurrentSolution, this->fluidModel()->startBlockSpaceIndexVector(),
                                    vecSolsLevelsets, startBlockSpaceIndexLevelsets );
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

    // TODO: build vecSolsLevelsets and startBlockSpaceIndexLevelsets only once
    std::vector<vector_ptrtype> vecSolsLevelsets;
    std::transform( this->levelsetModels().begin(), this->levelsetModels().end(), std::back_inserter( vecSolsLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->algebraicBlockVectorSolution()->vectorMonolithic(); } 
            );
    std::vector<size_type> startBlockSpaceIndexLevelsets;
    std::transform( this->levelsetModels().begin(), this->levelsetModels().end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } 
            );
    // Set current levelsetModel vector
    vecSolsLevelsets[lsModelIndex] = vecCurrentSolution;

    auto mctx = this->modelContext( 
            this->fluidModel()->algebraicBlockVectorSolution()->vectorMonolithic(), this->fluidModel()->startBlockSpaceIndexVector(),
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

    this->levelsetModel( lsModelIndex )->updateLinearPDE( data, mctx );

    this->log("MultiFluid","updateLinear_Levelset", "finish "+std::to_string(lsModelIndex)+sc);
}

} // namespace FeelModels
} // namespace Feel

#endif
