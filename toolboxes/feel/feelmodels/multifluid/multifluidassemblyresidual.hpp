/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef _MULTIFLUID_ASSEMBLYRESIDUAL_HPP
#define _MULTIFLUID_ASSEMBLYRESIDUAL_HPP 1

#include <feel/feelmodels/multifluid/multifluid.hpp>

namespace Feel {
namespace FeelModels {

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool useJacobianLinearTerms = data.useJacobianLinearTerms();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("MultiFluid","updateResidual", "start"+sc);
    this->timerTool("Solve").start();

    // TODO: build startBlockSpaceIndexLevelsets only once
    std::vector<size_type> startBlockSpaceIndexLevelsets;
    std::transform( M_levelsetModels.begin(), M_levelsetModels.end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } 
            );

    auto mctx = this->modelContext( XVec, M_fluidModel->startBlockSpaceIndexVector(), startBlockSpaceIndexLevelsets );

    M_fluidModel->updateResidual( data, mctx );
    for( levelset_model_ptrtype const& lsModel: M_levelsetModels )
        lsModel->updateResidual( data, mctx );

    //// Update interface forces
    //this->updateResidualInterfaceForces( data );

    //// Update inextensibility
    //this->updateResidualInextensibility( data );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","updateResidual","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    M_fluidModel->updateResidualDofElimination( data );
    for( levelset_model_ptrtype const& lsModel: M_levelsetModels )
        lsModel->updateResidualDofElimination( data );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidual_Fluid( DataUpdateResidual & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();

    // TODO: build startBlockSpaceIndexLevelsets only once
    std::vector<size_type> startBlockSpaceIndexLevelsets;
    std::transform( M_levelsetModels.begin(), M_levelsetModels.end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } 
            );

    auto mctx = this->modelContext( vecCurrentSolution, M_fluidModel->startBlockSpaceIndexVector(),
                                    M_algebraicBlockVectorSolutionLevelsets, startBlockSpaceIndexLevelsets );

    M_fluidModel->updateResidual( data, mctx );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidual_Levelset( size_type lsModelIndex, DataUpdateResidual & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();

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

    M_levelsetModels[lsModelIndex]->updateResidual( data, mctx );
}

} // namespace FeelModels
} // namespace Feel

#endif

