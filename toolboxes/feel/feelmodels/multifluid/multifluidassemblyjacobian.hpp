/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef _MULTIFLUID_ASSEMBLYJACOBIAN_HPP
#define _MULTIFLUID_ASSEMBLYJACOBIAN_HPP 1

#include <feel/feelmodels/multifluid/multifluid.hpp>

namespace Feel {
namespace FeelModels {

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    this->log("MultiFluid","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();

    // TODO: build startBlockSpaceIndexLevelsets only once
    std::vector<size_type> startBlockSpaceIndexLevelsets;
    std::transform( M_levelsetModels.begin(), M_levelsetModels.end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } 
            );

    auto mctx = this->modelContext( U, M_fluidModel->startBlockSpaceIndexVector(), startBlockSpaceIndexLevelsets );

    M_fluidModel->updateNewtonInitialGuess( data, mctx );
    for( levelset_model_ptrtype const& lsModel: M_levelsetModels )
        lsModel->updateNewtonInitialGuess( data, mctx );

    this->log("MultiFluid","updateNewtonInitialGuess","finish" );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc = (buildCstPart)?" (cst)":" (non cst)";
    this->log("MultiFluid", "updateJacobian", "start"+sc);
    this->timerTool("Solve").start();

    // TODO: build startBlockSpaceIndexLevelsets only once
    std::vector<size_type> startBlockSpaceIndexLevelsets;
    std::transform( M_levelsetModels.begin(), M_levelsetModels.end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } 
            );

    auto mctx = this->modelContext( XVec, M_fluidModel->startBlockSpaceIndexVector(), startBlockSpaceIndexLevelsets );

    M_fluidModel->updateJacobian( data, mctx );
    for( levelset_model_ptrtype const& lsModel: M_levelsetModels )
        lsModel->updateJacobian( data, mctx );

    //// Update interface forces
    //this->updateJacobianInterfaceForces( data );
    //// Update inextensibility
    //this->updateJacobianInextensibility( data );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","updateJacobian", fmt::format( "finish in {} s", timeElapsed ) );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    M_fluidModel->updateJacobianDofElimination( data );
    for( levelset_model_ptrtype const& lsModel: M_levelsetModels )
        lsModel->updateJacobianDofElimination( data );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobian_Fluid( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();

    // TODO: build startBlockSpaceIndexLevelsets only once
    std::vector<size_type> startBlockSpaceIndexLevelsets;
    std::transform( M_levelsetModels.begin(), M_levelsetModels.end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } 
            );

    auto mctx = this->modelContext( vecCurrentSolution, M_fluidModel->startBlockSpaceIndexVector(),
                                    M_algebraicBlockVectorSolutionLevelsets, startBlockSpaceIndexLevelsets );

    //bool currentStabDoGLSDoAssembly = M_fluidModel->stabilizationGLSDoAssembly();
    //M_fluidModel->setStabilizationGLSDoAssembly( true );
    M_fluidModel->updateJacobian( data, mctx );
    //M_fluidModel->setStabilizationGLSDoAssembly( currentStabDoGLSDoAssembly );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobian_Levelset( size_type lsModelIndex, DataUpdateJacobian & data ) const
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

    M_levelsetModels[lsModelIndex]->updateJacobian( data, mctx );
}

} // namespace FeelModels
} // namespace Feel

#endif
