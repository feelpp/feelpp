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

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobianInterfaceForces( DataUpdateJacobian & data ) const
{
    if( this->hasInterfaceForces() )
    {
        this->log("MultiFluid", "updateJacobianInterfaceForces", "start: update interface forces");
        this->timerTool("Solve").start();
        if( M_hasInterfaceForcesModel )
        {
            this->timerTool("Solve").start();
            for( auto const& lsInterfaceForces: M_levelsetInterfaceForcesModels )
            {
                for( auto const& force: lsInterfaceForces.second )
                {
                    if( force.second )
                    {
                        force.second->updateFluidInterfaceForcesJacobian( data );
                    }
                }
            }

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateJacobianInterfaceForces", fmt::format( "update interface (model) forces in {} s", timeElapsedInterfaceForces ) );
        }

        if( M_additionalInterfaceForcesModel.size() > 0 )
        {
            this->timerTool("Solve").start();
            for( auto const& f: M_additionalInterfaceForcesModel )
                f.second->updateFluidInterfaceForcesJacobian( data );

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateJacobianInterfaceForces", fmt::format( "update additional interface forces in {} s", timeElapsedInterfaceForces ) );
        }

        double timeElapsed = this->timerTool("Solve").stop();
        this->log( "MultiFluid", "updateJacobianInterfaceForces", 
                fmt::format( "interface forces updated in {} s", timeElapsed ) );
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobianInextensibility( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();

    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;

    if( this->M_enableInextensibility )
    {
        if( BuildNonCstPart )
        {
            auto mesh = this->mesh();
            auto XhV = this->fluidModel()->functionSpaceVelocity();

            auto myBfV = form2(
                    _test=XhV,_trial=XhV,_matrix=J,
                    //_pattern=size_type(Pattern::COUPLED),
                    _rowstart=this->fluidModel()->rowStartInMatrix(),
                    _colstart=this->fluidModel()->colStartInMatrix()
                    );

            auto u = XhV->element(XVec, this->fluidModel()->rowStartInVector());
            auto const& v = this->fluidModel()->fieldVelocity();
            auto Id = vf::Id<nDim, nDim>();

            for( size_type i = 0; i < M_levelsetModels.size(); ++i )
            {
                if( this->hasInextensibility(i) && this->inextensibilityMethod(i) == "penalty" )
                {
                    auto N = this->levelsetModel(i)->N();
                    auto NxN = idv(N)*trans(idv(N));
                    auto D = this->levelsetModel(i)->D();

                    this->timerTool("Solve").start();

                    if( BuildNonCstPart )
                    {
                        myBfV += integrate(
                                _range=elements(mesh),
                                _expr=this->M_inextensibilityGamma[i]*trace((Id-NxN)*gradt(u))*trace((Id-NxN)*grad(v))*idv(D)/h(),
                                _geomap=this->geomap()
                                );
                    }

                    double timeElapsedInextensibility_Penalty = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateJacobianInextensibility",
                            fmt::format( "assembly inextensibility (penalty) in {} s", timeElapsedInextensibility_Penalty ) );
                }
            }

            if( this->hasInextensibilityLM() )
            {
                CHECK( this->hasStartSubBlockSpaceIndex("inextensibility-lm") ) << " start dof index for inextensibility-lm is not present\n";
                this->timerTool("Solve").start();

                size_type startBlockIndexInextensibilityLM = this->startSubBlockSpaceIndex("inextensibility-lm");
                auto lambda = this->functionSpaceInextensibilityLM()->element();

                auto inextensibleLevelsetsExpr = Feel::FeelModels::globalLevelsetExpr( M_inextensibleLevelsets );
                auto inextensibleLevelsets = vf::project(
                        _space=this->M_levelsetSpaceManager->functionSpaceScalar(), 
                        _range=this->M_levelsetSpaceManager->rangeMeshElements(),
                        _expr=inextensibleLevelsetsExpr
                        );
                auto inextensibleLevelsetsDeltaExpr = Feel::FeelModels::levelsetDelta( _element=inextensibleLevelsets, _thickness=M_globalLevelsetThicknessInterface );
                auto N = trans(gradv(inextensibleLevelsets)) / sqrt( gradv(inextensibleLevelsets)*trans(gradv(inextensibleLevelsets)) );
                auto NxN = N*trans(N);

                form2( _trial=this->functionSpaceInextensibilityLM(), _test=this->fluidModel()->functionSpaceVelocity(), 
                        _matrix=J,
                        _rowstart=this->fluidModel()->rowStartInMatrix(),
                        _colstart=startBlockIndexInextensibilityLM ) +=
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=idt(lambda)*trace((Id-NxN)*grad(v))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );
                form2( _trial=this->fluidModel()->functionSpaceVelocity(), _test=this->functionSpaceInextensibilityLM(), 
                        _matrix=J,
                        _rowstart=startBlockIndexInextensibilityLM,
                        _colstart=this->fluidModel()->colStartInMatrix() ) +=
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=id(lambda)*trace((Id-NxN)*gradt(u))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );

                double timeElapsedInextensibility_LagrangeMult = this->timerTool("Solve").stop();
                this->log("MultiFluid","updateJacobianInextensibility",
                        fmt::format( "assembly inextensibility (lagrange-multiplier) in {} s", timeElapsedInextensibility_LagrangeMult ) );
            }
        }
    }
}

} // namespace FeelModels
} // namespace Feel

#endif
