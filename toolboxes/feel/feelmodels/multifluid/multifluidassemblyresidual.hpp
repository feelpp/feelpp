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
    this->log("MultiFluid","updateResidual", fmt::format( "finish in {} s", timeElapsed ) );
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

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidualInterfaceForces( DataUpdateResidual & data ) const
{
    if( this->hasInterfaceForces() )
    {
        this->log("MultiFluid", "updateResidualInterfaceForces", "start: update interface forces");
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
                        force.second->updateFluidInterfaceForcesResidual( data );
                    }
                }
            }

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateResidualInterfaceForces", fmt::format( "update interface (model) forces in {} s", timeElapsedInterfaceForces ) );
        }

        if( M_additionalInterfaceForcesModel.size() > 0 )
        {
            this->timerTool("Solve").start();
            for( auto const& f: M_additionalInterfaceForcesModel )
                f.second->updateFluidInterfaceForcesResidual( data );

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateResidualInterfaceForces", fmt::format( "update additional interface forces in {} s", timeElapsedInterfaceForces ) );
        }

        double timeElapsed = this->timerTool("Solve").stop();
        this->log( "MultiFluid", "updateResidualInterfaceForces", 
                fmt::format( "interface forces updated in {} s", timeElapsed ) );
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidualInextensibility( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    if( this->M_enableInextensibility )
    {
        if( BuildNonCstPart )
        {
            auto mesh = this->mesh();
            auto XhV = this->fluidModel()->functionSpaceVelocity();

            auto myLV = form1( 
                    _test=XhV,_vector=R,
                    //_pattern=size_type(Pattern::COUPLED),
                    _rowstart=this->fluidModel()->rowStartInVector()
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

                    if( !UseJacobianLinearTerms )
                    {
                        myLV += integrate(
                                _range=elements(mesh),
                                _expr=this->M_inextensibilityGamma[i]*trace((Id-NxN)*gradv(u))*trace((Id-NxN)*grad(v))*idv(D)/h(),
                                _geomap=this->geomap()
                                );
                    }

                    double timeElapsedInextensibility_Penalty = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateResidualInextensibility",
                            fmt::format( "assembly inextensibility (penalty) in {} s", timeElapsedInextensibility_Penalty ) );
                }
            }

            if( this->hasInextensibilityLM() )
            {
                CHECK( this->hasStartSubBlockSpaceIndex("inextensibility-lm") ) << " start dof index for inextensibility-lm is not present\n";
                this->timerTool("Solve").start();

                size_type startBlockIndexInextensibilityLM = this->startSubBlockSpaceIndex("inextensibility-lm");
                auto lambda = this->functionSpaceInextensibilityLM()->element(XVec,startBlockIndexInextensibilityLM);

                auto inextensibleLevelsetsExpr = Feel::FeelModels::globalLevelsetExpr( M_inextensibleLevelsets );
                auto inextensibleLevelsets = vf::project(
                        _space=this->M_levelsetSpaceManager->functionSpaceScalar(), 
                        _range=this->M_levelsetSpaceManager->rangeMeshElements(),
                        _expr=inextensibleLevelsetsExpr
                        );
                auto inextensibleLevelsetsDeltaExpr = Feel::FeelModels::levelsetDelta( _element=inextensibleLevelsets, _thickness=M_globalLevelsetThicknessInterface );
                auto N = trans(gradv(inextensibleLevelsets)) / sqrt( gradv(inextensibleLevelsets)*trans(gradv(inextensibleLevelsets)) );
                auto NxN = N*trans(N);

                form1( _test=this->fluidModel()->functionSpaceVelocity(), _vector=R,
                        _rowstart=this->fluidModel()->rowStartInVector() ) +=
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=idv(lambda)*trace((Id-NxN)*grad(v))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );
                form1( _test=this->functionSpaceInextensibilityLM(), _vector=R,
                        _rowstart=startBlockIndexInextensibilityLM ) += 
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=id(lambda)*trace((Id-NxN)*gradv(u))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );

                double timeElapsedInextensibility_LagrangeMult = this->timerTool("Solve").stop();
                this->log("MultiFluid","updateResidualInextensibility",
                        fmt::format( "assembly inextensibility (lagrange-multiplier) in {} s", timeElapsedInextensibility_LagrangeMult ) );
            }
        }
    }
}

} // namespace FeelModels
} // namespace Feel

#endif

