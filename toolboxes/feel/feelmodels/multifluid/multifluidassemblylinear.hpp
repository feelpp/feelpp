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
    this->log("MultiFluid", "updateLinearPDE", fmt::format( "finish in {} s", timeElapsed ) );
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

    M_levelsetModels[lsModelIndex]->updateLinearPDE( data, mctx );

    this->log("MultiFluid","updateLinear_Levelset", fmt::format( "finish {}{}", lsModelIndex, sc ) );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEInterfaceForces( DataUpdateLinear & data ) const
{
    // Update interface forces
    if( this->hasInterfaceForces() )
    {
        this->log("MultiFluid", "updateLinearPDEInterfaceForces", "start: update interface forces");
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
                        force.second->updateFluidInterfaceForcesLinearPDE( data );
                    }
                }
            }

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateLinearPDEInterfaceForces", fmt::format( "update interface (model) forces in {} s", timeElapsedInterfaceForces ) );
        }

        if( M_additionalInterfaceForcesModel.size() > 0 )
        {
            this->timerTool("Solve").start();
            for( auto const& f: M_additionalInterfaceForcesModel )
                f.second->updateFluidInterfaceForcesLinearPDE( data );

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateLinearPDEInterfaceForces", fmt::format( "update additional interface forces in {} s", timeElapsedInterfaceForces ) );
        }

        double timeElapsed = this->timerTool("Solve").stop();
        this->log( "MultiFluid", "updateLinearPDEInterfaceForces", 
                fmt::format( "interface forces updated in {} s", timeElapsed ) );
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEInextensibility( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;

    if( this->M_enableInextensibility )
    {
        if( BuildNonCstPart )
        {
            auto mesh = this->mesh();
            auto XhV = this->fluidModel()->functionSpaceVelocity();

            auto const& u = this->fluidModel()->fieldVelocity();
            auto const& v = u;
            auto Id = vf::Id<nDim, nDim>();

            auto myBfV = form2( 
                    _test=XhV, _trial=XhV, _matrix=A,
                    _rowstart=this->fluidModel()->rowStartInMatrix(),
                    _colstart=this->fluidModel()->colStartInMatrix()
                    );

            for( size_type i = 0; i < M_levelsetModels.size(); ++i )
            {
                if( this->hasInextensibility(i) && this->inextensibilityMethod(i) == "penalty" )
                {
                    auto N = this->levelsetModel(i)->N();
                    auto NxN = idv(N)*trans(idv(N));
                    auto D = this->levelsetModel(i)->D();

                    this->timerTool("Solve").start();

                    myBfV += integrate(
                            _range=elements(mesh),
                            _expr=this->M_inextensibilityGamma[i]*trace((Id-NxN)*gradt(u))*trace((Id-NxN)*grad(v))*idv(D)/h(),
                            _geomap=this->geomap()
                            );

                    double timeElapsedInextensibility_Penalty = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateLinearPDEInextensibility",
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
                auto inextensibleLevelsetsDeltaExpr = Feel::FeelModels::levelsetDelta( inextensibleLevelsets, M_globalLevelsetThicknessInterface );
                auto N = trans(gradv(inextensibleLevelsets)) / sqrt( gradv(inextensibleLevelsets)*trans(gradv(inextensibleLevelsets)) );
                auto NxN = N*trans(N);

                form2( _trial=this->functionSpaceInextensibilityLM(), _test=XhV, 
                        _matrix=A,
                        _rowstart=this->fluidModel()->rowStartInMatrix(),
                        _colstart=startBlockIndexInextensibilityLM ) +=
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=idt(lambda)*trace((Id-NxN)*grad(v))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );
                form2( _trial=XhV, _test=this->functionSpaceInextensibilityLM(), 
                        _matrix=A,
                        _rowstart=startBlockIndexInextensibilityLM,
                        _colstart=this->fluidModel()->colStartInMatrix() ) +=
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=id(lambda)*trace((Id-NxN)*gradt(u))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );

                double timeElapsedInextensibility_LagrangeMult = this->timerTool("Solve").stop();
                this->log("MultiFluid","updateLinearPDEInextensibility",
                        fmt::format( "assembly inextensibility (lagrange-multiplier) in {} s", timeElapsedInextensibility_LagrangeMult ) );
            }
        }
    }
}

} // namespace FeelModels
} // namespace Feel

#endif
