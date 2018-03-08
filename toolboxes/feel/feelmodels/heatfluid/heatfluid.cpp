/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2016-12-12

  Copyright (C) 2016 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file thermoelectric.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-12-12
 */

#include <feel/feelmodels/heatfluid/heatfluid.hpp>

#include <feel/feelvf/vf.hpp>
/*#include <feel/feelvf/form.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/operators.hpp>
 #include <feel/feelvf/operations.hpp>*/

#include <feel/feelmodels/modelmesh/createmesh.hpp>

namespace Feel
{
namespace FeelModels
{

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
HEATFLUID_CLASS_TEMPLATE_TYPE::HeatFluid( std::string const& prefix,
                                          bool buildMesh,
                                          WorldComm const& worldComm,
                                          std::string const& subPrefix,
                                          ModelBaseRepository const& modelRep )
    :
    super_type( prefix, worldComm, subPrefix, modelRep )
{
    this->log("HeatFluid","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatFluidConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatFluidSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatFluidPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatFluidTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"HeatFluid.info") );
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // build mesh
    if ( buildMesh )
        this->initMesh();
    //-----------------------------------------------------------------------------//
    this->log("HeatFluid","constructor", "finish");
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_useNaturalConvection = false;//"Forced-Convection";
    // M_solverNewtonInitialGuessUseLinearThermoElectric = boption(_prefix=this->prefix(),_name="solver-newton.initial-guess.use-linear-thermo-electric");
    // M_solverNewtonInitialGuessUseLinearHeatTransfer = boption(_prefix=this->prefix(),_name="solver-newton.initial-guess.use-linear-heat-transfer");
    // M_solverNewtonInitialGuessUseLinearElectric = boption(_prefix=this->prefix(),_name="solver-newton.initial-guess.use-linear-electric");
    // if ( M_solverNewtonInitialGuessUseLinearThermoElectric )
    // {
    //     M_solverNewtonInitialGuessUseLinearHeatTransfer = true;
    //     M_solverNewtonInitialGuessUseLinearElectric = true;
    // }
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("HeatFluid","initMesh", "start");
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("initMesh");
    this->log("HeatFluid","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()



HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
int
HEATFLUID_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = M_heatTransferModel->nBlockMatrixGraph() + M_fluidModel->nBlockMatrixGraph();
    return nBlock;
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
HEATFLUID_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);

    int nBlockFluid = M_fluidModel->nBlockMatrixGraph();
    int nBlockHeatTransfer = M_heatTransferModel->nBlockMatrixGraph();

    int startIndexBlockFluid = 0;
    int startIndexBlockHeat = nBlockFluid;

    auto blockMatFluid = M_fluidModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockFluid ;++tk1 )
        for (int tk2=0;tk2<nBlockFluid ;++tk2 )
            myblockGraph(startIndexBlockFluid+tk1,startIndexBlockFluid+tk2) = blockMatFluid(tk1,tk2);

    BlocksStencilPattern patCoupling1(1,M_fluidModel->functionSpace()->nSpaces,size_type(Pattern::ZERO));
    patCoupling1(0,0) = size_type(Pattern::COUPLED);
    myblockGraph(startIndexBlockHeat,startIndexBlockFluid) = stencil(_test=M_heatTransferModel->spaceTemperature(),
                                                                     _trial=M_fluidModel->functionSpace(),
                                                                     _pattern_block=patCoupling1,
                                                                     _diag_is_nonzero=false,_close=false)->graph();

    BlocksStencilPattern patCoupling2(M_fluidModel->functionSpace()->nSpaces,1,size_type(Pattern::ZERO));
    patCoupling2(0,0) = size_type(Pattern::COUPLED);
    myblockGraph(startIndexBlockFluid,startIndexBlockHeat) = stencil(_test=M_fluidModel->functionSpace(),
                                                                     _trial=M_heatTransferModel->spaceTemperature(),
                                                                     _pattern_block=patCoupling2,
                                                                     _diag_is_nonzero=false,_close=false)->graph();

    auto blockMatHeatTransfer = M_heatTransferModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockHeatTransfer ;++tk1 )
        for (int tk2=0;tk2<nBlockHeatTransfer ;++tk2 )
            myblockGraph(startIndexBlockHeat+tk1,startIndexBlockHeat+tk2) = blockMatHeatTransfer(tk1,tk2);

    myblockGraph.close();

    return myblockGraph;
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("HeatFluid","init", "start" );
    this->timerTool("Constructor").start();

    M_heatTransferModel.reset( new heattransfer_model_type(prefixvm(this->prefix(),"heat-transfer"), false, this->worldComm(),
                                                           this->subPrefix(), this->repository() ) );
    if ( !M_heatTransferModel->modelPropertiesPtr() )
        M_heatTransferModel->setModelProperties( this->modelPropertiesPtr() );
    M_heatTransferModel->loadMesh( this->mesh() );
    M_heatTransferModel->init( false );

    M_fluidModel.reset( new fluid_model_type(prefixvm(this->prefix(),"fluid"), false, this->worldComm(),
                                             this->subPrefix(), this->repository() ) );
    if ( !M_fluidModel->modelPropertiesPtr() )
        M_fluidModel->setModelProperties( this->modelPropertiesPtr() );
    M_fluidModel->setMesh( this->mesh() );
    M_fluidModel->init( false );

    if ( !this->isStationary() )
    {
        this->updateTime( this->timeStepBase()->time() );
        this->setTimeInitial( this->timeStepBase()->timeInitial() );
    }

    if ( !M_useNaturalConvection )
    {
        M_heatTransferModel->initAlgebraicFactory();
        M_fluidModel->initAlgebraicFactory();
    }

#if 0
    M_rangeMeshElements = ( M_heatTransferModel->thermalProperties()->isDefinedOnWholeMesh() && M_electricModel->electricProperties()->isDefinedOnWholeMesh() )?
        elements(this->mesh() ) :
        intersect( M_heatTransferModel->rangeMeshElements(), M_electricModel->rangeMeshElements() );
#endif

    // for ( auto const& rangeData : M_electricModel->electricProperties()->rangeMeshElementsByMaterial() )
    // {
    //     std::string const& matName = rangeData.first;
    //     if ( !M_heatTransferModel->thermalProperties()->hasMaterial( matName ) )
    //         continue;
    //     M_rangeMeshElementsByMaterial[matName] = rangeData.second;
    // }
    // post-process
    this->initPostProcess();

    // backend : use worldComm of Xh
    M_backendMonolithic = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldComm() );

    size_type currentStartIndex = 0;
    M_startBlockIndexFieldsInMatrix["fluid"] = currentStartIndex;
    currentStartIndex += 2;// TODO!!!!!!
    M_startBlockIndexFieldsInMatrix["heat-transfer"] = currentStartIndex;

    // vector solution
    int nBlockFluid = M_fluidModel->blockVectorSolution().size();
    int nBlockHeatTransfer = M_heatTransferModel->blockVectorSolution().size();
    int nBlock = nBlockFluid + nBlockHeatTransfer;
    M_blockVectorSolutionMonolithic.resize( nBlock );
    int indexBlock=0;
    for ( int tk1=0;tk1<nBlockFluid ;++tk1 )
        M_blockVectorSolutionMonolithic(indexBlock+tk1) = M_fluidModel->blockVectorSolution()(tk1);
    indexBlock += nBlockFluid;
    for ( int tk1=0;tk1<nBlockHeatTransfer ;++tk1 )
        M_blockVectorSolutionMonolithic(indexBlock+tk1) = M_heatTransferModel->blockVectorSolution()(tk1);
    indexBlock += nBlockHeatTransfer;
    // init petsc vector associated to the block
    M_blockVectorSolutionMonolithic.buildVector( this->backend() );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
    {
        if ( M_useNaturalConvection /*M_solverName == "Newton"*/ )
        {
            M_algebraicFactoryMonolithic.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );
        }
    }

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("HeatFluid","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("HeatFluid","initPostProcess", "start");
    this->timerTool("Constructor").start();

    std::string modelName = "heat-fluid";
    auto const& exportsFields = this->modelProperties().postProcess().exports( modelName ).fields();
    M_postProcessFieldExportedHeatTransfert = M_heatTransferModel->postProcessFieldExported( exportsFields, "heat-transfer" );
    M_postProcessFieldExportedFluid = M_fluidModel->postProcessFieldExported( exportsFields, "fluid" );

    if ( !M_postProcessFieldExportedHeatTransfert.empty() || !M_postProcessFieldExportedFluid.empty() )
    {
        std::string geoExportType="static";//change_coords_only, change, static
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export",
                               _geo=geoExportType,
                               _path=this->exporterPath() );

        if ( this->doRestart() && this->restartPath().empty() )
        {
            if ( M_exporter->doExport() )
                M_exporter->restart(this->timeInitial());
        }
    }

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("HeatFluid","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}



HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
boost::shared_ptr<std::ostringstream>
HEATFLUID_CLASS_TEMPLATE_TYPE::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << M_heatTransferModel->getInfo()->str();
    *_ostr << M_fluidModel->getInfo()->str();

    std::string myexporterType = M_exporter->type();
    int myexporterFreq = M_exporter->freq();

    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||---------------Info : HeatFluid---------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository();
    //*_ostr << "\n   Physical Model"
    //       << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient");
    //*_ostr << "\n   Numerical Solver"
    //       << "\n     -- solver   : " << M_solverName;
    *_ostr << "\n   Exporter"
           << "\n     -- type            : " << myexporterType
           << "\n     -- freq save       : " << myexporterFreq
           << "\n     -- fields [heat-transfer] : TODO"// << doExport_str
           << "\n     -- fields [electric] : TODO" //<< doExport_str
           << "\n   Processors"
           << "\n     -- number of proc : " << this->worldComm().globalSize()
           << "\n     -- current rank : " << this->worldComm().globalRank();

    if ( M_algebraicFactoryMonolithic )
        *_ostr << M_algebraicFactoryMonolithic->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->heatTransferModel()->updateTimeStep();
    this->fluidModel()->updateTimeStep();
    this->updateTime( this->timeStepBase()->time() );
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    if ( !M_exporter->doExport() ) return;

    this->log("HeatFluid","exportResults", "start");
    this->timerTool("PostProcessing").start();

    bool hasFieldToExportHeatTransfer = M_heatTransferModel->updateExportedFields( M_exporter,M_postProcessFieldExportedHeatTransfert,time );
    bool hasFieldToExportFluid = M_fluidModel->updateExportedFields( M_exporter,M_postProcessFieldExportedFluid,time );
    if ( hasFieldToExportHeatTransfer || hasFieldToExportFluid )
        M_exporter->save();

    M_heatTransferModel->exportResults( time );
    M_fluidModel->exportResults( time );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("HeatFluid","exportResults", "finish");
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    M_heatTransferModel->updateParameterValues();
    M_fluidModel->updateParameterValues();
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("HeatFluid","solve", "start");
    this->timerTool("Solve").start();

    if ( !M_useNaturalConvection )
    {
        M_fluidModel->solve();
        M_heatTransferModel->setFieldVelocityConvectionIsUsed( true );
        M_heatTransferModel->updateFieldVelocityConvection( idv(M_fluidModel->fieldVelocity()) );
        M_heatTransferModel->solve();
    }
    else
    {
        this->updateParameterValues();

    }
#if 0
    if ( M_solverName == "Linear" )
    {
        M_electricModel->setRowStartInMatrix( 0 );
        M_electricModel->setColStartInMatrix( 0 );
        M_electricModel->setRowStartInVector( 0 );
        M_electricModel->solve();
        M_heatTransferModel->solve();
    }
    else if ( M_solverName == "Newton" )
    {
        // initial guess
        if ( M_solverNewtonInitialGuessUseLinearElectric )
        {
            M_electricModel->setRowStartInMatrix( 0 );
            M_electricModel->setColStartInMatrix( 0 );
            M_electricModel->setRowStartInVector( 0 );
            M_electricModel->solve();
        }
        if ( M_solverNewtonInitialGuessUseLinearHeatTransfer )
            M_heatTransferModel->solve();

        // solve non linear monolithic system
        int nBlockHeatTransfer = M_heatTransferModel->nBlockMatrixGraph();
        M_electricModel->setRowStartInMatrix( nBlockHeatTransfer );
        M_electricModel->setColStartInMatrix( nBlockHeatTransfer );
        M_electricModel->setRowStartInVector( nBlockHeatTransfer );

        M_blockVectorSolutionMonolithic.updateVectorFromSubVectors();
        M_algebraicFactoryMonolithic->solve( "Newton", M_blockVectorSolutionMonolithic.vectorMonolithic() );
        M_blockVectorSolutionMonolithic.localize();

        M_electricModel->updateElectricField();
    }

#endif

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("HeatFluid","solve", (boost::format("finish in %1% s")%tElapsed).str() );

}



HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
    this->log("HeatFluid","updateNewtonInitialGuess","start" );
    M_heatTransferModel->updateNewtonInitialGuess( U );
    M_fluidModel->updateNewtonInitialGuess( U );
    this->log("HeatFluid","updateNewtonInitialGuess","finish" );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
#if 0
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("HeatFluid","updateJacobian", "start"+sc);
    size_type startBlockIndexTemperature = M_startBlockIndexFieldsInMatrix.find( "temperature" )->second;
    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;

    auto mesh = this->mesh();

    if ( !buildCstPart )
    {
        auto XhV = M_electricModel->spaceElectricPotential();
        auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );
        auto XhT = M_heatTransferModel->spaceTemperature();
        //auto const& t = M_heatTransferModel->fieldTemperature();
        auto const t = XhT->element(XVec, this->rowStartInVector()+startBlockIndexTemperature );

        auto mybfTT = form2( _test=XhT,_trial=XhT,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexTemperature,
                             _colstart=this->colStartInMatrix()+startBlockIndexTemperature );
        auto mybfTV = form2( _test=XhT,_trial=XhV,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexTemperature,
                             _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );
        auto mybfVV = form2( _test=XhV,_trial=XhV,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                             _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );
        auto mybfVT = form2( _test=XhV,_trial=XhT,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                             _colstart=this->colStartInMatrix()+startBlockIndexTemperature );

        for ( auto const& rangeData : M_rangeMeshElementsByMaterial )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            auto const& electricConductivity = M_electricModel->electricProperties()->electricConductivity( matName );
            if ( electricConductivity.isConstant() )
            {
                if ( M_modelUseJouleEffect )
                {
                    double sigma = electricConductivity.value();
                    mybfTV +=
                        integrate( _range=range,
                                   _expr= -sigma*2*inner(gradt(v),gradv(v))*id( t ),
                                   _geomap=this->geomap() );
                }
            }
            else
            {
                std::string symbolStr = "heat_transfer_T";
                auto sigma = electricConductivity.expr( symbolStr, idv(t) );
                if ( M_modelUseJouleEffect )
                {
                    mybfTV +=
                        integrate( _range=range,
                                   _expr= -sigma*2*inner(gradt(v),gradv(v))*id( t ),
                                   _geomap=this->geomap() );
                }

                if ( sigma.expression().hasSymbol( symbolStr ) )
                {
                    auto sigmaDiff = diff( electricConductivity.expr(),symbolStr,1,"",this->worldComm(),this->repository().expr());
                    auto sigmaDiffEval = expr( sigmaDiff, symbolStr, idv(t) );

                    if ( M_modelUseJouleEffect )
                    {
                        mybfTT +=
                            integrate( _range=range,
                                       _expr= -sigmaDiffEval*idt(t)*inner(gradv(v)/*,gradv(v)*/)*id( t ),
                                       _geomap=this->geomap() );
                    }

                    mybfVV +=
                        integrate( _range=range,
                                   _expr= sigma*inner(gradt(v),grad(v)),
                                   _geomap=this->geomap() );
                    mybfVT +=
                        integrate( _range=range,
                                   _expr= sigmaDiffEval*idt(t)*inner(gradv(v),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }
    }

    DataUpdateJacobian dataSubPhysics( data );
    dataSubPhysics.setDoBCStrongDirichlet( false );
    M_heatTransferModel->updateJacobian( dataSubPhysics );
    M_electricModel->updateJacobian( dataSubPhysics );

    if ( buildNonCstPart && doBCStrongDirichlet )
    {
        M_heatTransferModel->updateJacobianStrongDirichletBC( J,RBis );
        M_electricModel->updateJacobianStrongDirichletBC( J,RBis );
    }
    this->log("HeatFluid","updateJacobian", "finish"+sc);
#endif
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
#if 0
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool useJacobianLinearTerms = data.useJacobianLinearTerms();
    bool doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("HeatFluid","updateResidual", "start"+sc);

    size_type startBlockIndexTemperature = M_startBlockIndexFieldsInMatrix.find( "temperature" )->second;
    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;

    auto mesh = this->mesh();

    DataUpdateResidual dataSubPhysics( data );
    dataSubPhysics.setDoBCStrongDirichlet( false );
    M_heatTransferModel->updateResidual( dataSubPhysics );
    M_electricModel->updateResidual( dataSubPhysics );

    if ( !buildCstPart )
    {
        auto XhV = M_electricModel->spaceElectricPotential();
        auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );
        auto XhT = M_heatTransferModel->spaceTemperature();
        auto const t = XhT->element(XVec, this->rowStartInVector()+startBlockIndexTemperature );
        auto mylfT = form1( _test=XhT, _vector=R,
                            _rowstart=this->rowStartInVector()+startBlockIndexTemperature );
        auto mylfV = form1( _test=XhV, _vector=R,
                            _rowstart=this->rowStartInVector()+startBlockIndexElectricPotential );

        for ( auto const& rangeData : M_rangeMeshElementsByMaterial )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            auto const& electricConductivity = M_electricModel->electricProperties()->electricConductivity( matName );
            if ( electricConductivity.isConstant() )
            {
                if ( M_modelUseJouleEffect )
                {
                    double sigma = electricConductivity.value();
                    mylfT +=
                        integrate( _range=range,
                                   _expr= -sigma*inner(gradv(v)/*,gradv(v)*/)*id( t ),
                                   _geomap=this->geomap() );
                }
            }
            else
            {
                std::string symbolStr = "heat_transfer_T";
                auto sigma = electricConductivity.expr( symbolStr, idv(t) );
                if ( M_modelUseJouleEffect )
                {
                    mylfT +=
                        integrate( _range=range,
                                   _expr= -sigma*inner(gradv(v)/*,gradv(v)*/)*id( t ),
                                   _geomap=this->geomap() );
                }
                if ( sigma.expression().hasSymbol( symbolStr ) )
                {
                    mylfV +=
                        integrate( _range=range,
                                   _expr= sigma*inner(gradv(v),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }
    }

    if ( !buildCstPart && doBCStrongDirichlet &&
         ( M_heatTransferModel->hasMarkerDirichletBCelimination() || M_electricModel->hasMarkerDirichletBCelimination() ) )
    {
        R->close();
        M_heatTransferModel->updateResidualStrongDirichletBC( R );
        M_electricModel->updateResidualStrongDirichletBC( R );
    }
    this->log("HeatFluid","updateResidual", "finish"+sc);
#endif
}


} // end namespace FeelModels
} // end namespace Feel
