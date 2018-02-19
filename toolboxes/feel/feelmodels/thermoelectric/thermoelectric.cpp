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

#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>

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

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::ThermoElectric( std::string const& prefix,
                                                    bool buildMesh,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    ModelBaseRepository const& modelRep )
    :
    super_type( prefix, worldComm, subPrefix, modelRep )
{
    this->log("ThermoElectric","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoElectricConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoElectricSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoElectricPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoElectricTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"ThermoElectric.info") );
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // build mesh
    if ( buildMesh )
        this->createMesh();
    //-----------------------------------------------------------------------------//
    this->log("ThermoElectric","constructor", "finish");
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{

}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::createMesh()
{
    this->log("ThermoElectric","createMesh", "start");
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("createMesh");
    this->log("ThermoElectric","createMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = M_heatTransferModel->nBlockMatrixGraph() + M_electricModel->nBlockMatrixGraph();
    return nBlock;
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);

    int indexBlock=0;

    int nBlockHeatTransfer = M_heatTransferModel->nBlockMatrixGraph();
    auto blockMatHeatTransfer = M_heatTransferModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockHeatTransfer ;++tk1 )
        for (int tk2=0;tk2<nBlockHeatTransfer ;++tk2 )
            myblockGraph(indexBlock+tk1,indexBlock+tk2) = blockMatHeatTransfer(tk1,tk2);

    BlocksStencilPattern patCoupling1(1,nBlockHeatTransfer,size_type(Pattern::ZERO));
    patCoupling1(0,0) = size_type(Pattern::COUPLED);
    myblockGraph(indexBlock,indexBlock+nBlockHeatTransfer) = stencil(_test=M_heatTransferModel->spaceTemperature(),
                                                                     _trial=M_electricModel->spaceElectricPotential(),
                                                                     _pattern_block=patCoupling1,
                                                                     _diag_is_nonzero=false,_close=false)->graph();

    indexBlock += nBlockHeatTransfer;

    int nBlockElectric = M_electricModel->nBlockMatrixGraph();
    auto blockMatElectric = M_electricModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockElectric ;++tk1 )
        for (int tk2=0;tk2<nBlockElectric ;++tk2 )
            myblockGraph(indexBlock+tk1,indexBlock+tk2) = blockMatElectric(tk1,tk2);

    myblockGraph.close();

    return myblockGraph;
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("ThermoElectric","init", "start" );
    this->timerTool("Constructor").start();

    M_modelName = this->modelProperties().models().model("thermo-electric").equations();

    M_solverName = "Linear";
    if ( M_modelName == "ThermoElectric-linear" )
    {
        M_modelName = "ThermoElectric";
        M_solverName = "Linear";
    }
    else if ( M_modelName == "ThermoElectric-nonlinear" )
    {
        M_modelName = "ThermoElectric";
        M_solverName = "Newton";
    }
    CHECK( M_modelName == "ThermoElectric" ) << "invalid model name : " << M_modelName << "\n";


    bool useSubSolver = (M_solverName == "Linear");

    M_heatTransferModel.reset( new heattransfer_model_type(prefixvm(this->prefix(),"thermo"), false, this->worldComm(),
                                                           this->subPrefix(), this->repository() ) );
    M_heatTransferModel->loadMesh( this->mesh() );
    M_heatTransferModel->init( useSubSolver );
    if ( M_solverName == "Linear" )
    {
        M_heatTransferModel->algebraicFactory()->addFunctionLinearPreAssemblyNonCst = boost::bind( &self_type::updateLinearPreAssemblyJouleLaw,
                                                                                                   boost::ref( *this ), _1, _2 );
    }

    M_electricModel.reset( new electric_model_type(prefixvm(this->prefix(),"electric"), false, this->worldComm(),
                                                   this->subPrefix(), this->repository() ) );
    M_electricModel->setMesh( this->mesh() );
    M_electricModel->init( useSubSolver );


    M_rangeMeshElements = ( M_heatTransferModel->thermalProperties()->isDefinedOnWholeMesh() && M_electricModel->electricProperties()->isDefinedOnWholeMesh() )?
        elements(this->mesh() ) :
        intersect( M_heatTransferModel->rangeMeshElements(), M_electricModel->rangeMeshElements() );

    // post-process
    this->initPostProcess();

    // backend : use worldComm of Xh
    M_backendMonolithic = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldComm() );

    size_type currentStartIndex = 0;
    M_startBlockIndexFieldsInMatrix["temperature"] = currentStartIndex;
    currentStartIndex += M_heatTransferModel->nBlockMatrixGraph();
    M_startBlockIndexFieldsInMatrix["potential-electric"] = currentStartIndex;

    // vector solution
    int nBlock = this->nBlockMatrixGraph();
    M_blockVectorSolutionMonolithic.resize( nBlock );
    int indexBlock=0;

    int nBlockHeatTransfer = M_heatTransferModel->nBlockMatrixGraph();
    for ( int tk1=0;tk1<nBlockHeatTransfer ;++tk1 )
        M_blockVectorSolutionMonolithic(indexBlock+tk1) = M_heatTransferModel->blockVectorSolution()(tk1);
    indexBlock += nBlockHeatTransfer;
    int nBlockElectric = M_electricModel->nBlockMatrixGraph();
    for ( int tk1=0;tk1<nBlockElectric ;++tk1 )
        M_blockVectorSolutionMonolithic(indexBlock+tk1) = M_electricModel->blockVectorSolution()(tk1);
    indexBlock += nBlockElectric;

    // init petsc vector associated to the block
    M_blockVectorSolutionMonolithic.buildVector( this->backend() );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
    {
        if ( !useSubSolver )
        {
            // matrix graph of non zero
            typename model_algebraic_factory_type::graph_ptrtype graph( new typename model_algebraic_factory_type::graph_type( this->buildBlockMatrixGraph() ) );
            // tool for assembly and solver
            M_algebraicFactoryMonolithic.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend(),
                                                                                  graph, graph->mapRow().indexSplit() ) );
        }
    }

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("ThermoElectric","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("ThermoElectric","initPostProcess", "start");
    this->timerTool("Constructor").start();

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

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("ThermoElectric","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
boost::shared_ptr<std::ostringstream>
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << M_heatTransferModel->getInfo()->str();
    *_ostr << M_electricModel->getInfo()->str();

    std::string myexporterType = M_exporter->type();
    int myexporterFreq = M_exporter->freq();

    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||----------Info : ThermoElectric---------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository();
    //*_ostr << "\n   Physical Model"
    //       << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient");
    *_ostr << "\n   Numerical Solver"
           << "\n     -- solver   : " << M_solverName;
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



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    if ( !M_exporter->doExport() ) return;

    this->log("ThermoElectric","exportResults", "start");
    this->timerTool("PostProcessing").start();

    bool hasFieldToExportHeatTransfer = M_heatTransferModel->updateExportedFields( M_exporter,time );
    bool hasFieldToExportElectric = M_electricModel->updateExportedFields( M_exporter,time );
    if ( hasFieldToExportHeatTransfer || hasFieldToExportElectric )
        M_exporter->save();

    M_heatTransferModel->exportMeasures( time );
    M_electricModel->exportMeasures( time );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("ThermoElectric","exportResults", "finish");
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    M_heatTransferModel->updateParameterValues();
    M_electricModel->updateParameterValues();
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("ThermoElectric","solve", "start");
    this->timerTool("Solve").start();

    this->updateParameterValues();
    //std::string theThermoElectricModel = this->modelProperties().model();
    //if ( theThermoElectricModel == "ThermoElectric-linear" )
    if ( M_solverName == "Linear" )
    {
        M_electricModel->setRowStartInMatrix( 0 );
        M_electricModel->setColStartInMatrix( 0 );
        M_electricModel->setRowStartInVector( 0 );

        //auto mySolutionVectorLinear = M_backendElectricModel->toBackendVectorPtr( this->fieldElectricPotential() );
        //M_algebraicFactoryElectricModel->solve( "LinearSystem", mySolutionVectorLinear );
        M_electricModel->solve();
        M_heatTransferModel->solve();
    }
    else if ( M_solverName == "Newton" )
    {
        int nBlockHeatTransfer = M_heatTransferModel->nBlockMatrixGraph();
        M_electricModel->setRowStartInMatrix( nBlockHeatTransfer );
        M_electricModel->setColStartInMatrix( nBlockHeatTransfer );
        M_electricModel->setRowStartInVector( nBlockHeatTransfer );

        M_blockVectorSolutionMonolithic.updateVectorFromSubVectors();
        M_algebraicFactoryMonolithic->solve( "Newton", M_blockVectorSolutionMonolithic.vectorMonolithic() );
        M_blockVectorSolutionMonolithic.localize();

        M_electricModel->updateElectricField();
    }


    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("ThermoElectric","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
#if 0
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("ThermoElectric","updateLinearPDE", "start"+sc);
    boost::mpi::timer thetimer;

    auto mesh = this->mesh();
    auto XhV = this->spaceElectricPotential();
    auto const& v = this->fieldElectricPotential();

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix() ,
                                              _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=XhV, _vector=F,
                               _rowstart=this->rowStartInVector() );

    //--------------------------------------------------------------------------------------------------//

    auto sigma = idv(M_electricProperties->fieldElectricConductivity());
    if ( buildCstPart )
    {
        bilinearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= sigma*inner(gradt(v),grad(v)),
                       _geomap=this->geomap() );
    }

    // update source term
    this->updateSourceTermLinearPDE(F, buildCstPart);

    // update bc
    this->updateWeakBCLinearPDE(A,F,buildCstPart);

    if ( !buildCstPart && _doBCStrongDirichlet)
    {
        this->updateBCStrongDirichletLinearPDE(A,F);
    }
#endif
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearPreAssemblyJouleLaw( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const
{
    this->log("ThermoElectric","updateLinearPreAssemblyJouleLaw","start" );

    auto mesh = this->mesh();
    auto XhV = M_electricModel->spaceElectricPotential();
    auto const& v = M_electricModel->fieldElectricPotential();
    auto XhT = M_heatTransferModel->spaceTemperature();
    auto const& t = M_heatTransferModel->fieldTemperature();
    auto sigma = idv(M_electricModel->electricProperties()->fieldElectricConductivity());
    form1( _test=XhT,_vector=F,
           _rowstart=M_heatTransferModel->rowStartInVector() ) +=
        integrate( _range=M_rangeMeshElements,
                       _expr= sigma*inner(gradv(v),gradv(v))*id(t),
                       _geomap=this->geomap() );

    this->log("ThermoElectric","updateLinearPreAssemblyJouleLaw","finish" );
}



THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
    this->log("ThermoElectric","updateNewtonInitialGuess","start" );
    M_heatTransferModel->updateNewtonInitialGuess( U );
    M_electricModel->updateNewtonInitialGuess( U );
    this->log("ThermoElectric","updateNewtonInitialGuess","finish" );
}
THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("ThermoElectric","updateJacobian", "start"+sc);
    size_type startBlockIndexTemperature = M_startBlockIndexFieldsInMatrix.find( "temperature" )->second;
    size_type startBlockIndexElectricPotential = M_startBlockIndexFieldsInMatrix.find( "potential-electric" )->second;

    auto mesh = this->mesh();

    if ( !buildCstPart )
    {
        auto XhV = M_electricModel->spaceElectricPotential();
        auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );
        auto XhT = M_heatTransferModel->spaceTemperature();
        auto const& t = M_heatTransferModel->fieldTemperature();

        auto sigma = idv(M_electricModel->electricProperties()->fieldElectricConductivity());

        form2( _test=XhT,_trial=XhV,_matrix=J,
               _pattern=size_type(Pattern::COUPLED),
               _rowstart=this->rowStartInMatrix()+startBlockIndexTemperature,
               _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential  ) +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -sigma*2*inner(gradt(v),gradv(v))*id( t ),
                       _geomap=this->geomap() );
    }

    DataUpdateJacobian dataSubPhysics( data );
    dataSubPhysics.setDoBCStrongDirichlet( false );
    M_heatTransferModel->updateJacobian( dataSubPhysics );
    M_electricModel->updateJacobian( dataSubPhysics );

    if ( buildNonCstPart && doBCStrongDirichlet )
    {
        M_heatTransferModel->updateBCStrongDirichletJacobian( J,RBis );
        M_electricModel->updateBCStrongDirichletJacobian( J,RBis );
    }
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool useJacobianLinearTerms = data.useJacobianLinearTerms();
    bool doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("ThermoElectric","updateResidual", "start"+sc);

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
        auto sigma = idv(M_electricModel->electricProperties()->fieldElectricConductivity());
        auto myLinearFormThermo = form1( _test=XhT, _vector=R,
                                         _rowstart=this->rowStartInVector() );

        myLinearFormThermo +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -sigma*inner(gradv(v),gradv(v))*id( t ),
                       _geomap=this->geomap() );
    }

    if ( !buildCstPart && doBCStrongDirichlet &&
         ( M_heatTransferModel->hasMarkerDirichletBCelimination() || M_electricModel->hasMarkerDirichletBCelimination() ) )
    {
        R->close();
        M_heatTransferModel->updateBCDirichletStrongResidual( R );
        M_electricModel->updateBCDirichletStrongResidual( R );
    }
}


} // end namespace FeelModels
} // end namespace Feel
