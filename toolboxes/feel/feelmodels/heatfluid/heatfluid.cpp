/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2018-03-06

  Copyright (C) 2018 Feel++ Consortium

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

#include <feel/feelmodels/heatfluid/heatfluid.hpp>

#include <feel/feelvf/vf.hpp>
/*#include <feel/feelvf/form.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/operators.hpp>
 #include <feel/feelvf/operations.hpp>*/

#include <feel/feelmodels/modelmesh/createmesh.hpp>
#include <feel/feelmodels/modelcore/utils.hpp>

namespace Feel
{
namespace FeelModels
{

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
HEATFLUID_CLASS_TEMPLATE_TYPE::HeatFluid( std::string const& prefix,
                                          std::string const& keyword,
                                          worldcomm_ptr_t const& worldComm,
                                          std::string const& subPrefix,
                                          ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelPhysics<mesh_type::nDim>( "heat-fluid" ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep )
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

    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    this->log("HeatFluid","constructor", "finish");
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_useNaturalConvection = boption(_prefix=this->prefix(),_name="use-natural-convection");//"Forced-Convection";

    M_BoussinesqRefTemperature = doption(_name="Boussinesq.ref-temperature",_prefix=this->prefix());
    std::string gravityStr;
    if ( Environment::vm().count(prefixvm(this->prefix(),"gravity-force").c_str()) )
        gravityStr = soption(_name="gravity-force",_prefix=this->prefix());
    else if (nDim == 2 )
        gravityStr = "{0,-9.80665}";
    else if (nDim == 3 )
        gravityStr = "{0,0,-9.80665}";
    M_gravityForce = expr<nDim,1,2>( gravityStr,"",this->worldComm(),this->repository().expr() );

    M_useSemiImplicitTimeScheme = boption(_name="use-semi-implicit-time-scheme",_prefix=this->prefix());
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
    int nBlock = M_heatModel->nBlockMatrixGraph() + M_fluidModel->nBlockMatrixGraph();
    return nBlock;
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
HEATFLUID_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);

    int nBlockFluid = M_fluidModel->nBlockMatrixGraph();
    int nBlockHeat = M_heatModel->nBlockMatrixGraph();

    int startIndexBlockFluid = 0;
    int startIndexBlockHeat = nBlockFluid;

    auto blockMatFluid = M_fluidModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockFluid ;++tk1 )
        for (int tk2=0;tk2<nBlockFluid ;++tk2 )
            myblockGraph(startIndexBlockFluid+tk1,startIndexBlockFluid+tk2) = blockMatFluid(tk1,tk2);

    myblockGraph(startIndexBlockHeat,startIndexBlockFluid) = stencil(_test=M_heatModel->spaceTemperature(),
                                                                     _trial=M_fluidModel->functionSpaceVelocity(),
                                                                     _diag_is_nonzero=false,_close=false)->graph();

    myblockGraph(startIndexBlockFluid,startIndexBlockHeat) = stencil(_test=M_fluidModel->functionSpaceVelocity(),
                                                                     _trial=M_heatModel->spaceTemperature(),
                                                                     _diag_is_nonzero=false,_close=false)->graph();
    if ( M_fluidModel->stabilizationGLS() )
    {
        myblockGraph(startIndexBlockFluid+1,startIndexBlockHeat) = stencil(_test=M_fluidModel->functionSpaceVelocity(),
                                                                           _trial=M_heatModel->spaceTemperature(),
                                                                           _diag_is_nonzero=false,_close=false)->graph();
    }

    auto blockMatHeat = M_heatModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockHeat ;++tk1 )
        for (int tk2=0;tk2<nBlockHeat ;++tk2 )
            myblockGraph(startIndexBlockHeat+tk1,startIndexBlockHeat+tk2) = blockMatHeat(tk1,tk2);

    myblockGraph.close();

    return myblockGraph;
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("HeatFluid","init", "start" );
    this->timerTool("Constructor").start();

    M_heatModel = std::make_shared<heat_model_type>(prefixvm(this->prefix(),"heat"), "heat", this->worldCommPtr(),
                                                    this->subPrefix(), this->repository() );
    M_fluidModel = std::make_shared<fluid_model_type>(prefixvm(this->prefix(),"fluid"), "fluid", this->worldCommPtr(),
                                                      this->subPrefix(), this->repository() );

    if ( this->physics().empty() )
    {
        typename ModelPhysics<mesh_type::nDim>::subphysic_description_type subPhyicsDesc;
        subPhyicsDesc[M_heatModel->physicType()] = std::make_tuple( M_heatModel->keyword(),M_heatModel );
        subPhyicsDesc[M_fluidModel->physicType()] = std::make_tuple( M_fluidModel->keyword(), M_fluidModel );
        this->initPhysics( this->keyword(), this->modelProperties().models(), subPhyicsDesc );
    }

    // physical properties
    if ( !M_materialsProperties )
    {
        auto paramValues = this->modelProperties().parameters().toParameterValues();
        this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    if ( !M_mesh )
        this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    // init heat toolbox
    M_heatModel->setPhysics( this->physics( M_heatModel->physicType() ), M_heatModel->keyword() );
    M_heatModel->setManageParameterValues( false );
    if ( !M_heatModel->modelPropertiesPtr() )
    {
        M_heatModel->setModelProperties( this->modelPropertiesPtr() );
        M_heatModel->setManageParameterValuesOfModelProperties( false );
    }
    M_heatModel->setMesh( this->mesh() );
    M_heatModel->setMaterialsProperties( M_materialsProperties );
    std::string velConvExprStr = (boost::format( nDim==2? "{%1%_U_0,%1%_U_1}:%1%_U_0:%1%_U_1":"{%1%_U_0,%1%_U_1,%1%_U_2}:%1%_U_0:%1%_U_1:%1%_U_2" )%M_fluidModel->keyword() ).str();
    auto velConvExpr = expr<nDim,1>( velConvExprStr, "",this->worldComm(),this->repository().expr() );
    for ( std::string const& matName : M_materialsProperties->physicToMaterials( this->physicsAvailableFromCurrentType() ) )
        M_heatModel->setVelocityConvectionExpr( matName,velConvExpr );
    M_heatModel->init( false );

    // init fluid toolbox
    M_fluidModel->setPhysics( this->physics( M_fluidModel->physicType() ), M_fluidModel->keyword() );
    M_fluidModel->setManageParameterValues( false );
    if ( !M_fluidModel->modelPropertiesPtr() )
    {
        M_fluidModel->setModelProperties( this->modelPropertiesPtr() );
        M_fluidModel->setManageParameterValuesOfModelProperties( false );
    }
    M_fluidModel->setMesh( this->mesh() );
    M_fluidModel->init( false );

    if ( !this->isStationary() )
    {
        this->updateTime( this->timeStepBase()->time() );
        this->setTimeInitial( this->timeStepBase()->timeInitial() );
    }
    else
        M_useSemiImplicitTimeScheme = false;

    if ( !M_useNaturalConvection || M_useSemiImplicitTimeScheme )
    {
        // init velocity convection (TODO : the initial velocity from the fluid)
        //M_heatModel->setFieldVelocityConvectionIsUsed( true );
        //M_heatModel->updateForUseFunctionSpacesVelocityConvection();
        M_heatModel->initAlgebraicFactory();
        M_fluidModel->initAlgebraicFactory();
        M_heatModel->algebraicFactory()->setFunctionLinearAssembly( boost::bind( &self_type::updateLinear_Heat,
                                                                                 boost::ref( *this ), _1 ) );
        M_heatModel->algebraicFactory()->setFunctionResidualAssembly( boost::bind( &self_type::updateResidual_Heat,
                                                                                   boost::ref( *this ), _1 ) );
        M_heatModel->algebraicFactory()->setFunctionJacobianAssembly( boost::bind( &self_type::updateJacobian_Heat,
                                                                                   boost::ref( *this ), _1 ) );

        if ( M_useSemiImplicitTimeScheme )
        {
            M_fluidModel->setSolverName("Oseen");
            M_fluidModel->setStabilizationGLSDoAssembly( false );
            M_fluidModel->algebraicFactory()->addFunctionLinearAssembly( boost::bind( &self_type::updateLinearFluidSolver,
                                                                                      boost::ref( *this ), _1 ) );
            M_fluidModel->algebraicFactory()->addFunctionResidualAssembly( boost::bind( &self_type::updateResidualFluidSolver,
                                                                                        boost::ref( *this ), _1 ) );

            if ( M_fluidModel->timeStepping() == "Theta" && M_fluidModel->stabilizationGLS() )
            {
                CHECK( M_heatModel->timeStepping() == "Theta" ) << "not implemented";
                std::string timeSteppingPreviousSolutionDataKey_heat = prefixvm( this->heatModel()->prefix(),"time-stepping.previous-solution");
                this->fluidModel()->algebraicFactory()->dataInfos().addVectorInfo( timeSteppingPreviousSolutionDataKey_heat,
                                                                                  this->heatModel()->algebraicFactory()->dataInfos().vectorInfo( timeSteppingPreviousSolutionDataKey_heat ) );
            }
        }
    }
    else
    {
        M_fluidModel->setStabilizationGLSDoAssembly( false );
        // if ( M_useSemiImplicitTimeScheme )
        //     M_fluidModel->setSolverName("Oseen");
    }

    // post-process
    this->initPostProcess();

    // update constant parameters
    this->updateParameterValues();

    // backend
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() );

    // block vector solution
    auto blockVectorSolutionFluid = M_fluidModel->blockVectorSolution();
    auto blockVectorSolutionHeat = M_heatModel->blockVectorSolution();
    int nBlockFluid = blockVectorSolutionFluid.size();
    int nBlockHeat = blockVectorSolutionHeat.size();
    int nBlock = nBlockFluid + nBlockHeat;
    M_blockVectorSolution.resize( nBlock );
    int indexBlock=0;
    int numberOfBlockSpaceFluid = 0;
    for ( int k=0;k<nBlockFluid ;++k )
    {
        M_blockVectorSolution(indexBlock+k) = blockVectorSolutionFluid(k);
        numberOfBlockSpaceFluid += blockVectorSolutionFluid(k)->map().numberOfDofIdToContainerId();
    }
    indexBlock += nBlockFluid;
    for ( int k=0;k<nBlockHeat ;++k )
        M_blockVectorSolution(indexBlock+k) = blockVectorSolutionHeat(k);
    indexBlock += nBlockHeat;

    size_type currentStartBlockSpaceIndex = 0;
    this->setStartSubBlockSpaceIndex( "fluid", currentStartBlockSpaceIndex );
    currentStartBlockSpaceIndex += numberOfBlockSpaceFluid;
    this->setStartSubBlockSpaceIndex( "heat", currentStartBlockSpaceIndex );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
    {
        if ( M_useNaturalConvection && !M_useSemiImplicitTimeScheme )
        {
            // init monolithic vector associated to the block vector
            M_blockVectorSolution.buildVector( this->backend() );

            M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );
            if ( M_fluidModel->hasOperatorPCD() )
                M_algebraicFactory->preconditionerTool()->attachOperatorPCD( "pcd", M_fluidModel->operatorPCD() );

            if ( ( M_fluidModel->timeStepping() == "Theta" ) || ( M_heatModel->timeStepping() == "Theta" ) )
            {
                M_timeStepThetaSchemePreviousContrib = this->backend()->newVector(M_blockVectorSolution.vectorMonolithic()->mapPtr() );
                M_algebraicFactory->addVectorResidualAssembly( M_timeStepThetaSchemePreviousContrib, 1.0, "Theta-Time-Stepping-Previous-Contrib", true );
                M_algebraicFactory->addVectorLinearRhsAssembly( M_timeStepThetaSchemePreviousContrib, -1.0, "Theta-Time-Stepping-Previous-Contrib", false );
                if ( M_fluidModel->stabilizationGLS() || M_heatModel->stabilizationGLS() )
                {
                    auto vecPreviousSolution =  this->backend()->newVector( M_blockVectorSolution.vectorMonolithic()->mapPtr() );
                    M_algebraicFactory->dataInfos().addVectorInfo( prefixvm( this->prefix(),"time-stepping.previous-solution"), vecPreviousSolution );
                    M_algebraicFactory->dataInfos().addVectorInfo( prefixvm( M_fluidModel->prefix(),"time-stepping.previous-solution"), vecPreviousSolution );
                    M_algebraicFactory->dataInfos().addVectorInfo( prefixvm( M_heatModel->prefix(),"time-stepping.previous-solution"), vecPreviousSolution );
                }
            }
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

    // need to not include export fields of material of subphysics
    std::set<std::string> ppExportsAllFieldsAvailableHeat = Feel::FeelModels::detail::set_difference( this->heatModel()->postProcessExportsAllFieldsAvailable(),
                                                                                                      this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->heatModel()->physicsAvailable() ) );
    std::set<std::string> ppExportsAllFieldsAvailablefluid = Feel::FeelModels::detail::set_difference( this->fluidModel()->postProcessExportsAllFieldsAvailable(),
                                                                                                       this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->fluidModel()->physicsAvailable() ) );
   std::set<std::string> ppExportsAllFieldsAvailable;
   for ( auto const& s : ppExportsAllFieldsAvailableHeat )
       ppExportsAllFieldsAvailable.insert( prefixvm( this->heatModel()->keyword(), s) );
   for ( auto const& s : ppExportsAllFieldsAvailablefluid )
       ppExportsAllFieldsAvailable.insert( prefixvm( this->fluidModel()->keyword(), s) );
    this->setPostProcessExportsAllFieldsAvailable( ppExportsAllFieldsAvailable );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    super_type::initPostProcess();

    if ( !this->postProcessExportsFields().empty() )
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
std::shared_ptr<std::ostringstream>
HEATFLUID_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << M_heatModel->getInfo()->str();
    *_ostr << M_fluidModel->getInfo()->str();

    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||---------------Info : HeatFluid---------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository();
    *_ostr << "\n   Physical Model"
           << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient" )
           << "\n     -- natural convection  : " << std::string( (M_useNaturalConvection)? "ON":"OFF" );
    //*_ostr << "\n   Numerical Solver"
    //       << "\n     -- solver   : " << M_solverName;
    if ( M_exporter )
    {
        *_ostr << "\n   Exporter"
               << "\n     -- type            : " << M_exporter->type()
               << "\n     -- freq save       : " << M_exporter->freq();
        std::string fieldExported;
        for ( std::string const& fieldName : this->postProcessExportsFields() )
            fieldExported=(fieldExported.empty())? fieldName : fieldExported + " - " + fieldName;
        *_ostr << "\n     -- fields : " << fieldExported;
    }

    *_ostr << "\n   Processors"
           << "\n     -- number of proc : " << this->worldComm().globalSize()
           << "\n     -- current rank : " << this->worldComm().globalRank();

    if ( M_algebraicFactory )
        *_ostr << M_algebraicFactory->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    this->updateTimeStepCurrentResidual();
    this->heatModel()->startTimeStep();
    this->fluidModel()->startTimeStep();
    this->updateTime( this->fluidModel()->time() );
    this->updateParameterValues();
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->updateTimeStepCurrentResidual();
    this->heatModel()->updateTimeStep();
    this->fluidModel()->updateTimeStep();
    this->updateTime( this->fluidModel()->time() );
    this->updateParameterValues();
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateTimeStepCurrentResidual()
{
    if ( M_heatModel->isStationary() && M_fluidModel->isStationaryModel() )
        return;
    if ( !M_algebraicFactory )
        return;
    if ( ( M_fluidModel->timeStepping() == "Theta" ) || ( M_heatModel->timeStepping() == "Theta" ) )
    {
        M_timeStepThetaSchemePreviousContrib->zero();
        M_blockVectorSolution.updateVectorFromSubVectors();
        ModelAlgebraic::DataUpdateResidual dataResidual( M_blockVectorSolution.vectorMonolithic(), M_timeStepThetaSchemePreviousContrib, true, false );
        dataResidual.addInfo( prefixvm( this->prefix(), "time-stepping.evaluate-residual-without-time-derivative" ) );
        if ( M_fluidModel->timeStepping() == "Theta" )
            dataResidual.addInfo( prefixvm( M_fluidModel->prefix(), "time-stepping.evaluate-residual-without-time-derivative" ) );
        if ( M_heatModel->timeStepping() == "Theta" )
            dataResidual.addInfo( prefixvm( M_heatModel->prefix(), "time-stepping.evaluate-residual-without-time-derivative" ) );

        M_algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", false );
        M_algebraicFactory->evaluateResidual( dataResidual );
        M_algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", true );

        if ( M_fluidModel->stabilizationGLS() || M_heatModel->stabilizationGLS() )
        {
            auto & dataInfos = M_algebraicFactory->dataInfos();
            *dataInfos.vectorInfo( prefixvm( this->prefix(),"time-stepping.previous-solution") ) = *M_blockVectorSolution.vectorMonolithic();
        }
    }
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("HeatFluid","exportResults", "start");
    this->timerTool("PostProcessing").start();

    auto mfields = this->modelFields();
    auto symbolExpr = this->symbolsExpr( mfields );
    M_heatModel->exportResults( time, symbolExpr );
    M_fluidModel->exportResults( time, symbolExpr );

    //auto fields = hana::concat( M_heatModel->allFields( M_heatModel->keyword() ), M_fluidModel->allFields( M_fluidModel->keyword() ) );
    auto exprExport = hana::concat( M_materialsProperties->exprPostProcessExports( this->mesh(),this->physicsAvailable(),symbolExpr ),
                                    hana::concat( M_heatModel->exprPostProcessExportsToolbox( symbolExpr,M_heatModel->keyword() ),
                                                  M_fluidModel->exprPostProcessExports( symbolExpr,M_fluidModel->keyword() ) ) );
    this->executePostProcessExports( M_exporter, time, mfields, symbolExpr, exprExport );

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
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );

    this->setParameterValues( paramValues );

    // M_heatModel->updateParameterValues();
    // M_fluidModel->updateParameterValues();
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }
    M_heatModel->setParameterValues( paramValues );
    M_fluidModel->setParameterValues( paramValues );
    //M_fluidModel->updateParameterValues();
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
#if 0
        M_heatModel->setFieldVelocityConvectionIsUsed( true );
        for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
        {
            auto const& range = rangeData.second;
            M_heatModel->updateFieldVelocityConvection( range,idv(M_fluidModel->fieldVelocity()) );
        }
#endif
        M_heatModel->solve();
    }
    else
    {
        if ( M_useSemiImplicitTimeScheme )
        {
#if 0
            M_heatModel->setFieldVelocityConvectionIsUsed( true );
            auto const& uConvection = *M_fluidModel->fieldConvectionVelocityExtrapolatedPtr();
            for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
            {
                auto const& range = rangeData.second;
                M_heatModel->updateFieldVelocityConvection( range,idv(uConvection) );
            }
#endif
            M_heatModel->solve();
            M_fluidModel->solve();
        }
        else
        {
            M_fluidModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("fluid") );
            M_heatModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("heat") );

            M_blockVectorSolution.updateVectorFromSubVectors();

            M_algebraicFactory->solve( "Newton", M_blockVectorSolution.vectorMonolithic() );

            M_blockVectorSolution.localize();
        }
    }

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
HEATFLUID_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner( DataUpdateLinear & data ) const
{
    M_heatModel->updateInHousePreconditioner( data );
    M_fluidModel->updateInHousePreconditioner( data );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner( DataUpdateJacobian & data ) const
{
    M_heatModel->updateInHousePreconditioner( data );
    M_fluidModel->updateInHousePreconditioner( data );
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::postSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    M_heatModel->postSolveNewton( rhs, sol );
    M_fluidModel->postSolveNewton( rhs, sol );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::postSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    M_heatModel->postSolvePicard( rhs, sol );
    M_fluidModel->postSolvePicard( rhs, sol );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::postSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    M_heatModel->postSolveLinear( rhs, sol );
    M_fluidModel->postSolveLinear( rhs, sol );
}




} // end namespace FeelModels
} // end namespace Feel
