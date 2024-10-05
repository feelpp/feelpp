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
    M_useNaturalConvection = boption(_prefix=this->prefix(),_name="use-natural-convection");
    M_solverName = soption(_prefix=this->prefix(),_name="solver");

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

    if ( this->modelProperties().jsonData().contains("Meshes") )
        super_type::super_model_meshes_type::setup( this->modelProperties().jsonData().at("Meshes"), {this->keyword()} );
    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    super_type::super_model_meshes_type::modelMesh( this->keyword() ).setFunctionApplyRemesh(
        [this]( typename super_type::super_model_meshes_type::mesh_base_ptrtype mold,
                typename super_type::super_model_meshes_type::mesh_base_ptrtype mnew ) { this->applyRemesh( std::dynamic_pointer_cast<mesh_type>( mold ),
                                                                                                            std::dynamic_pointer_cast<mesh_type>( mnew ) ); }
                                                                                             );

    CHECK( this->mesh() ) << "mesh generation fail";

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
    int indexBlockVelocity = startIndexBlockFluid + this->fluidModel()->startSubBlockSpaceIndex("velocity");
    int indexBlockTemperature = startIndexBlockHeat + this->heatModel()->startSubBlockSpaceIndex("temperature");

    auto blockMatFluid = M_fluidModel->buildBlockMatrixGraph();
    for (int tk1=0;tk1<nBlockFluid ;++tk1 )
        for (int tk2=0;tk2<nBlockFluid ;++tk2 )
            myblockGraph(startIndexBlockFluid+tk1,startIndexBlockFluid+tk2) = blockMatFluid(tk1,tk2);

    myblockGraph(indexBlockTemperature,indexBlockVelocity) = stencil(_test=M_heatModel->spaceTemperature(),
                                                                     _trial=M_fluidModel->functionSpaceVelocity(),
                                                                     _diag_is_nonzero=false,_close=false)->graph();

    myblockGraph(indexBlockVelocity,indexBlockTemperature) = stencil(_test=M_fluidModel->functionSpaceVelocity(),
                                                                     _trial=M_heatModel->spaceTemperature(),
                                                                     _diag_is_nonzero=false,_close=false)->graph();
    if ( M_fluidModel->stabilizationGLS() && ( this->fluidModel()->stabilizationGLSType() == "pspg" || this->fluidModel()->stabilizationGLSType() == "supg-pspg" || this->fluidModel()->stabilizationGLSType() == "gls" ) )
    {
        int indexBlockPressure = startIndexBlockFluid + this->fluidModel()->startSubBlockSpaceIndex("pressure");
        myblockGraph(indexBlockPressure,indexBlockTemperature) = stencil(_test=M_fluidModel->functionSpacePressure(),
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
HEATFLUID_CLASS_TEMPLATE_TYPE::updatePhysics( typename super_physics_type::PhysicsTreeNode & physicsTree, ModelModels const& models )
{
    physicsTree.addChild( M_heatModel, models );
    physicsTree.addChild( M_fluidModel, models );

    physicsTree.updateMaterialSupportFromChildren( "intersect" );
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("HeatFluid","init", "start" );
    this->timerTool("Constructor").start();

    this->initModelProperties();

    M_heatModel = std::make_shared<heat_model_type>(prefixvm(this->prefix(),"heat"), "heat", this->worldCommPtr(),
                                                    this->subPrefix(), this->repository() );
    M_fluidModel = std::make_shared<fluid_model_type>(prefixvm(this->prefix(),"fluid"), "fluid", this->worldCommPtr(),
                                                      this->subPrefix(), this->repository() );

    // physics
    this->initPhysics(
        this->shared_from_this(),
        [this]( typename super_physics_type::PhysicsTree & physicsTree ) {
            physicsTree.updatePhysics( this->shared_from_this(), this->modelProperties().models() );
            //CHECK( M_fluidModel && M_heatModel ) << "aiai";
            physicsTree.updatePhysics( M_heatModel, this->modelProperties().models() );
            physicsTree.updatePhysics( M_fluidModel, this->modelProperties().models() );
        } );

    // physical properties
    if ( !M_materialsProperties )
    {
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    // init heat toolbox
    M_heatModel->setManageParameterValues( false );
    if ( !M_heatModel->modelPropertiesPtr() )
    {
        M_heatModel->setModelProperties( this->modelPropertiesPtr() );
        M_heatModel->setManageParameterValuesOfModelProperties( false );
    }
    M_heatModel->setModelMeshAsShared( this->modelMesh() );
    M_heatModel->setMaterialsProperties( M_materialsProperties );
    M_heatModel->init( false );

    // init fluid toolbox
    M_fluidModel->setManageParameterValues( false );
    if ( !M_fluidModel->modelPropertiesPtr() )
    {
        M_fluidModel->setModelProperties( this->modelPropertiesPtr() );
        M_fluidModel->setManageParameterValuesOfModelProperties( false );
    }
    M_fluidModel->setModelMeshAsShared( this->modelMesh() );
    M_fluidModel->setMaterialsProperties( M_materialsProperties );
    M_fluidModel->init( false );

    if ( !this->isStationary() )
    {
        this->updateTime( this->timeStepBase()->time() );
        this->setTimeInitial( this->timeStepBase()->timeInitial() );
    }
    else
        M_useSemiImplicitTimeScheme = false;

    // defined the convection fluid velocity into heat toolbox
    std::string fluidVelocitySymbolUsed = M_useSemiImplicitTimeScheme? "beta_u" : "U"; // TODO : get this info from fluid toolbox
    std::string velConvExprStr = (boost::format( nDim==2? "{%1%_%2%_0,%1%_%2%_1}:%1%_%2%_0:%1%_%2%_1":"{%1%_%2%_0,%1%_%2%_1,%1%_%2%_2}:%1%_%2%_0:%1%_%2%_1:%1%_%2%_2" )%M_fluidModel->keyword() %fluidVelocitySymbolUsed ).str();
    nl::json j_heatConvectionVelocity = { {"expr",velConvExprStr} };
    for ( auto const& [physicId,physicObj] : this->physicsFromCurrentType() )
    {
        for ( auto subphysic : physicObj->subphysicsFromType( M_heatModel->physicType() ) )
        {
            auto physicHeatObj = std::dynamic_pointer_cast<ModelPhysicHeat<nDim>>(subphysic);
            CHECK( physicHeatObj ) << "invalid cast";
            physicHeatObj->setupConvection( j_heatConvectionVelocity );
        }
    }


    if ( M_useNaturalConvection )
        M_fluidModel->setStabilizationGLSDoAssembly( false );

    // post-process
    this->initPostProcess();

    // update constant parameters
    this->updateParameterValues();

    // update initial conditions
    this->updateInitialConditions( this->symbolsExpr() );

    // solver
    if ( M_solverName == "automatic" )
    {
        M_solverName = "Newton";
    }

    // algebraic model
    this->initAlgebraicModel();

    // mesh adaptation at event after_init
    using mesh_adaptation_type = typename super_type::super_model_meshes_type::mesh_adaptation_type;
    this->template updateMeshAdaptation<mesh_type>( this->keyword(),
                                                    mesh_adaptation_type::createEvent<mesh_adaptation_type::Event::Type::after_init>(),
                                                    this->symbolsExpr() );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
        this->initAlgebraicFactory();


    this->setIsUpdatedForUse( true );

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("HeatFluid","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::applyRemesh( mesh_ptrtype oldMesh, mesh_ptrtype newMesh, std::shared_ptr<RemeshInterpolation> remeshInterp )
{
    this->log("HeatFluid","applyRemesh", "start");
    //mesh_ptrtype oldMesh = this->mesh();

    // material prop
    this->materialsProperties()->removeMesh( oldMesh );
    this->materialsProperties()->addMesh( newMesh );

    // model mesh
    super_type::super_model_meshes_type::applyRemesh( this->keyword(), newMesh );

    // subtoolboxes
    M_heatModel->applyRemesh( oldMesh, newMesh, remeshInterp );
    M_fluidModel->applyRemesh( oldMesh, newMesh, remeshInterp );

    // reset algebraic data/tools
    this->removeAllAlgebraicDataAndTools();
    this->initAlgebraicModel();

    this->initAlgebraicFactory(); // TODO : Theta time scheme
    this->log("HeatFluid","applyRemesh", "finish");
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
#if 0
        std::string geoExportType="static";//change_coords_only, change, static
#else
        bool useStaticExporter = boption(_name="exporter.use-static-mesh",_prefix=this->prefix());
        std::string geoExportType = useStaticExporter? "static":"change";
#endif
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

    auto se = this->symbolsExpr();
    this->template initPostProcessMeshes<mesh_type>( se );

    // start or restart the export of measures
    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasures().restart( this->timeInitial() );
    }


    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("HeatFluid","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::initAlgebraicModel()
{
    // backend
    this->initAlgebraicBackend();

    // block vector solution
    auto const& blockVectorSolutionFluid = *(M_fluidModel->algebraicBlockVectorSolution());
    auto const& blockVectorSolutionHeat = *(M_heatModel->algebraicBlockVectorSolution());
    int nBlockFluid = blockVectorSolutionFluid.size();
    int nBlockHeat = blockVectorSolutionHeat.size();
    int nBlock = nBlockFluid + nBlockHeat;
    auto bvs = this->initAlgebraicBlockVectorSolution( nBlock );
    int indexBlock=0;
    int numberOfBlockSpaceFluid = 0;
    for ( int k=0;k<nBlockFluid ;++k )
    {
        bvs->operator()(indexBlock+k) = blockVectorSolutionFluid(k);
        numberOfBlockSpaceFluid += blockVectorSolutionFluid(k)->map().numberOfDofIdToContainerId();
    }
    indexBlock += nBlockFluid;
    for ( int k=0;k<nBlockHeat ;++k )
        bvs->operator()(indexBlock+k) = blockVectorSolutionHeat(k);
    indexBlock += nBlockHeat;

    size_type currentStartBlockSpaceIndex = 0;
    this->setStartSubBlockSpaceIndex( "fluid", currentStartBlockSpaceIndex );
    currentStartBlockSpaceIndex += numberOfBlockSpaceFluid;
    this->setStartSubBlockSpaceIndex( "heat", currentStartBlockSpaceIndex );

    // init monolithic vector associated to the block vector
    bvs->buildVector( this->backend() );

    if ( ( M_fluidModel->timeStepping() == "Theta" || M_heatModel->timeStepping() == "Theta" ) &&
         ( M_fluidModel->stabilizationGLS() || M_heatModel->stabilizationGLS() ) )
        M_timeStepThetaSchemePreviousSolution = this->backend()->newVector( bvs->vectorMonolithic()->mapPtr() );
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    if ( M_useNaturalConvection && !M_useSemiImplicitTimeScheme )
    {
        auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
        this->setAlgebraicFactory( algebraicFactory );
        if ( M_fluidModel->hasOperatorPCD() )
            algebraicFactory->preconditionerTool()->attachOperatorPCD( "pcd", M_fluidModel->operatorPCD() );
        if ( M_fluidModel->hasOperatorPMM() )
            algebraicFactory->preconditionerTool()->attachOperatorPMM( "pmm", M_fluidModel->operatorPMM() );

        if ( ( M_fluidModel->timeStepping() == "Theta" ) || ( M_heatModel->timeStepping() == "Theta" ) )
        {
            M_timeStepThetaSchemePreviousContrib = this->backend()->newVector( this->algebraicBlockVectorSolution()->vectorMonolithic()->mapPtr() );
            algebraicFactory->addVectorResidualAssembly( M_timeStepThetaSchemePreviousContrib, 1.0, "Theta-Time-Stepping-Previous-Contrib", true );
            algebraicFactory->addVectorLinearRhsAssembly( M_timeStepThetaSchemePreviousContrib, -1.0, "Theta-Time-Stepping-Previous-Contrib", false );
            if ( M_fluidModel->stabilizationGLS() || M_heatModel->stabilizationGLS() )
            {
                algebraicFactory->dataInfos().addVectorInfo( "time-stepping.previous-solution", M_timeStepThetaSchemePreviousSolution );
            }
        }
    }

    if ( !M_useNaturalConvection || M_useSemiImplicitTimeScheme )
    {
        if ( !M_heatModel->algebraicFactory() )
            M_heatModel->initAlgebraicFactory();
        if ( !M_fluidModel->algebraicFactory() )
            M_fluidModel->initAlgebraicFactory();
        M_heatModel->algebraicFactory()->setFunctionLinearAssembly( std::bind( &self_type::updateLinear_Heat,
                                                                               std::ref( *this ), std::placeholders::_1 ) );
        M_heatModel->algebraicFactory()->setFunctionResidualAssembly( std::bind( &self_type::updateResidual_Heat,
                                                                                 std::ref( *this ), std::placeholders::_1 ) );
        M_heatModel->algebraicFactory()->setFunctionJacobianAssembly( std::bind( &self_type::updateJacobian_Heat,
                                                                                 std::ref( *this ), std::placeholders::_1 ) );

        M_fluidModel->algebraicFactory()->setFunctionLinearAssembly( std::bind( &self_type::updateLinear_Fluid,
                                                                                std::ref( *this ), std::placeholders::_1 ) );
        M_fluidModel->algebraicFactory()->setFunctionResidualAssembly( std::bind( &self_type::updateResidual_Fluid,
                                                                                  std::ref( *this ), std::placeholders::_1 ) );
        M_fluidModel->algebraicFactory()->setFunctionJacobianAssembly( std::bind( &self_type::updateJacobian_Fluid,
                                                                                  std::ref( *this ), std::placeholders::_1 ) );

        M_fluidModel->algebraicFactory()->setFunctionUpdateInHousePreconditionerLinear( std::bind( static_cast<void(self_type::*)(ModelAlgebraic::DataUpdateLinear&) const>(&self_type::updateInHousePreconditioner_Fluid),
                                                                                                   std::ref( *this ), std::placeholders::_1 ) );
        M_fluidModel->algebraicFactory()->setFunctionUpdateInHousePreconditionerJacobian( std::bind( static_cast<void(self_type::*)(ModelAlgebraic::DataUpdateJacobian&) const>(&self_type::updateInHousePreconditioner_Fluid),
                                                                                                     std::ref( *this ), std::placeholders::_1 ) );
    }

}


#if 0
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

    if ( this->algebraicFactory() )
        *_ostr << this->algebraicFactory()->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}
#endif

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );

    // Algebraic Solver
    if ( this->algebraicFactory() )
        this->algebraicFactory()->updateInformationObject( p["Algebraic Solver"] );

    p["Toolbox Heat"] = M_heatModel->journalSection().to_string();
    p["Toolbox Fluid"] = M_fluidModel->journalSection().to_string();
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
HEATFLUID_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfo->add( "Materials Properties", this->materialsProperties()->tabulateInformations(jsonInfo.at("Materials Properties"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    if ( jsonInfo.contains( "Algebraic Solver" ) )
        tabInfo->add( "Algebraic Solver", model_algebraic_factory_type::tabulateInformations( jsonInfo.at("Algebraic Solver"), tabInfoProp ) );

    // generate sub toolboxes info
    if ( M_heatModel && jsonInfo.contains( "Toolbox Heat" ) )
    {
        nl::json::json_pointer jsonPointerHeat( jsonInfo.at( "Toolbox Heat" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerHeat ) )
        {
            auto tabInfos_heat = M_heatModel->tabulateInformations( JournalManager::journalData().at( jsonPointerHeat ), tabInfoProp );
            TabulateInformationsSections::cast( tabInfos_heat )->erase( "Materials Properties" );
            tabInfo->add( "Toolbox Heat", tabInfos_heat );
        }
    }
    if ( M_fluidModel && jsonInfo.contains( "Toolbox Fluid" ) )
    {
        nl::json::json_pointer jsonPointerFluid( jsonInfo.at( "Toolbox Fluid" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerFluid ) )
        {
            auto tabInfos_fluid = M_fluidModel->tabulateInformations( JournalManager::journalData().at( jsonPointerFluid ), tabInfoProp );
            TabulateInformationsSections::cast( tabInfos_fluid )->erase( "Materials Properties" );
            tabInfo->add( "Toolbox Fluid", tabInfos_fluid );
        }
    }

    return tabInfo;
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    this->fluidModel()->startTimeStepPreProcess();
    this->updateTimeStepCurrentResidual();
    this->heatModel()->startTimeStep();
    this->fluidModel()->startTimeStep( false );
    this->updateTime( this->fluidModel()->time() );
    this->updateParameterValues();
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
#if 1 // AT END OF TIME STEP
    using mesh_adaptation_type = typename super_type::super_model_meshes_type::mesh_adaptation_type;
    this->template updateMeshAdaptation<mesh_type>( this->keyword(),
                                                    mesh_adaptation_type::createEvent<mesh_adaptation_type::Event::Type::each_time_step>( this->time(),this->timeStepBase()->iteration() ),
                                                    this->symbolsExpr() );
#endif
    this->updateTimeStepCurrentResidual();
    this->heatModel()->updateTimeStep();
    this->fluidModel()->updateTimeStep();
    this->updateTime( this->fluidModel()->time() );
    this->updateParameterValues();

#if 0 // AT BEGIN OF TIME STEP
    using mesh_adaptation_type = typename super_type::super_model_meshes_type::mesh_adaptation_type;
    this->template updateMeshAdaptation<mesh_type>( this->keyword(),
                                                    mesh_adaptation_type::createEvent<mesh_adaptation_type::Event::Type::each_time_step>( this->time(),this->timeStepBase()->iteration() ),
                                                    this->symbolsExpr() );
#endif

}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateTimeStepCurrentResidual()
{
    if ( M_heatModel->isStationary() && M_fluidModel->isStationaryModel() )
        return;

    auto algebraicFactory = this->algebraicFactory();
    if ( ( M_fluidModel->timeStepping() == "Theta" ) || ( M_heatModel->timeStepping() == "Theta" ) )
    {
        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
        if ( algebraicFactory )
        {
            M_timeStepThetaSchemePreviousContrib->zero();
            ModelAlgebraic::DataUpdateResidual dataResidual( this->algebraicBlockVectorSolution()->vectorMonolithic(), M_timeStepThetaSchemePreviousContrib, true, false );
            dataResidual.addInfo( prefixvm( this->prefix(), "time-stepping.evaluate-residual-without-time-derivative" ) );
            if ( M_fluidModel->timeStepping() == "Theta" )
                dataResidual.addInfo( prefixvm( M_fluidModel->prefix(), "time-stepping.evaluate-residual-without-time-derivative" ) );
            if ( M_heatModel->timeStepping() == "Theta" )
                dataResidual.addInfo( prefixvm( M_heatModel->prefix(), "time-stepping.evaluate-residual-without-time-derivative" ) );

            M_fluidModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("fluid") );
            M_heatModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("heat") );

            algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", false );
            algebraicFactory->evaluateResidual( dataResidual );
            algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", true );
        }

        if ( M_fluidModel->stabilizationGLS() || M_heatModel->stabilizationGLS() )
        {
            *M_timeStepThetaSchemePreviousSolution = *this->algebraicBlockVectorSolution()->vectorMonolithic();

            if ( algebraicFactory )
            {
                auto & dataInfos = algebraicFactory->dataInfos();
                // *dataInfos.vectorInfo( "time-stepping.previous-solution" ) = *M_blockVectorSolution.vectorMonolithic();
                dataInfos.addParameterValuesInfo( "time-stepping.previous-parameter-values", M_currentParameterValues );
            }
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
    auto se = this->symbolsExpr( mfields );
    M_heatModel->exportResults( time, se );
    M_fluidModel->exportResults( time, se );

    //auto fields = hana::concat( M_heatModel->allFields( M_heatModel->keyword() ), M_fluidModel->allFields( M_fluidModel->keyword() ) );
    auto exprExport = hana::concat( M_materialsProperties->exprPostProcessExports( this->mesh(),this->physicsAvailable(),se ),
                                    hana::concat( M_heatModel->exprPostProcessExportsToolbox( se,M_heatModel->keyword() ),
                                                  M_fluidModel->exprPostProcessExports( se,M_fluidModel->keyword() ) ) );
    if ( M_exporter && M_exporter->exporterGeometry() == EXPORTER_GEOMETRY_CHANGE ) // TODO mv this code
        M_exporter->defaultTimeSet()->setMesh( this->mesh() );
    this->executePostProcessExports( M_exporter, time, mfields, se, exprExport );


    auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
    auto defaultRangeMeshElements = mom->isDefinedOnWholeMesh( this->physicsAvailableFromCurrentType() )?
        elements(this->mesh()) : markedelements(this->mesh(), mom->markers( this->physicsAvailableFromCurrentType() ));
    model_measures_quantities_empty_t mquantities;
    // execute common post process and save measures
    super_type::executePostProcessMeasures( time, this->mesh(), defaultRangeMeshElements, se, mfields, mquantities );


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
    for ( auto [physicName,physicData] : this->physics/*FromCurrentType*/() )
        physicData->updateParameterValues( paramValues );

    this->setParameterValues( paramValues );
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto const& [param,val] : paramValues )
        M_currentParameterValues[param] = val;

    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->modelProperties().initialConditions().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        physicData->setParameterValues( paramValues );

    M_heatModel->setParameterValues( paramValues );
    M_fluidModel->setParameterValues( paramValues );
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
        M_heatModel->solve();
        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
    }
    else
    {
        if ( M_useSemiImplicitTimeScheme )
        {
            M_heatModel->solve();
            M_fluidModel->solve();
            this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
        }
        else
        {
            this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();

            M_fluidModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("fluid") );
            M_heatModel->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("heat") );

            this->algebraicFactory()->solve( M_solverName, this->algebraicBlockVectorSolution()->vectorMonolithic() );

            this->algebraicBlockVectorSolution()->localize();
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
    vector_ptrtype const& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( vecCurrentSolution, this->heatModel(), this->fluidModel() );
    M_heatModel->updateInHousePreconditioner( data );
    M_fluidModel->updateInHousePreconditioner( data, mctx );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner( DataUpdateJacobian & data ) const
{
    vector_ptrtype const& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( vecCurrentSolution, this->heatModel(), this->fluidModel() );
    M_heatModel->updateInHousePreconditioner( data );
    M_fluidModel->updateInHousePreconditioner( data, mctx );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner_Fluid( DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( this->heatModel()->algebraicBlockVectorSolution()->vectorMonolithic(), 0,
                                    vecCurrentSolution, this->fluidModel()->startBlockSpaceIndexVector() );
    M_fluidModel->updateInHousePreconditioner( data, mctx );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner_Fluid( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( this->heatModel()->algebraicBlockVectorSolution()->vectorMonolithic(), 0,
                                    vecCurrentSolution, this->fluidModel()->startBlockSpaceIndexVector() );
    M_fluidModel->updateInHousePreconditioner( data, mctx );
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
