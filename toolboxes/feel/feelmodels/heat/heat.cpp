/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/heat/heat.hpp>

#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>

namespace Feel
{
namespace FeelModels
{

HEAT_CLASS_TEMPLATE_DECLARATIONS
HEAT_CLASS_TEMPLATE_TYPE::Heat( std::string const& prefix,
                                std::string const& keyword,
                                worldcomm_ptr_t const& worldComm,
                                std::string const& subPrefix,
                                ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelPhysics<nDim>( "heat" ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep )
{
    this->log("Heat","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".HeatTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    this->log("Heat","constructor", "finish");

}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_useExtendedDoftable = boption(_name="use-extended-doftable",_prefix=this->prefix());

    M_stabilizationGLS = boption(_name="stabilization-gls",_prefix=this->prefix());
    M_stabilizationGLSType = soption(_name="stabilization-gls.type",_prefix=this->prefix());

    M_stabilizationGLS_checkConductivityDependencyOnCoordinates = boption(_name="stabilization-gls.check-conductivity-dependency-on-coordinates",_prefix=this->prefix());

    // time stepping
    M_timeStepping = soption(_name="time-stepping",_prefix=this->prefix());
    M_timeStepThetaValue = doption(_name="time-stepping.theta.value",_prefix=this->prefix());

    M_solverName = soption(_name="solver",_prefix=this->prefix());
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("Heat","initMesh", "start");
    this->timerTool("Constructor").start();

    if ( this->modelProperties().jsonData().contains("Meshes") )
        super_type::super_model_meshes_type::setup( this->modelProperties().jsonData().at("Meshes"), {this->keyword()} );
    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( this->mesh() ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("initMesh");
    this->log("Heat","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initMaterialProperties()
{
    this->log("Heat","initMaterialProperties", "start" );
    this->timerTool("Constructor").start();

    if ( !M_materialsProperties )
    {
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    double tElpased = this->timerTool("Constructor").stop("initMaterialProperties");
    this->log("Heat","initMaterialProperties",(boost::format("finish in %1% s")%tElpased).str() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    this->log("Heat","initFunctionSpaces", "start" );
    this->timerTool("Constructor").start();

    auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
    // functionspace
    if ( mom->isDefinedOnWholeMesh( this->physicsAvailableFromCurrentType() ) )
    {
        M_rangeMeshElements = elements(this->mesh());
        M_Xh = space_temperature_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(), _extended_doftable=M_useExtendedDoftable );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), mom->markers( this->physicsAvailableFromCurrentType() ));
        M_Xh = space_temperature_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),_range=M_rangeMeshElements, _extended_doftable=M_useExtendedDoftable );
    }

    //M_fieldTemperature.reset( new element_temperature_type(M_Xh,"temperature"));
    M_fieldTemperature =  M_Xh->elementPtr( "temperature" );

    double tElpased = this->timerTool("Constructor").stop("initFunctionSpaces");
    this->log("Heat","initFunctionSpaces",(boost::format("finish in %1% s")%tElpased).str() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
HEAT_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    myblockGraph(0,0) = stencil(_test=this->spaceTemperature(),
                                _trial=this->spaceTemperature() )->graph();
    return myblockGraph;
}
#if 0
HEAT_CLASS_TEMPLATE_DECLARATIONS
typename HEAT_CLASS_TEMPLATE_TYPE::size_type
HEAT_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    size_type res = this->spaceTemperature()->nLocalDofWithGhost();
    return res;
}
#endif

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("Heat","init", "start" );
    this->timerTool("Constructor").start();

    this->initModelProperties();

    // physics
    this->initPhysics( this->shared_from_this(), this->modelProperties().models() );

    this->initMaterialProperties();

    this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    for ( auto & [physicId,physicObj] : this->physicsFromCurrentType() )
        std::static_pointer_cast<ModelPhysicHeat<nDim>>(physicObj)->updateForUse( this->materialsProperties(), this->mesh() );

    this->initFunctionSpaces();

    this->initBoundaryConditions();

    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();

    // stabilization gls
    if ( M_stabilizationGLS )
    {
        typedef StabilizationGLSParameter<mesh_type, nOrderTemperature> stab_gls_parameter_impl_type;
        M_stabilizationGLSParameter.reset( new stab_gls_parameter_impl_type( this->mesh(),prefixvm(this->prefix(),"stabilization-gls.parameter") ) );
        M_stabilizationGLSParameter->init();
    }

    // update constant parameters into
    this->updateParameterValues();

    // update initial conditions
    this->updateInitialConditions( this->symbolsExpr() );

    // post-process
    this->initPostProcess();

    // automatic solver selection
    if ( M_solverName == "automatic" )
    {
        auto mfields = this->modelFields();
        auto se = this->symbolsExpr( mfields );
        auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( 0/*rowStartInVector*/ ) );
        auto trialSymbolNames = tse.names();
        bool isNonLinear = false;
        for ( std::string tsName : trialSymbolNames )
        {
            if ( this->materialsProperties()->hasThermalConductivityDependingOnSymbol( tsName ) )
            {
                isNonLinear = true;
                break;
            }
            for ( auto const& [bcName,bcData] : M_boundaryConditions->heatFlux() )
            {
                auto neumannExpr = bcData->expr();
                if ( neumannExpr.hasSymbolDependency( tsName, se ) )
                {
                    isNonLinear = true;
                    break;
                }
            }

            if ( isNonLinear )
                break;
        }
        M_solverName = isNonLinear? "Newton" : "Linear";
    }


    this->initAlgebraicModel();

    // algebraic solver
    if ( buildModelAlgebraicFactory )
        this->initAlgebraicFactory();

    this->setIsUpdatedForUse( true );

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("Heat","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initAlgebraicModel()
{
    // backend
    this->initAlgebraicBackend();

    // subspaces index
    size_type currentStartIndex = 0;
    this->setStartSubBlockSpaceIndex( "temperature", currentStartIndex++ );

    this->updateAlgebraicDofEliminationIds();

     // vector solution
    auto bvs = this->initAlgebraicBlockVectorSolution( 1 );
    bvs->operator()(0) = this->fieldTemperaturePtr();
    // init petsc vector associated to the block
    bvs->buildVector( this->backend() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::applyRemesh( mesh_ptrtype const& newMesh )
{
    mesh_ptrtype oldMesh = this->mesh();

    // material prop
    this->materialsProperties()->removeMesh( oldMesh );
    this->materialsProperties()->addMesh( newMesh );

    this->setMesh( newMesh );

    // function space and fields
    space_temperature_ptrtype old_Xh = M_Xh;
    element_temperature_ptrtype old_fieldTemperature = M_fieldTemperature;
    this->initFunctionSpaces();

    // createInterpolationOp
    auto opI_temperature = opInterpolation(_domainSpace=old_Xh,
                                           _imageSpace=M_Xh,
                                           _range=M_rangeMeshElements );
    auto matrixInterpolation_temperature = opI_temperature->matPtr();
    matrixInterpolation_temperature->multVector( *old_fieldTemperature, *M_fieldTemperature );

    // time stepping
    if ( M_bdfTemperature )
        M_bdfTemperature->applyRemesh( M_Xh, matrixInterpolation_temperature );

    // TODO : stabilization gls

    // TODO : post process ??

    // reset algebraic data/tools
    this->removeAllAlgebraicDataAndTools();
    this->initAlgebraicModel();

    this->initAlgebraicFactory(); // TODO : Theta time scheme

}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    this->log("Heat","initTimeStep", "start" );
    this->timerTool("Constructor").start();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix

    int bdfOrder = 1;
    if ( M_timeStepping == "BDF" )
        bdfOrder = ioption(_prefix=this->prefix(),_name="bdf.order");
    int nConsecutiveSave = std::max( 3, bdfOrder ); // at least 3 is required when restart with theta scheme

    M_bdfTemperature = this->createBdf( this->spaceTemperature(),"temperature", bdfOrder, nConsecutiveSave, myFileFormat );

    if (!this->doRestart())
    {
        // up current time
        this->updateTime( M_bdfTemperature->timeInitial() );
    }
    else
    {
        // start time step
        double tir = M_bdfTemperature->restart();
        // load a previous solution as current solution
        *this->fieldTemperaturePtr() = M_bdfTemperature->unknown(0);
        // up initial time
        this->setTimeInitial( tir );
        // up current time
        this->updateTime( tir );
    }

    double tElapsed = this->timerTool("Constructor").stop("initTimeStep");
    this->log("Heat","initTimeStep", (boost::format("finish in %1% s") %tElapsed).str() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("Heat","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->setPostProcessExportsAllFieldsAvailable( {"temperature","velocity-convection"} );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    this->setPostProcessSaveAllFieldsAvailable( {"temperature" } );
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

        // restart exporter
        if ( M_exporter->doExport() && this->doRestart() && this->restartPath().empty() )
            M_exporter->restart(this->timeInitial());
    }


    if ( this->modelProperties().postProcess().hasJsonProperties( this->keyword() ) )
    {
        auto const& j_pp = this->modelProperties().postProcess().jsonProperties( this->keyword() );
        std::string ppTypeMeasures = "Measures";
        if ( j_pp.contains( ppTypeMeasures ) )
        {
            auto j_pp_measures = j_pp.at( ppTypeMeasures );
            for ( auto const& [j_pp_measureskey,j_pp_measuresval] : j_pp_measures.items() )
            {
                if ( j_pp_measureskey == "Normal-Heat-Flux" )
                {
                    for ( auto const& [j_pp_measures_nhfkey,j_pp_measures_nhfval] : j_pp_measuresval.items() )
                    {
                        auto indexesAllCases = ModelIndexes::generateAllCases( j_pp_measures_nhfval );
                        for ( auto const& indexes : indexesAllCases )
                        {
                            ModelMeasuresNormalFluxGeneric ppFlux;
                            ppFlux.setup( j_pp_measures_nhfval, indexes.replace( j_pp_measures_nhfkey ), indexes );
                            if ( !ppFlux.markers().empty() )
                                M_postProcessMeasuresNormalHeatFlux[ppFlux.name()] = ppFlux;
                        }
                    }
                }
            }
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

    double tElpased = this->timerTool("Constructor").stop("initPostProcess");
    this->log("Heat","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
    this->setAlgebraicFactory( algebraicFactory );

    if ( M_timeStepping == "Theta" )
    {
        M_timeStepThetaSchemePreviousContrib = this->backend()->newVector( this->algebraicBlockVectorSolution()->vectorMonolithic()->mapPtr() );
        algebraicFactory->addVectorResidualAssembly( M_timeStepThetaSchemePreviousContrib, 1.0, "Theta-Time-Stepping-Previous-Contrib", true );
        algebraicFactory->addVectorLinearRhsAssembly( M_timeStepThetaSchemePreviousContrib, -1.0, "Theta-Time-Stepping-Previous-Contrib", false );
        if ( M_stabilizationGLS )
            algebraicFactory->dataInfos().addVectorInfo( "time-stepping.previous-solution", this->backend()->newVector( this->algebraicBlockVectorSolution()->vectorMonolithic()->mapPtr() ) );
    }

}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    super_physics_type::updateInformationObjectFromCurrentType( p["Physics"] );

    // Physics
    nl::json subPt;
    subPt.emplace( "time mode", std::string( (this->isStationary())?"Stationary":"Transient") );
    p["Physics2"] = subPt;

    // Boundary Conditions
    M_boundaryConditions->updateInformationObject( p["Boundary Conditions"] );

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );

    // FunctionSpace
    subPt.clear();
    subPt["Temperature"] = M_Xh->journalSection().to_string();
    p.emplace( "Function Spaces",  subPt );
    if ( M_stabilizationGLS )
    {
        subPt.clear();
        subPt.emplace( "type", M_stabilizationGLSType );
        if ( M_stabilizationGLSParameter )
            subPt.emplace( "paramter method", M_stabilizationGLSParameter->method() );
        p["Finite element stabilization"] = subPt;
    }

    this->modelFields().updateInformationObject( p["Fields"] );

    if ( !this->isStationary() )
    {
        subPt.clear();
        subPt.emplace( "initial time", this->timeStepBase()->timeInitial() );
        subPt.emplace( "final time", this->timeStepBase()->timeFinal() );
        subPt.emplace( "time step", this->timeStepBase()->timeStep() );
        subPt.emplace( "type", M_timeStepping );
        p["Time Discretization"] = subPt;
    }


    // Algebraic Solver
    if ( this->algebraicFactory() )
    {
        this->algebraicFactory()->updateInformationObject( p["Algebraic Solver"] );
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
HEAT_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics") )
        tabInfo->add( "Physics", super_physics_type::tabulateInformations( jsonInfo.at("Physics"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics2") )
    {
        Feel::Table tabInfoPhysics;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoPhysics, jsonInfo.at("Physics2"), tabInfoProp );
        tabInfo->add( "Physics2", TabulateInformations::New( tabInfoPhysics, tabInfoProp ) );
    }

    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfo->add( "Materials Properties", this->materialsProperties()->tabulateInformations(jsonInfo.at("Materials Properties"), tabInfoProp ) );

    if ( jsonInfo.contains("Boundary Conditions") )
        tabInfo->add( "Boundary Conditions", HeatBoundaryConditions::tabulateInformations( jsonInfo.at("Boundary Conditions"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    if ( jsonInfo.contains("Function Spaces") )
    {
        auto const& jsonInfoFunctionSpaces = jsonInfo.at("Function Spaces");
        auto tabInfoFunctionSpaces = TabulateInformationsSections::New( tabInfoProp );

        nl::json::json_pointer jsonPointerSpaceTemperature( jsonInfoFunctionSpaces.at( "Temperature" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerSpaceTemperature ) )
            tabInfoFunctionSpaces->add( "Temperature", TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( JournalManager::journalData().at( jsonPointerSpaceTemperature ), tabInfoProp ) );

        tabInfo->add( "Function Spaces", tabInfoFunctionSpaces );
    }

    // fields
    if ( jsonInfo.contains("Fields") )
        tabInfo->add( "Fields", TabulateInformationTools::FromJSON::tabulateInformationsModelFields( jsonInfo.at("Fields"), tabInfoProp ) );

    if ( jsonInfo.contains("Time Discretization") )
    {
        Feel::Table tabInfoTimeDiscr;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoTimeDiscr, jsonInfo.at("Time Discretization"), tabInfoProp );
        tabInfo->add( "Time Discretization", TabulateInformations::New( tabInfoTimeDiscr, tabInfoProp ) );
    }

    if ( jsonInfo.contains("Finite element stabilization") )
    {
        Feel::Table tabInfoStab;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoStab, jsonInfo.at("Finite element stabilization"), tabInfoProp );
        tabInfo->add( "Finite element stabilization", TabulateInformations::New( tabInfoStab, tabInfoProp ) );
    }

    if ( jsonInfo.contains( "Algebraic Solver" ) )
        tabInfo->add( "Algebraic Solver", model_algebraic_factory_type::tabulateInformations( jsonInfo.at("Algebraic Solver"), tabInfoProp ) );

    return tabInfo;
}

#if 0
HEAT_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
HEAT_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
#if 0
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||-----------------Info : Heat------------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository();
    *_ostr << "\n   Physical Model"
           << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient");
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
            if ( this->hasVelocityConvectionExpr( matName ) )
                *_ostr << "\n     -- convection-velocity [" << matName << "] : " <<  str( this->velocityConvectionExpr( matName ).expression() );
    *_ostr << "\n   Boundary conditions"
           << M_bcDirichletMarkerManagement.getInfoDirichletBC()
           << M_bcNeumannMarkerManagement.getInfoNeumannBC()
           << M_bcRobinMarkerManagement.getInfoRobinBC();
    *_ostr << this->materialsProperties()->getInfoMaterialParameters()->str();
#if 0
    *_ostr << "\n   Mesh Discretization"
           << "\n     -- mesh filename      : " << this->meshFile()
           << "\n     -- number of element : " << this->mesh()->numGlobalElements()
           << "\n     -- order             : " << nOrderGeo;
#endif
    *_ostr << "\n   Space Temperature Discretization"
           << "\n     -- order         : " << nOrderPoly
           << "\n     -- number of dof : " << M_Xh->nDof() << " (" << M_Xh->nLocalDof() << ")";
    if ( !this->isStationary() )
    {
        *_ostr << "\n   Time Discretization"
               << "\n     -- initial time : " << this->timeStepBase()->timeInitial()
               << "\n     -- final time   : " << this->timeStepBase()->timeFinal()
               << "\n     -- time step    : " << this->timeStepBase()->timeStep()
               << "\n     -- type : " << M_timeStepping;
    }
    if ( M_stabilizationGLS )
    {
        *_ostr << "\n   Finite element stabilization"
               << "\n     -- type : " << M_stabilizationGLSType;
        if ( M_stabilizationGLSParameter )
            *_ostr << "\n     -- paramter method : " << M_stabilizationGLSParameter->method();
    }
    if ( this->algebraicFactory() )
        *_ostr << this->algebraicFactory()->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";
#endif
    return _ostr;
}
#endif

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );
    for ( auto [physicName,physicData] : this->physics/*FromCurrentType*/() )
        physicData->updateParameterValues( paramValues );

    this->updateParameterValues_postProcess( paramValues, prefixvm("postprocess",this->keyword(),"_" ) );

    this->setParameterValues( paramValues );
}
HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    this->log("Heat","setParameterValues", "start");

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

    M_boundaryConditions->setParameterValues( paramValues );

    this->log("Heat","setParameterValues", "finish");
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    M_boundaryConditions = std::make_shared<boundary_conditions_type>( this->shared_from_this() );
    if ( !this->modelProperties().boundaryConditions().hasSection( this->keyword() ) )
        return;
    M_boundaryConditions->setup( this->modelProperties().boundaryConditions().section( this->keyword() ) );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateAlgebraicDofEliminationIds()
{
    for ( auto const& [bcName,bcData] : M_boundaryConditions->temperatureImposed() )
        bcData->updateDofEliminationIds( *this, "temperature", this->spaceTemperature() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Heat","solve", "start");
    this->timerTool("Solve").start();

    this->setStartBlockSpaceIndex( 0 );

    this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();

    this->algebraicFactory()->solve( M_solverName, this->algebraicBlockVectorSolution()->vectorMonolithic() );

    this->algebraicBlockVectorSolution()->localize();

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("Heat","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}



HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    this->exportResults( time, mfields, se, this->exprPostProcessExports( se ) );
}
#if 0
HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::executePostProcessMeasures( double time )
{
    auto mfields = this->modelFields();
    this->executePostProcessMeasures( time, mfields, this->symbolsExpr( mfields ) );
}
#endif

HEAT_CLASS_TEMPLATE_DECLARATIONS
bool
HEAT_CLASS_TEMPLATE_TYPE::checkResults() const
{
    const_cast<self_type*>(this)->updateParameterValues();
    auto se = this->symbolsExpr();
    return super_type::checkResults( se );
}


HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    this->log("Heat","startTimeStep", "start");

    // some time stepping require to compute residual without time derivative
    this->updateTimeStepCurrentResidual();

    // start time step
    if (!this->doRestart())
        M_bdfTemperature->start( M_bdfTemperature->unknowns() );
     // up current time
    this->updateTime( M_bdfTemperature->time() );

    // update all expressions in bc or in house prec
    this->updateParameterValues();

    this->log("Heat","startTimeStep", "finish");
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->log("Heat","updateTimeStep", "start");
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    // some time stepping require to compute residual without time derivative
    this->updateTimeStepCurrentResidual();

    bool rebuildCstAssembly = false;
    if ( M_timeStepping == "BDF" )
    {
        int previousTimeOrder = this->timeStepBdfTemperature()->timeOrder();
        M_bdfTemperature->next( this->fieldTemperature() );
        int currentTimeOrder = this->timeStepBdfTemperature()->timeOrder();
        rebuildCstAssembly = previousTimeOrder != currentTimeOrder && this->timeStepBase()->strategy() == TS_STRATEGY_DT_CONSTANT;
        this->updateTime( this->timeStepBdfTemperature()->time() );
    }
    else if ( M_timeStepping == "Theta" )
    {
        M_bdfTemperature->next( this->fieldTemperature() );
        this->updateTime( this->timeStepBdfTemperature()->time() );
    }

    if ( rebuildCstAssembly )
        this->setNeedToRebuildCstPart(true);

    this->updateParameterValues();

    this->timerTool("TimeStepping").stop("updateTimeStep");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("Heat","updateTimeStep", "finish");
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateTimeStepCurrentResidual()
{
    if ( this->isStationary() )
        return;

    auto algebraicFactory = this->algebraicFactory();
    if ( !algebraicFactory )
        return;

    if ( M_timeStepping == "Theta" )
    {
        M_timeStepThetaSchemePreviousContrib->zero();
        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
        ModelAlgebraic::DataUpdateResidual dataResidual( this->algebraicBlockVectorSolution()->vectorMonolithic(), M_timeStepThetaSchemePreviousContrib, true, false );
        dataResidual.addInfo( prefixvm( this->prefix(), "time-stepping.evaluate-residual-without-time-derivative" ) );
        this->setStartBlockSpaceIndex( 0 );
        algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", false );
        algebraicFactory->evaluateResidual( dataResidual );
        algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", true );

        if ( M_stabilizationGLS )
        {
            auto & dataInfos = algebraicFactory->dataInfos();
            *dataInfos.vectorInfo( "time-stepping.previous-solution" ) = *this->algebraicBlockVectorSolution()->vectorMonolithic();
            dataInfos.addParameterValuesInfo( "time-stepping.previous-parameter-values", M_currentParameterValues );
        }
    }
}


} // end namespace FeelModels
} // end namespace Feel
