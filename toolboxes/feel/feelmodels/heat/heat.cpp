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
    // M_fieldVelocityConvectionIsUsed = boption(_name="use_velocity-convection",_prefix=this->prefix()) ||
    //     Environment::vm().count(prefixvm(this->prefix(),"velocity-convection").c_str());
    // M_fieldVelocityConvectionIsIncompressible = boption(_name="velocity-convection_is_incompressible",_prefix=this->prefix());

    M_stabilizationGLS = boption(_name="stabilization-gls",_prefix=this->prefix());
    M_stabilizationGLSType = soption(_name="stabilization-gls.type",_prefix=this->prefix());

    // time stepping
    M_timeStepping = soption(_name="time-stepping",_prefix=this->prefix());
    M_timeStepThetaValue = doption(_name="time-stepping.theta.value",_prefix=this->prefix());
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("Heat","initMesh", "start");
    this->timerTool("Constructor").start();

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

    if ( Environment::vm().count(prefixvm(this->prefix(),"velocity-convection").c_str()) )
    {
        auto theExpr = expr<nDim,1>( soption(_prefix=this->prefix(),_name="velocity-convection"),
                                     "",this->worldComm(),this->repository().expr() );
        for ( std::string const& matName : M_materialsProperties->physicToMaterials( this->physicsAvailableFromCurrentType() ) )
            this->setVelocityConvectionExpr( matName,theExpr );
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
        M_Xh = space_temperature_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), mom->markers( this->physicsAvailableFromCurrentType() ));
        M_Xh = space_temperature_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
    }

    M_fieldTemperature.reset( new element_temperature_type(M_Xh,"temperature"));

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

HEAT_CLASS_TEMPLATE_DECLARATIONS
typename HEAT_CLASS_TEMPLATE_TYPE::size_type
HEAT_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    size_type res = this->spaceTemperature()->nLocalDofWithGhost();
    return res;
}


HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("Heat","init", "start" );
    this->timerTool("Constructor").start();

    if ( this->physics().empty() )
        this->initPhysics( this->keyword(), this->modelProperties().models() );

    this->initMaterialProperties();

    if ( !this->mesh() )
        this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

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

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh->worldCommPtr() );

    // subspaces index
    size_type currentStartIndex = 0;
    this->setStartSubBlockSpaceIndex( "temperature", currentStartIndex++ );

     // vector solution
    M_blockVectorSolution.resize( 1 );
    M_blockVectorSolution(0) = this->fieldTemperaturePtr();

    // init petsc vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );

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
        std::string geoExportType="static";//change_coords_only, change, static
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export",
                               _geo=geoExportType,
                               _path=this->exporterPath() );

        // restart exporter
        if ( M_exporter->doExport() && this->doRestart() && this->restartPath().empty() )
            M_exporter->restart(this->timeInitial());
    }

    pt::ptree ptree = this->modelProperties().postProcess().pTree( this->keyword() );
    //  heat flux measures
    std::string ppTypeMeasures = "Measures";
    for( auto const& ptreeLevel0 : ptree )
    {
        std::string ptreeLevel0Name = ptreeLevel0.first;
        if ( ptreeLevel0Name != ppTypeMeasures ) continue;
        for( auto const& ptreeLevel1 : ptreeLevel0.second )
        {
            std::string ptreeLevel1Name = ptreeLevel1.first;

            if ( ptreeLevel1Name == "Normal-Heat-Flux" )
            {
                for( auto const& ptreeNormalHeatFlux : ptreeLevel1.second )
                {
                    auto indexesAllCases = ModelIndexes::generateAllCases( ptreeNormalHeatFlux.second );
                    for ( auto const& indexes : indexesAllCases )
                    {
                        ModelMeasuresNormalFluxGeneric ppFlux;
                        ppFlux.setup( ptreeNormalHeatFlux.second, indexes.replace( ptreeNormalHeatFlux.first ), indexes );
                        if ( !ppFlux.markers().empty() )
                            M_postProcessMeasuresNormalHeatFlux[ppFlux.name()] = ppFlux;
                    }
                }
            }
        }
    }

    // point measures
    auto fieldNamesWithSpaceTemperature = std::make_pair( std::set<std::string>({"temperature"}), this->spaceTemperature() );
    auto fieldNamesWithSpaces = hana::make_tuple( fieldNamesWithSpaceTemperature );
    M_measurePointsEvaluation = std::make_shared<measure_points_evaluation_type>( fieldNamesWithSpaces );
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( this->keyword() ) )
        M_measurePointsEvaluation->init( evalPoints );

    // start or restart the export of measures
    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasuresIO().restart( "time", this->timeInitial() );
        else
            this->postProcessMeasuresIO().setMeasure( "time", this->timeInitial() ); //just for have time in first column
    }

    double tElpased = this->timerTool("Constructor").stop("initPostProcess");
    this->log("Heat","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );

    if ( M_timeStepping == "Theta" )
    {
        M_timeStepThetaSchemePreviousContrib = this->backend()->newVector(M_blockVectorSolution.vectorMonolithic()->mapPtr() );
        M_algebraicFactory->addVectorResidualAssembly( M_timeStepThetaSchemePreviousContrib, 1.0, "Theta-Time-Stepping-Previous-Contrib", true );
        M_algebraicFactory->addVectorLinearRhsAssembly( M_timeStepThetaSchemePreviousContrib, -1.0, "Theta-Time-Stepping-Previous-Contrib", false );
        if ( M_stabilizationGLS )
            M_algebraicFactory->dataInfos().addVectorInfo( "time-stepping.previous-solution", this->backend()->newVector( M_blockVectorSolution.vectorMonolithic()->mapPtr() ) );
    }

}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateInformationObject( pt::ptree & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.get_child_optional( "Environment" ) )
        return;

    pt::ptree subPt;
    super_type::super_model_base_type::updateInformationObject( subPt );
    p.put_child( "Environment", subPt );
    subPt.clear();
    super_type::super_model_meshes_type::updateInformationObject( subPt );
    p.put_child( "Meshes", subPt );

    // Physics
    pt::ptree subPt2;
    subPt.clear();
    subPt.put( "time mode", std::string( (this->isStationary())?"Stationary":"Transient") );
    //subPt.put( "velocity-convection",  std::string( (this->fieldVelocityConvectionIsUsedAndOperational())?"Yes":"No" ) );
    p.put_child( "Physics", subPt );

    // Boundary Conditions
    subPt.clear();
    subPt2.clear();
    M_bcDirichletMarkerManagement.updateInformationObjectDirichletBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    subPt2.clear();
    M_bcNeumannMarkerManagement.updateInformationObjectNeumannBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    subPt2.clear();
    M_bcRobinMarkerManagement.updateInformationObjectRobinBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    p.put_child( "Boundary Conditions",subPt );

    // Materials properties
    if ( this->materialsProperties() )
    {
        subPt.clear();
        this->materialsProperties()->updateInformationObject( subPt );
        p.put_child( "Materials Properties", subPt );
    }

    // FunctionSpace
    subPt.clear();
    subPt2.clear();
    M_Xh->updateInformationObject( subPt2 );
    //subPt.put_child( "FunctionSpace Temperature",  M_Xh->journalSectionName() );
    subPt.put_child( "Temperature", subPt2 );
    p.put_child( "Function Spaces",  subPt );
    //if ( this->fieldVelocityConvectionIsUsedAndOperational() )
    //    p.put( "FunctionSpace Velocity Convection", M_XhVelocityConvection->journalSectionName() );
    if ( M_stabilizationGLS )
    {
        subPt.clear();
        subPt.put( "type", M_stabilizationGLSType );
        if ( M_stabilizationGLSParameter )
            subPt.put( "paramter method", M_stabilizationGLSParameter->method() );
        p.put_child( "Finite element stabilization", subPt );
    }

    if ( !this->isStationary() )
    {
        subPt.clear();
        subPt.put( "initial time", this->timeStepBase()->timeInitial() );
        subPt.put( "final time", this->timeStepBase()->timeFinal() );
        subPt.put( "time step", this->timeStepBase()->timeStep() );
        subPt.put( "type", M_timeStepping );
        p.put_child( "Time Discretization", subPt );
    }


    // Algebraic Solver
    if ( M_algebraicFactory )
    {
        subPt.clear();
        M_algebraicFactory->updateInformationObject( subPt );
        p.put_child( "Algebraic Solver", subPt );
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
tabulate::Table
HEAT_CLASS_TEMPLATE_TYPE::tabulateInformation( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    std::vector<std::pair<std::string,tabulate::Table>> tabInfoSections;

    if ( jsonInfo.contains("Environment") )
        tabInfoSections.push_back( std::make_pair( "Environment", super_type::super_model_base_type::tabulateInformation( jsonInfo.at("Environment"), tabInfoProp ) ) );

    if ( jsonInfo.contains("Physics") )
    {
        tabulate::Table tabInfoPhysics;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoPhysics, jsonInfo.at("Physics"), tabInfoProp );
        tabInfoSections.push_back( std::make_pair( "Physics",  tabInfoPhysics ) );
    }

    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfoSections.push_back( std::make_pair( "Materials Properties", this->materialsProperties()->tabulateInformation(jsonInfo.at("Materials Properties"), tabInfoProp ) ) );

    tabInfoSections.push_back( std::make_pair( "Boundary conditions",  tabulate::Table{} ) );

    if ( jsonInfo.contains("Meshes") )
       tabInfoSections.push_back( std::make_pair( "Meshes", super_type::super_model_meshes_type::tabulateInformation( jsonInfo.at("Meshes"), tabInfoProp ) ) );

    if ( jsonInfo.contains("Function Spaces") )
    {
        auto const& jsonInfoFunctionSpaces = jsonInfo.at("Function Spaces");
        tabulate::Table tabInfoFunctionSpaces;
        tabInfoFunctionSpaces.add_row({"Temperature"});
        tabInfoFunctionSpaces.add_row({ TabulateInformationTools::FromJSON::tabulateFunctionSpace( jsonInfoFunctionSpaces.at( "Temperature" ), tabInfoProp ) });
        tabInfoSections.push_back( std::make_pair( "Function Spaces",  tabInfoFunctionSpaces ) );
    }

    if ( jsonInfo.contains("Time Discretization") )
    {
        tabulate::Table tabInfoTimeDiscr;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoTimeDiscr, jsonInfo.at("Time Discretization"), tabInfoProp );
        tabInfoSections.push_back( std::make_pair( "Time Discretization",  tabInfoTimeDiscr ) );
    }

    if ( jsonInfo.contains("Finite element stabilization") )
    {
        tabulate::Table tabInfoStab;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoStab, jsonInfo.at("Finite element stabilization"), tabInfoProp );
        tabInfoSections.push_back( std::make_pair( "Finite element stabilization",  tabInfoStab ) );
    }

    if ( jsonInfo.contains( "Algebraic Solver" ) )
        tabInfoSections.push_back( std::make_pair( "Algebraic Solver", model_algebraic_factory_type::tabulateInformation( jsonInfo.at("Algebraic Solver"), tabInfoProp ) ) );

    return TabulateInformationTools::createSections( tabInfoSections, (boost::format("Toolbox Heat : %1%")%this->keyword()).str() );
}


HEAT_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
HEAT_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
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
    if ( M_algebraicFactory )
        *_ostr << M_algebraicFactory->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}

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

    M_bcDirichlet.setParameterValues( paramValues );
    M_bcNeumann.setParameterValues( paramValues );
    M_bcRobin.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );
    for ( auto & [matName, velConvExpr] : M_exprVelocityConvection )
        velConvExpr.setParameterValues( paramValues );

    this->log("Heat","setParameterValues", "finish");
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    M_bcDirichletMarkerManagement.clearMarkerDirichletBC();
    M_bcNeumannMarkerManagement.clearMarkerNeumannBC();
    M_bcRobinMarkerManagement.clearMarkerRobinBC();

    this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
        M_bcDirichletMarkerManagement.addMarkerDirichletBC("elimination", name(d), markers(d) );
    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "Neumann" );
    for( auto const& d : this->M_bcNeumann )
        M_bcNeumannMarkerManagement.addMarkerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d),markers(d));

    this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "temperature", "Robin" );
    for( auto const& d : this->M_bcRobin )
        M_bcRobinMarkerManagement.addMarkerRobinBC( name(d),markers(d) );

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "VolumicForces" );

    auto mesh = this->mesh();
    auto XhTemperature = this->spaceTemperature();
    std::set<std::string> temperatureMarkers;

    // strong Dirichlet bc on temperature from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) );
        temperatureMarkers.insert( listMark.begin(), listMark.end() );
    }
    auto meshMarkersTemperatureByEntities = detail::distributeMarkerListOnSubEntity( mesh, temperatureMarkers );

    // on topological faces
    auto const& listMarkedFacesTemperature = std::get<0>( meshMarkersTemperatureByEntities );
    if ( !listMarkedFacesTemperature.empty() )
        this->updateDofEliminationIds( "temperature", XhTemperature, markedfaces( mesh,listMarkedFacesTemperature ) );

}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Heat","solve", "start");
    this->timerTool("Solve").start();

    this->setStartBlockSpaceIndex( 0 );

    M_blockVectorSolution.updateVectorFromSubVectors();

    if ( this->materialsProperties()->hasThermalConductivityDependingOnSymbol( "heat_T" ) )
        M_algebraicFactory->solve( "Newton", M_blockVectorSolution.vectorMonolithic() );
    else
        M_algebraicFactory->solve( "LinearSystem", M_blockVectorSolution.vectorMonolithic() );

    M_blockVectorSolution.localize();

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
    if ( !M_algebraicFactory )
        return;
    if ( M_timeStepping == "Theta" )
    {
        M_timeStepThetaSchemePreviousContrib->zero();
        M_blockVectorSolution.updateVectorFromSubVectors();
        ModelAlgebraic::DataUpdateResidual dataResidual( M_blockVectorSolution.vectorMonolithic(), M_timeStepThetaSchemePreviousContrib, true, false );
        dataResidual.addInfo( prefixvm( this->prefix(), "time-stepping.evaluate-residual-without-time-derivative" ) );
        this->setStartBlockSpaceIndex( 0 );
        M_algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", false );
        M_algebraicFactory->evaluateResidual( dataResidual );
        M_algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", true );

        if ( M_stabilizationGLS )
        {
            auto & dataInfos = M_algebraicFactory->dataInfos();
            *dataInfos.vectorInfo( "time-stepping.previous-solution" ) = *M_blockVectorSolution.vectorMonolithic();
            dataInfos.addParameterValuesInfo( "time-stepping.previous-parameter-values", M_currentParameterValues );
        }
    }
}


} // end namespace FeelModels
} // end namespace Feel