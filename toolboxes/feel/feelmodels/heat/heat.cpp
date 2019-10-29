/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/heat/heat.hpp>
//#include <feel/feelfilters/loadgmshmesh.hpp>
//#include <feel/feelfilters/geotool.hpp>

//#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
//#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>

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
    M_thermalProperties( new thermalproperties_type( prefix, this->repository().expr() ) )
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
    M_fieldVelocityConvectionIsUsed = boption(_name="use_velocity-convection",_prefix=this->prefix()) ||
        Environment::vm().count(prefixvm(this->prefix(),"velocity-convection").c_str());
    M_fieldVelocityConvectionIsIncompressible = boption(_name="velocity-convection_is_incompressible",_prefix=this->prefix());

    M_doExportAll = boption(_name="do_export_all",_prefix=this->prefix());
    M_doExportVelocityConvection = boption(_name="do_export_velocity-convection",_prefix=this->prefix());

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

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("initMesh");
    this->log("Heat","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initMaterialProperties()
{
    this->log("Heat","initMaterialProperties", "start" );
    this->timerTool("Constructor").start();

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().materials().setParameterValues( paramValues );
    M_thermalProperties->updateForUse( M_mesh, this->modelProperties().materials() );

    double tElpased = this->timerTool("Constructor").stop("initMaterialProperties");
    this->log("Heat","initMaterialProperties",(boost::format("finish in %1% s")%tElpased).str() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    this->log("Heat","initFunctionSpaces", "start" );
    this->timerTool("Constructor").start();

    // functionspace
    if ( M_thermalProperties->isDefinedOnWholeMesh() )
    {
        M_rangeMeshElements = elements(M_mesh);
        M_Xh = space_temperature_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    }
    else
    {
        M_rangeMeshElements = markedelements(M_mesh, M_thermalProperties->markers());
        M_Xh = space_temperature_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
    }

    M_fieldTemperature.reset( new element_temperature_type(M_Xh,"temperature"));

    if ( this->fieldVelocityConvectionIsUsed() )
        this->updateForUseFunctionSpacesVelocityConvection();

    double tElpased = this->timerTool("Constructor").stop("initFunctionSpaces");
    this->log("Heat","initFunctionSpaces",(boost::format("finish in %1% s")%tElpased).str() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateForUseFunctionSpacesVelocityConvection()
{
    if ( !M_XhVelocityConvection )
    {
        if ( M_thermalProperties->isDefinedOnWholeMesh() )
            M_XhVelocityConvection = space_velocityconvection_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
        else
            M_XhVelocityConvection = space_velocityconvection_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm(), _range=M_rangeMeshElements );
    }

    if ( !M_fieldVelocityConvection )
    {
        M_fieldVelocityConvection.reset( new element_velocityconvection_type(M_XhVelocityConvection,"VelocityConvection"));
        // load the field velocity convection from a math expr
        if ( Environment::vm().count(prefixvm(this->prefix(),"velocity-convection").c_str()) )
        {
            M_exprVelocityConvection = expr<nDim,1>( soption(_prefix=this->prefix(),_name="velocity-convection"),
                                                     "",this->worldComm(),this->repository().expr() );
            this->updateFieldVelocityConvection();
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateFieldVelocityConvection( bool onlyExprWithTimeSymbol )
{
    if ( M_exprVelocityConvection.get_ptr() == 0 )
        return;

    if ( onlyExprWithTimeSymbol && !M_exprVelocityConvection->expression().hasSymbol("t") )
        return;

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    M_exprVelocityConvection->setParameterValues( paramValues );
    M_fieldVelocityConvection->on(_range=M_rangeMeshElements,_expr=*M_exprVelocityConvection);
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

    if ( !M_mesh )
        this->initMesh();

    this->initMaterialProperties();

    this->initFunctionSpaces();

    this->initBoundaryConditions();

    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();

    this->initInitialConditions();

    // stabilization gls
    if ( M_stabilizationGLS )
    {
        typedef StabilizationGLSParameter<mesh_type, nOrderTemperature> stab_gls_parameter_impl_type;
        M_stabilizationGLSParameter.reset( new stab_gls_parameter_impl_type( this->mesh(),prefixvm(this->prefix(),"stabilization-gls.parameter") ) );
        M_stabilizationGLSParameter->init();
    }

    // post-process
    this->initPostProcess();

    // update fields
    this->updateFields( this->symbolsExpr() );

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh->worldCommPtr() );

    // subspaces index
    size_type currentStartIndex = 0;
    this->setStartSubBlockSpaceIndex( "temperature", currentStartIndex++ );

     // vector solution
    M_blockVectorSolution.resize( 1 );
    M_blockVectorSolution(0) = this->fieldTemperaturePtr();

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
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        std::string suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    fs::path saveTsDir = fs::path(this->rootRepository())/fs::path( prefixvm(this->prefix(),prefixvm(this->subPrefix(),"ts")) );

    int bdfOrder = 1;
    if ( M_timeStepping == "BDF" )
        bdfOrder = ioption(_prefix=this->prefix(),_name="bdf.order");
    int nConsecutiveSave = std::max( 3, bdfOrder ); // at least 3 is required when restart with theta scheme 
    M_bdfTemperature = bdf( _space=this->spaceTemperature(),
                            _name="temperature"+suffixName,
                            _prefix=this->prefix(),
                            _order=bdfOrder,
                            // don't use the fluid.bdf {initial,final,step}time but the general bdf info, the order will be from fluid.bdf
                            _initial_time=this->timeInitial(),
                            _final_time=this->timeFinal(),
                            _time_step=this->timeStep(),
                            _restart=this->doRestart(),
                            _restart_path=this->restartPath(),
                            _restart_at_last_save=this->restartAtLastSave(),
                            _save=this->tsSaveInFile(), _format=myFileFormat, _freq=this->tsSaveFreq(),
                            _n_consecutive_save=nConsecutiveSave );
    M_bdfTemperature->setPathSave( ( saveTsDir/"temperature" ).string() );

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
HEAT_CLASS_TEMPLATE_TYPE::initInitialConditions()
{
    if ( !this->doRestart() )
    {
        std::vector<element_temperature_ptrtype> icTemperatureFields;
        if ( this->isStationary() )
            icTemperatureFields = { this->fieldTemperaturePtr() };
        else
            icTemperatureFields = M_bdfTemperature->unknowns();
        this->updateInitialConditions( "temperature", M_rangeMeshElements, this->symbolsExpr(), icTemperatureFields );

        if ( !this->isStationary() )
            *this->fieldTemperaturePtr() = M_bdfTemperature->unknown(0);

        if ( Environment::vm().count( prefixvm(this->prefix(),"initial-solution.temperature").c_str() ) )
        {
            auto myexpr = expr( soption(_prefix=this->prefix(),_name="initial-solution.temperature"),
                                "",this->worldComm(),this->repository().expr() );
            this->fieldTemperaturePtr()->on(_range=M_rangeMeshElements,_expr=myexpr);
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("Heat","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->setPostProcessExportsAllFieldsAvailable( {"temperature","velocity-convection","thermal-conductivity","density","pid"} );
    this->setPostProcessSaveAllFieldsAvailable( {"temperature","velocity-convection","thermal-conductivity","density"} );
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
                // get list of marker
                std::set<std::string> markerSet;
                std::string markerUnique = ptreeLevel1.second.template get_value<std::string>();
                if ( markerUnique.empty() )
                {
                    for (auto const& ptreeMarker : ptreeLevel1.second )
                    {
                        std::string marker = ptreeMarker.second.template get_value<std::string>();
                        markerSet.insert( marker );
                    }
                }
                else
                {
                    markerSet.insert( markerUnique );
                }
                // save forces measure for each marker
                for ( std::string const& marker : markerSet )
                {
                    ModelMeasuresForces myPpForces;
                    myPpForces.addMarker( marker );
                    myPpForces.setName( marker );
                    //std::cout << "add ppHeatFlux with name " << marker<<"\n";
                    std::string name = myPpForces.name();
                    M_postProcessMeasuresForces.push_back( myPpForces );
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
    // init petsc vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );

    M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );

    if ( M_timeStepping == "Theta" )
    {
        M_timeStepThetaSchemePreviousContrib = this->backend()->newVector(M_blockVectorSolution.vectorMonolithic()->mapPtr() );
        M_algebraicFactory->addVectorResidualAssembly( M_timeStepThetaSchemePreviousContrib, 1.0, "Theta-Time-Stepping-Previous-Contrib", true );
        M_algebraicFactory->addVectorLinearRhsAssembly( M_timeStepThetaSchemePreviousContrib, -1.0, "Theta-Time-Stepping-Previous-Contrib", false );
        if ( M_stabilizationGLS )
            M_algebraicFactory->dataInfos().addVectorInfo( prefixvm( this->prefix(),"time-stepping.previous-solution"), this->backend()->newVector( M_blockVectorSolution.vectorMonolithic()->mapPtr() ) );
    }

}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateInformationObject( pt::ptree & p )
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.get_child_optional( "Prefix" ) )
        return;
    p.put( "Prefix", this->prefix() );
    p.put( "Root Repository", this->rootRepository() );

    // Physical Model
    pt::ptree subPt, subPt2;
    subPt.put( "time mode", std::string( (this->isStationary())?"Stationary":"Transient") );
    subPt.put( "velocity-convection",  std::string( (this->fieldVelocityConvectionIsUsedAndOperational())?"Yes":"No" ) );
    p.put_child( "Physical Model", subPt );

    // Boundary Conditions
    subPt.clear();
    subPt2.clear();
    this->updateInformationObjectDirichletBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    subPt2.clear();
    this->updateInformationObjectNeumannBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    subPt2.clear();
    this->updateInformationObjectRobinBC( subPt2 );
    for( const auto& ptIter : subPt2 )
        subPt.put_child( ptIter.first, ptIter.second );
    p.put_child( "Boundary Conditions",subPt );

    // Materials parameters
    subPt.clear();
    this->thermalProperties()->updateInformationObject( subPt );
    p.put_child( "Materials parameters", subPt );

    // Mesh and FunctionSpace
    subPt.clear();
    subPt.put("filename", this->meshFile());
    M_mesh->putInformationObject( subPt );
    p.put( "Mesh",  M_mesh->journalSectionName() );
    p.put( "FunctionSpace Temperature",  M_Xh->journalSectionName() );
    if ( this->fieldVelocityConvectionIsUsedAndOperational() )
        p.put( "FunctionSpace Velocity Convection", M_XhVelocityConvection->journalSectionName() );
    if ( M_stabilizationGLS )
    {
        subPt.clear();
        subPt.put( "type", M_stabilizationGLSType );
        if ( M_stabilizationGLSParameter )
            subPt.put( "paramter method", M_stabilizationGLSParameter->method() );
        p.put_child( "Finite element stabilization", subPt );
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
           << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient")
           << "\n     -- velocity-convection : " << std::string( (this->fieldVelocityConvectionIsUsedAndOperational())?"Yes":"No" );
    *_ostr << "\n   Boundary conditions"
           << this->getInfoDirichletBC()
           << this->getInfoNeumannBC()
           << this->getInfoRobinBC();
    *_ostr << this->thermalProperties()->getInfoMaterialParameters()->str();
    *_ostr << "\n   Mesh Discretization"
           << "\n     -- mesh filename      : " << this->meshFile()
           << "\n     -- number of element : " << M_mesh->numGlobalElements()
           << "\n     -- order             : " << nOrderGeo;
    *_ostr << "\n   Space Temperature Discretization"
           << "\n     -- order         : " << nOrderPoly
           << "\n     -- number of dof : " << M_Xh->nDof() << " (" << M_Xh->nLocalDof() << ")";
    if ( this->fieldVelocityConvectionIsUsedAndOperational() )
        *_ostr << "\n   Space Velocity Convection Discretization"
               << "\n     -- number of dof : " << M_XhVelocityConvection->nDof() << " (" << M_XhVelocityConvection->nLocalDof() << ")";
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
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().parameters().setParameterValues( paramValues );

    this->thermalProperties()->setParameterValues( paramValues );
    M_bcDirichlet.setParameterValues( paramValues );
    M_bcNeumann.setParameterValues( paramValues );
    M_bcRobin.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerRobinBC();

    this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC("elimination", name(d), markers(d) );
    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "Neumann" );
    for( auto const& d : this->M_bcNeumann )
        this->addMarkerNeumannBC(NeumannBCShape::SCALAR,name(d),markers(d));

    this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "temperature", "Robin" );
    for( auto const& d : this->M_bcRobin )
        this->addMarkerRobinBC( name(d),markers(d) );

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "VolumicForces" );

    auto mesh = this->mesh();
    auto XhTemperature = this->spaceTemperature();
    std::set<std::string> temperatureMarkers;

    // strong Dirichlet bc on temperature from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = this->markerDirichletBCByNameId( "elimination",name(d) );
        temperatureMarkers.insert( listMark.begin(), listMark.end() );
    }
    auto meshMarkersTemperatureByEntities = detail::distributeMarkerListOnSubEntity( mesh, temperatureMarkers );

    // on topological faces
    auto const& listMarkedFacesTemperature = std::get<0>( meshMarkersTemperatureByEntities );
    if ( !listMarkedFacesTemperature.empty() )
    {
        auto therange = markedfaces(mesh,listMarkedFacesTemperature );
        auto dofsToAdd = XhTemperature->dofs( therange );
        XhTemperature->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
        this->dofEliminationIdsAll("temperature",MESH_FACES).insert( dofsToAdd.begin(), dofsToAdd.end() );
        auto dofsMultiProcessToAdd = XhTemperature->dofs( therange, ComponentType::NO_COMPONENT, true );
        this->dofEliminationIdsMultiProcess("temperature",MESH_FACES).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Heat","solve", "start");
    this->timerTool("Solve").start();

    this->updateParameterValues();

    this->setStartBlockSpaceIndex( 0 );

    M_blockVectorSolution.updateVectorFromSubVectors();

    if ( this->thermalProperties()->hasThermalConductivityDependingOnSymbol( "heat_T" ) )
        M_algebraicFactory->solve( "Newton", M_blockVectorSolution.vectorMonolithic() );
    else
        M_algebraicFactory->solve( "LinearSystem", M_blockVectorSolution.vectorMonolithic() );

    M_blockVectorSolution.localize();

    this->updateFields( this->symbolsExpr() );

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
    this->exportResults( time, this->symbolsExpr() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::executePostProcessMeasures( double time )
{
    this->executePostProcessMeasures( time, this->allFields(), this->symbolsExpr() );
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

    // update velocity convection id symbolic expr exist and  depend only of time
    this->updateFieldVelocityConvection( true );

    // maybe rebuild cst jacobian or linear
    if ( M_algebraicFactory && rebuildCstAssembly )
    {
        if (!this->rebuildCstPartInLinearSystem())
        {
            this->log("Heat","updateTimeStep", "do rebuildCstLinearPDE" );
            M_algebraicFactory->rebuildCstLinearPDE(this->blockVectorSolution().vectorMonolithic());
        }
    }

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
        M_algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", false );
        M_algebraicFactory->evaluateResidual( dataResidual );
        M_algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", true );

        if ( M_stabilizationGLS )
        {
            auto & dataInfos = M_algebraicFactory->dataInfos();
            *dataInfos.vectorInfo( prefixvm( this->prefix(),"time-stepping.previous-solution") ) = *M_blockVectorSolution.vectorMonolithic();
            if ( this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                std::string convectionOseenEntry = prefixvm( this->prefix(),"time-stepping.previous-convection-velocity-field" );
                if ( !dataInfos.hasVectorInfo( convectionOseenEntry ) )
                    dataInfos.addVectorInfo( convectionOseenEntry, this->backend()->newVector( this->fieldVelocityConvection().mapPtr() ) );
                *dataInfos.vectorInfo( convectionOseenEntry ) = this->fieldVelocityConvection();
            }
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    this->updateLinearPDE( data, this->symbolsExpr() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    if ( M_bcDirichlet.empty() ) return;

    this->log("Heat","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto mesh = this->mesh();
    auto u = this->spaceTemperature()->element( U, this->rowStartInVector() );

    for( auto const& d : M_bcDirichlet )
    {
        auto theExpr = expression(d,this->symbolsExpr());
        u.on(_range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
             _expr=theExpr );
    }

    // update info for synchronization
    this->updateDofEliminationIds( "temperature", data );

    this->log("Heat","updateNewtonInitialGuess","finish" );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto const t = this->spaceTemperature()->element(XVec, this->rowStartInVector());
    this->updateJacobian( data, this->symbolsExpr(t) );
}
HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto const t = this->spaceTemperature()->element(XVec, this->rowStartInVector());
    this->updateResidual( data, this->symbolsExpr(t) );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( !this->hasMarkerDirichletBCelimination() ) return;
    //if ( this->M_bcDirichlet.empty() ) return;

    this->log("Heat","updateResidualDofElimination","start" );

    this->updateDofEliminationIds( "temperature", data );

    this->log("Heat","updateResidualDofElimination","finish" );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    if ( !this->hasMarkerDirichletBCelimination() ) return;
    //if ( this->M_bcDirichlet.empty() ) return;

    this->log("Heat","updateJacobianDofElimination","start" );

    this->updateDofEliminationIds( "temperature", data );

    this->log("Heat","updateJacobianDofElimination","finish" );
}


HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    if ( !this->hasMarkerDirichletBCelimination() ) return;
    //if ( this->M_bcDirichlet.empty() ) return;

    this->log("Heat","updateLinearPDEDofElimination","start" );

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = this->fieldTemperature();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for( auto const& d : this->M_bcDirichlet )
    {
        auto theExpr = expression(d,this->symbolsExpr());
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
                _element=u,_rhs=F,_expr=theExpr );
    }

    this->log("Heat","updateLinearPDEDofElimination","finish" );
}

} // end namespace FeelModels
} // end namespace Feel
