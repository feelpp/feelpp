/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/heat/heat.hpp>
//#include <feel/feelfilters/loadgmshmesh.hpp>
//#include <feel/feelfilters/geotool.hpp>

//#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
//#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>

#include <feel/feelmodels/modelcore/modelmeasuresnormevaluation.hpp>

namespace Feel
{
namespace FeelModels
{

HEAT_CLASS_TEMPLATE_DECLARATIONS
HEAT_CLASS_TEMPLATE_TYPE::Heat( std::string const& prefix,
                                bool buildMesh,
                                WorldComm const& worldComm,
                                std::string const& subPrefix,
                                ModelBaseRepository const& modelRep )
    :
    super_type( prefix, worldComm, subPrefix, modelRep ),
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

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"Heat.info") );
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // build mesh
    if ( buildMesh )
        this->initMesh();
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

    M_fieldTemperature.reset( new element_temperature_type(M_Xh,"U"));

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
size_type
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

    // load an initial solution from a math expr
    auto initialSolution = this->modelProperties().initialConditions().getScalarFields( "temperature", "" );
    for( auto const& d : initialSolution )
    {
        auto rangeElt = (marker(d).empty())? M_rangeMeshElements : markedelements(this->mesh(),marker(d));
        this->fieldTemperaturePtr()->on(_range=rangeElt,_expr=expression(d),_geomap=this->geomap());
    }
    if ( Environment::vm().count( prefixvm(this->prefix(),"initial-solution.temperature").c_str() ) )
    {
        auto myexpr = expr( soption(_prefix=this->prefix(),_name="initial-solution.temperature"),
                            "",this->worldComm(),this->repository().expr() );
        this->fieldTemperaturePtr()->on(_range=M_rangeMeshElements,_expr=myexpr);
    }

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

    // post-process
    this->initPostProcess();

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh->worldComm() );

    // vector solution
    M_blockVectorSolution.resize( 1 );
    M_blockVectorSolution(0) = this->fieldTemperaturePtr();
    // init petsc vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
        this->initAlgebraicFactory();

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

    std::string suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_bdfTemperature = bdf( _space=this->spaceTemperature(),
                            _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"temperature"+suffixName)),
                            _prefix=this->prefix(),
                            // don't use the fluid.bdf {initial,final,step}time but the general bdf info, the order will be from fluid.bdf
                            _initial_time=this->timeInitial(),
                            _final_time=this->timeFinal(),
                            _time_step=this->timeStep(),
                            _restart=this->doRestart(),
                            _restart_path=this->restartPath(),
                            _restart_at_last_save=this->restartAtLastSave(),
                            _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );

    M_bdfTemperature->setPathSave( (fs::path(this->rootRepository()) /
                                    fs::path( prefixvm(this->prefix(), (boost::format("bdf_o_%1%_dt_%2%")%this->timeStep() %M_bdfTemperature->bdfOrder()).str() ) ) ).string() );

    // start or restart time step scheme
    if (!this->doRestart())
    {
        // start time step
        M_bdfTemperature->start(this->fieldTemperature());
        // up current time
        this->updateTime( M_bdfTemperature->time() );
    }
    else
    {
        // start time step
        M_bdfTemperature->restart();
        // load a previous solution as current solution
        *this->fieldTemperaturePtr() = M_bdfTemperature->unknown(0);
        // up initial time
        this->setTimeInitial( M_bdfTemperature->timeInitial() );
        // up current time
        this->updateTime( M_bdfTemperature->time() );
    }

    double tElapsed = this->timerTool("Constructor").stop("initTimeStep");
    this->log("Heat","initTimeStep", (boost::format("finish in %1% s") %tElapsed).str() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
std::set<std::string>
HEAT_CLASS_TEMPLATE_TYPE::postProcessFieldExported( std::set<std::string> const& ifields, std::string const& prefix ) const
{
    std::set<std::string> res;
    for ( auto const& o : ifields )
    {
        if ( o == prefixvm(prefix,"temperature") || o == prefixvm(prefix,"all") )
            res.insert( "temperature" );
        if ( o == prefixvm(prefix,"velocity-convection") || o == prefixvm(prefix,"all") )
            res.insert( "velocity-convection" );
        if ( o == prefixvm(prefix,"thermal-conductivity") || o == prefixvm(prefix,"all") )
            res.insert( "thermal-conductivity" );
        if ( o == prefixvm(prefix,"density") || o == prefixvm(prefix,"all") )
            res.insert( "density" );
        if ( o == prefixvm(prefix,"pid") || o == prefixvm(prefix,"all") )
            res.insert( "pid" );
    }
    return res;
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("Heat","initPostProcess", "start");
    this->timerTool("Constructor").start();

    std::string modelName = "heat";
    M_postProcessFieldExported = this->postProcessFieldExported( this->modelProperties().postProcess().exports( modelName ).fields() );

    if ( !M_postProcessFieldExported.empty() )
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

    pt::ptree ptree = this->modelProperties().postProcess().pTree( modelName );
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
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( modelName ) )
    {
        auto const& ptPos = evalPoints.pointPosition();
        node_type ptCoord(3);
        for ( int c=0;c<3;++c )
            ptCoord[c]=ptPos.value()(c);

        auto const& fields = evalPoints.fields();
        for ( std::string const& field : fields )
        {
            if ( field == "temperature" )
            {
                if ( !M_postProcessMeasuresContextTemperature )
                    M_postProcessMeasuresContextTemperature.reset( new context_temperature_type( spaceTemperature()->context() ) );
                int ctxId = M_postProcessMeasuresContextTemperature->nPoints();
                M_postProcessMeasuresContextTemperature->add( ptCoord );
                std::string ptNameExport = (boost::format("temperature_%1%")%ptPos.name()).str();
                this->postProcessMeasuresEvaluatorContext().add("temperature", ctxId, ptNameExport );
            }
        }
    }


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
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
boost::shared_ptr<std::ostringstream>
HEAT_CLASS_TEMPLATE_TYPE::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
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
        this->addMarkerDirichletBC("elimination", marker(d) );
    this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "Neumann" );
    for( auto const& d : this->M_bcNeumann )
        this->addMarkerNeumannBC(NeumannBCShape::SCALAR,marker(d));

    this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "temperature", "Robin" );
    for( auto const& d : this->M_bcRobin )
        this->addMarkerRobinBC( marker(d) );

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "VolumicForces" );

    auto mesh = this->mesh();
    auto XhTemperature = this->spaceTemperature();
    auto & dofsWithValueImposedTemperature = M_dofsWithValueImposed["temperature"];
    dofsWithValueImposedTemperature.clear();
    std::set<std::string> temperatureMarkers;

    // strong Dirichlet bc on temperature from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = this->markerDirichletBCByNameId( "elimination",marker(d) );
        temperatureMarkers.insert( listMark.begin(), listMark.end() );
    }
    auto meshMarkersTemperatureByEntities = detail::distributeMarkerListOnSubEntity( mesh, temperatureMarkers );

    // on topological faces
    auto const& listMarkedFacesTemperature = std::get<0>( meshMarkersTemperatureByEntities );
    for ( auto const& faceWrap : markedfaces(mesh,listMarkedFacesTemperature ) )
    {
        auto const& face = unwrap_ref( faceWrap );
        auto facedof = XhTemperature->dof()->faceLocalDof( face.id() );
        for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
            dofsWithValueImposedTemperature.insert( it->index() );
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
    this->log("Heat","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->exportFields( time );

    this->exportMeasures( time );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("Heat","exportResults", "finish");
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::exportFields( double time )
{
    bool hasFieldToExport = this->updateExportedFields( M_exporter, M_postProcessFieldExported, time );
    if ( hasFieldToExport )
        M_exporter->save();
}
HEAT_CLASS_TEMPLATE_DECLARATIONS
bool
HEAT_CLASS_TEMPLATE_TYPE::updateExportedFields( export_ptrtype exporter, std::set<std::string> const& fields, double time )
{
    if ( !exporter ) return false;
    if ( !exporter->doExport() ) return false;

    bool hasFieldToExport = false;
    if ( fields.find( "temperature" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"temperature"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"temperature")),
                                     this->fieldTemperature() );
        hasFieldToExport = true;
    }
    if ( fields.find( "pid" ) != fields.end() )
    {
        exporter->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
        hasFieldToExport = true;
    }

    if ( fields.find( "velocity-convection" ) != fields.end() )
    {
        if ( ( M_doExportVelocityConvection || M_doExportAll ) && this->fieldVelocityConvectionIsOperational() )
        {
            exporter->step( time )->add( prefixvm(this->prefix(),"velocity-convection"),
                                         prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity-convection")),
                                         this->fieldVelocityConvection() );
            hasFieldToExport = true;
        }
    }
    if ( fields.find( "thermal-conductivity" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"thermal-conductivity"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"thermal-conductivity")),
                                     this->thermalProperties()->fieldThermalConductivity() );
        hasFieldToExport = true;
    }
    if ( fields.find( "density" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"density"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"density")),
                                     this->thermalProperties()->fieldRho() );
        hasFieldToExport = true;
    }
    return hasFieldToExport;
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::exportMeasures( double time )
{
    bool hasMeasure = false;
    std::string modelName = "heat";

    // compute measures
    for ( auto const& ppForces : M_postProcessMeasuresForces )
    {
        CHECK( ppForces.meshMarkers().size() == 1 ) << "TODO";
        auto const& u = this->fieldTemperature();
        auto kappa = idv(this->thermalProperties()->fieldThermalConductivity());
        double heatFlux = integrate(_range=markedfaces(this->mesh(),ppForces.meshMarkers() ),
                                    _expr=kappa*gradv(u)*N() ).evaluate()(0,0);
        std::string name = ppForces.name();
        this->postProcessMeasuresIO().setMeasure("NormalHeatFlux_"+name,heatFlux);
        hasMeasure = true;
    }

    // point measures
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( modelName ) )
    {
        auto const& ptPos = evalPoints.pointPosition();
        if ( !ptPos.hasExpression() )
            continue;
        node_type ptCoord(3);
        for ( int c=0;c<3;++c )
            ptCoord[c]=ptPos.value()(c);

        auto const& fields = evalPoints.fields();
        for ( std::string const& field : fields )
        {
            if ( field == "temperature" )
            {
                std::string ptNameExport = (boost::format("temperature_%1%")%ptPos.name()).str();
                int ptIdInCtx = this->postProcessMeasuresEvaluatorContext().ctxId("temperature",ptNameExport);
                if ( ptIdInCtx >= 0 )
                    M_postProcessMeasuresContextTemperature->replace( ptIdInCtx, ptCoord );
            }
        }
    }
    if ( M_postProcessMeasuresContextTemperature && this->postProcessMeasuresEvaluatorContext().has("temperature") )
    {
        auto evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextTemperature,
                                                _expr=idv(this->fieldTemperature()) );
        for ( int ctxId=0;ctxId<M_postProcessMeasuresContextTemperature->nPoints();++ctxId )
        {
            if ( !this->postProcessMeasuresEvaluatorContext().has( "temperature", ctxId ) ) continue;
            std::string ptNameExport = this->postProcessMeasuresEvaluatorContext().name( "temperature",ctxId );
            this->postProcessMeasuresIO().setMeasure( ptNameExport, evalAtNodes( ctxId ) );
            hasMeasure = true;
        }
    }

    for ( auto const& ppNorm : this->modelProperties().postProcess().measuresNorm( modelName ) )
    {
        std::string const& field = ppNorm.field();
        auto range = ppNorm.markers().empty()? M_rangeMeshElements : markedelements(this->mesh(),ppNorm.markers() );
        std::map<std::string,double> resPpNorms;
        if ( field == "temperature" )
            measureNormEvaluation( range, this->fieldTemperature(), ppNorm, resPpNorms );
        for ( auto const& resPpNorm : resPpNorms )
        {
            this->postProcessMeasuresIO().setMeasure( resPpNorm.first, resPpNorm.second );
            hasMeasure = true;
        }
    }

    if ( hasMeasure )
    {
        if ( !this->isStationary() )
            this->postProcessMeasuresIO().setMeasure( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
    }
}


HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateBdf()
{
    this->log("Heat","updateBdf", "start");
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    int previousTimeOrder = this->timeStepBdfTemperature()->timeOrder();

    M_bdfTemperature->next( this->fieldTemperature() );

    int currentTimeOrder = this->timeStepBdfTemperature()->timeOrder();

    this->updateTime( this->timeStepBdfTemperature()->time() );

    // update velocity convection id symbolic expr exist and  depend only of time
    this->updateFieldVelocityConvection( true );

    // maybe rebuild cst jacobian or linear
    if ( M_algebraicFactory &&
         previousTimeOrder!=currentTimeOrder &&
         this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        if (!this->rebuildCstPartInLinearSystem())
        {
            if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".Heat","updateBdf", "do rebuildCstLinearPDE",
                                                       this->worldComm(),this->verboseAllProc());
            M_algebraicFactory->rebuildCstLinearPDE(this->blockVectorSolution().vectorMonolithic());
        }
    }

    this->timerTool("TimeStepping").stop("updateBdf");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("Heat","updateBdf", "finish");
}



} // end namespace FeelModels
} // end namespace Feel
