/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/thermodyn/thermodynbase.hpp>
//#include <feel/feelfilters/loadgmshmesh.hpp>
//#include <feel/feelfilters/geotool.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

namespace Feel
{
namespace FeelModels
{

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::ThermoDynamicsBase( std::string const& prefix,
                                                            bool buildMesh,
                                                            WorldComm const& worldComm,
                                                            std::string const& subPrefix,
                                                            std::string const& rootRepository )
    :
    super_type( prefix, worldComm, subPrefix, rootRepository ),
    M_thermalProperties( new thermalproperties_type( prefix ) )
{
    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoDynamicsConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoDynamicsSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoDynamicsPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".ThermoDynamicsTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::build()
{
    this->log("ThermoDynamics","build", "start" );

    // create or reload mesh
    this->createMesh();
    // functionSpaces and elements
    this->createFunctionSpaces();
    // bdf time schema
    this->createTimeDiscretisation();
    //export
    this->createExporters();

    this->log("ThermoDynamics","build", "finish" );
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::loadMesh( mesh_ptrtype mesh )
{
    this->log("ThermoDynamics","loadMesh", "start" );

    // create or reload mesh
    if (this->doRestart()) this->createMesh();
    else this->M_mesh = mesh;
    // functionSpaces and elements
    this->createFunctionSpaces();
    // bdf time schema
    this->createTimeDiscretisation();
    //export
    this->createExporters();

    this->log("ThermoDynamics","loadMesh", "finish" );
}



THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_fieldVelocityConvectionIsUsed = boption(_name="use_velocity-convection",_prefix=this->prefix()) ||
        Environment::vm().count(prefixvm(this->prefix(),"velocity-convection").c_str());
    M_fieldVelocityConvectionIsIncompressible = boption(_name="velocity-convection_is_incompressible",_prefix=this->prefix());

    M_doExportAll = boption(_name="do_export_all",_prefix=this->prefix());
    M_doExportVelocityConvection = boption(_name="do_export_velocity-convection",_prefix=this->prefix());
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::createMesh()
{
    this->log("ThermoDynamics","createMesh", "start");
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("createMesh");
    this->log("ThermoDynamics","createMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::createFunctionSpaces()
{
    this->log("ThermoDynamics","createFunctionSpaces", "start" );
    this->timerTool("Constructor").start();

    // functionspace
    M_Xh = space_temperature_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );

    M_fieldTemperature.reset( new element_temperature_type(M_Xh,"U"));

    if ( this->fieldVelocityConvectionIsUsed() )
        this->updateForUseFunctionSpacesVelocityConvection();

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh->worldComm() );

    double tElpased = this->timerTool("Constructor").stop("createSpaces");
    this->log("ThermoDynamics","createFunctionSpaces",(boost::format("finish in %1% s")%tElpased).str() );
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::updateForUseFunctionSpacesVelocityConvection()
{
    if ( !M_XhVelocityConvection )
        M_XhVelocityConvection = space_velocityconvection_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    if ( !M_fieldVelocityConvection )
    {
        M_fieldVelocityConvection.reset( new element_velocityconvection_type(M_XhVelocityConvection,"VelocityConvection"));
        // load the field velocity convection from a math expr
        if ( Environment::vm().count(prefixvm(this->prefix(),"velocity-convection").c_str()) )
        {
            std::string pathGinacExpr = this->directoryLibSymbExpr() + "/velocity-convection";
            M_exprVelocityConvection /*auto myexpr*/ = expr<nDim,1>( soption(_prefix=this->prefix(),_name="velocity-convection"),
                                                                     this->modelProperties().parameters().toParameterValues(), pathGinacExpr );
            this->updateFieldVelocityConvection();
            //M_fieldVelocityConvection->on(_range=elements(this->mesh()),_expr=*M_exprVelocityConvection/*myexpr*/);
        }
    }
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::updateFieldVelocityConvection( bool onlyExprWithTimeSymbol )
{
    if ( M_exprVelocityConvection.get_ptr() == 0 )
        return;

    if ( onlyExprWithTimeSymbol && !M_exprVelocityConvection->expression().hasSymbol("t") )
        return;

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    M_exprVelocityConvection->setParameterValues( paramValues );
    M_fieldVelocityConvection->on(_range=elements(this->mesh()),_expr=*M_exprVelocityConvection);
}


THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::createTimeDiscretisation()
{
    this->log("ThermoDynamics","createTimeDiscretisation", "start" );
    this->timerTool("Constructor").start();

    std::string suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_bdfTemperature = bdf( _vm=Environment::vm(), _space=this->spaceTemperature(),
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

    double tElpased = this->timerTool("Constructor").stop("createSpaces");
    this->log("ThermoDynamics","createTimeDiscretisation",(boost::format("finish in %1% s")%tElpased).str() );
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::createExporters()
{
    this->log("ThermoDynamics","createExporters", "start");
    this->timerTool("Constructor").start();

    std::string geoExportType="static";//change_coords_only, change, static
    M_exporter = exporter( _mesh=this->mesh(),
                           _name="Export",
                           _geo=geoExportType,
                           _path=this->exporterPath() );

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("ThermoDynamics","createExporters",(boost::format("finish in %1% s")%tElpased).str() );
}



THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    myblockGraph(0,0) = stencil(_test=this->spaceTemperature(),
                                _trial=this->spaceTemperature() )->graph();
    return myblockGraph;
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
size_type
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    size_type res = this->spaceTemperature()->nLocalDofWithGhost();
    return res;
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory, model_algebraic_factory_type::appli_ptrtype const& app )
{
    this->log("ThermoDynamics","init", "start" );
    this->timerTool("Constructor").start();

    M_XhScalarP0 = space_scalar_P0_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    M_thermalProperties->initFromSpace( M_XhScalarP0 );
    M_thermalProperties->updateFromModelMaterials( this->modelProperties().materials() );

    if ( this->fieldVelocityConvectionIsUsed() )
        this->updateForUseFunctionSpacesVelocityConvection();

    // load an initial solution from a math expr
    if ( Environment::vm().count(prefixvm(this->prefix(),"initial-solution.temperature").c_str()) )
    {
        std::string pathGinacExpr = this->directoryLibSymbExpr() + "/initial-solution.temperature";
        auto myexpr = expr( soption(_prefix=this->prefix(),_name="initial-solution.temperature"),
                            this->modelProperties().parameters().toParameterValues(), pathGinacExpr );
        this->fieldTemperature()->on(_range=elements(this->mesh()),_expr=myexpr);
    }

    // vector solution
    M_blockVectorSolution.resize( 1 );
    M_blockVectorSolution(0) = this->fieldTemperature();
    // init petsc vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );
    // start or restart time step scheme and exporter
    if (!this->doRestart())
    {
        // start time step
        M_bdfTemperature->start(*this->fieldTemperature());
        // up current time
        this->updateTime( M_bdfTemperature->time() );
    }
    else
    {
        // start time step
        M_bdfTemperature->restart();
        // load a previous solution as current solution
        *this->fieldTemperature() = M_bdfTemperature->unknown(0);
        // up initial time
        this->setTimeInitial( M_bdfTemperature->timeInitial() );
        // up current time
        this->updateTime( M_bdfTemperature->time() );
    }

    // post-process
    this->initPostProcess();

    // algebraic solver
    if ( buildModelAlgebraicFactory )
    {
        // matrix graph of non zero
        auto graph = this->buildBlockMatrixGraph()(0,0);
        // tool for assembly and solver
        M_algebraicFactory.reset( new model_algebraic_factory_type(app,this->backend(),
                                                                   graph, graph->mapRow().indexSplit() ) );
    }

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("ThermoDynamics","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    // restart exporter
    if (this->doRestart() )
        this->restartExporters();

    bool hasMeasure = false;

    //  heat flux measures
    auto const& ptree = this->modelProperties().postProcess().pTree();
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
                    this->postProcessMeasuresIO().setMeasure("NormalHeatFlux_"+name,0.);
                    hasMeasure = true;
                }
            }
        }
    }

    // point measures
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint() )
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
                this->postProcessMeasuresIO().setMeasure( ptNameExport, 0. );
                hasMeasure = true;
            }
        }
    }


    if ( hasMeasure )
    {
        if ( !this->isStationary() )
            this->postProcessMeasuresIO().setParameter( "time", this->timeInitial() );
        // start or restart measure file
        if (!this->doRestart())
            this->postProcessMeasuresIO().start();
        else if ( !this->isStationary() )
            this->postProcessMeasuresIO().restart( "time", this->timeInitial() );
    }

}


THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::restartExporters()
{
    if (this->doRestart() && this->restartPath().empty() )
    {
        if ( M_exporter->doExport() )
            M_exporter->restart(this->timeInitial());
    }
}


THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
boost::shared_ptr<std::ostringstream>
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||----------Info : ThermoDynamics---------------||"
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
           << "\n     -- msh filename      : " << this->mshfileStr()
           << "\n     -- number of element : " << M_mesh->numGlobalElements()
           << "\n     -- order             : " << nOrderGeo;
    *_ostr << "\n   Space Temperature Discretization"
           << "\n     -- order         : " << nOrderPoly
           << "\n     -- number of dof : " << M_Xh->nDof() << " (" << M_Xh->nLocalDof() << ")";
    if ( this->fieldVelocityConvectionIsUsedAndOperational() )
        *_ostr << "\n   Space Velocity Convection Discretization"
               << "\n     -- number of dof : " << M_XhVelocityConvection->nDof() << " (" << M_XhVelocityConvection->nLocalDof() << ")";
    if ( M_algebraicFactory )
        *_ostr << M_algebraicFactory->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}


THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("ThermoDynamics","solve", "start");
    this->timerTool("Solve").start();

    M_algebraicFactory->linearSolver(this->blockVectorSolution().vector());
    M_blockVectorSolution.localize();

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("ThermoDynamics","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}



THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    if ( !M_exporter->doExport() ) return;

    this->log("ThermoDynamics","exportResults", "start");
    this->timerTool("PostProcessing").start();

    M_exporter->step( time )->add( prefixvm(this->prefix(),"temperature"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"temperature")),
                                   *this->fieldTemperature() );
    if ( ( M_doExportVelocityConvection || M_doExportAll ) && this->fieldVelocityConvectionIsOperational() )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"velocity-convection"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity-convection")),
                                       *this->fieldVelocityConvection() );
    }
    M_exporter->save();

    this->exportMeasures( time );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("ThermoDynamics","exportResults", "finish");
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::exportMeasures( double time )
{
    bool hasMeasure = false;

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
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint() )
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


    if ( hasMeasure )
    {
        if ( !this->isStationary() )
            this->postProcessMeasuresIO().setParameter( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
    }
}


THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::updateBdf()
{
    this->log("ThermoDynamics","updateBdf", "start");
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    int previousTimeOrder = this->timeStepBdfTemperature()->timeOrder();

    M_bdfTemperature->next( *this->fieldTemperature() );

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
            if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","updateBdf", "do rebuildCstLinearPDE",
                                                       this->worldComm(),this->verboseAllProc());
            M_algebraicFactory->rebuildCstLinearPDE(this->blockVectorSolution().vector());
        }
    }

    this->timerTool("TimeStepping").stop("updateBdf");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("ThermoDynamics","updateBdf", "finish");
}









THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("ThermoDynamics","updateLinearPDE", "start"+sc);
    boost::mpi::timer thetimer;

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = *this->fieldTemperature();
    auto const& v = *this->fieldTemperature();

    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=this->rowStartInVector() );

    //--------------------------------------------------------------------------------------------------//

    if ( buildCstPart )
    {
        //double kappa = this->thermalProperties()->cstThermalConductivity();
        auto kappa = idv(this->thermalProperties()->fieldThermalConductivity());
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= kappa*inner(gradt(u),grad(v)),
                       _geomap=this->geomap() );
    }

    if ( this->fieldVelocityConvectionIsUsedAndOperational() && !buildCstPart )
    {
        //double thecoeff = this->thermalProperties()->cstRho()*this->thermalProperties()->cstHeatCapacity();
        auto thecoeff = idv(this->thermalProperties()->fieldRho())*idv(this->thermalProperties()->fieldHeatCapacity());
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= thecoeff*(gradt(u)*idv(this->fieldVelocityConvection()))*id(v),
                       _geomap=this->geomap() );


        //viscous dissipation
        if ( true )
        {
            double mu = 1.;
            auto defv = sym(gradv( this->fieldVelocityConvection() ) );
            auto defv2 = inner(defv,defv);

            if ( !this->fieldVelocityConvectionIsIncompressible() )
            {
#if 0
                bilinearForm_PatternCoupled +=
                    integrate( _range=elements(mesh),
                               _expr= thecoeff*(idt(u)*divv(this->fieldVelocityConvection()))*id(v),
                               _geomap=this->geomap() );
#endif
                myLinearForm +=
                    integrate( _range=elements(mesh),
                               _expr= 2*mu*defv2*id(v),
                               _geomap=this->geomap() );
            }
            else
            {
                auto incomp2 = pow( divv( this->fieldVelocityConvection() ),2 );
                myLinearForm +=
                    integrate( _range=elements(mesh),
                               _expr= 2*mu*(defv2-(1./3.)*incomp2)*id(v),
                               _geomap=this->geomap() );
            }
        }
    }

    if (!this->isStationary())
    {
        bool buildNonCstPart=!buildCstPart;
        bool BuildNonCstPart_Form2TransientTerm = buildNonCstPart;
        bool BuildNonCstPart_Form1TransientTerm = buildNonCstPart;
        if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        {
            BuildNonCstPart_Form2TransientTerm = buildCstPart;
        }

        if (BuildNonCstPart_Form2TransientTerm)
        {
            //double thecoeff = this->thermalProperties()->cstRho()*this->thermalProperties()->cstHeatCapacity()*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
            auto thecoeff = idv(this->thermalProperties()->fieldRho())*idv(this->thermalProperties()->fieldHeatCapacity())*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= thecoeff*idt(u)*id(v),
                           _geomap=this->geomap() );
        }

        if (BuildNonCstPart_Form1TransientTerm)
        {
            //double thecoeff = this->thermalProperties()->cstRho()*this->thermalProperties()->cstHeatCapacity();
            auto thecoeff = idv(this->thermalProperties()->fieldRho())*idv(this->thermalProperties()->fieldHeatCapacity());
            auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
            myLinearForm +=
                integrate( _range=elements(mesh),
                           _expr= thecoeff*idv(rhsTimeStep)*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//

    // update source term
    this->updateSourceTermLinearPDE(F, buildCstPart);

    // update bc
    this->updateWeakBCLinearPDE(A,F,buildCstPart);
    if ( !buildCstPart && _doBCStrongDirichlet)
        this->updateBCStrongDirichletLinearPDE(A,F);


    double timeElapsed = thetimer.elapsed();
    this->log("ThermoDynamics","updateLinearPDE",
              "finish in "+(boost::format("%1% s") % timeElapsed).str() );
}


} // end namespace FeelModels
} // end namespace Feel
