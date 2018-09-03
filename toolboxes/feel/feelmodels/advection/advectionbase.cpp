/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
*/

#include <feel/feelmodels/advection/advectionbase.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
#include <feel/feelmodels/advection/advectionbasestabilisation.cpp>

namespace Feel {
namespace FeelModels {

//----------------------------------------------------------------------------//
// Static member initialization
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
const std::map<std::string, AdvectionStabMethod> 
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::AdvectionStabMethodIdMap = {
    {"NONE", AdvectionStabMethod::NONE},
    {"GALS", AdvectionStabMethod::GALS},
    {"gls", AdvectionStabMethod::GALS},
    {"CIP", AdvectionStabMethod::CIP},
    {"SUPG", AdvectionStabMethod::SUPG},
    {"supg", AdvectionStabMethod::SUPG},
    {"SGS", AdvectionStabMethod::SGS}
};

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Construction
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::AdvectionBase( 
        std::string const& prefix,
        WorldComm const& worldComm,
        std::string const& subPrefix,
        std::string const& rootRepository )
:
    super_type( prefix, worldComm, subPrefix, rootRepository ),
    M_isBuilt(false),
    M_isUpdatedForUse(false),
    M_diffusionReactionModel( new diffusionreaction_model_type( prefix ) ),
    M_gamma1(std::pow(nOrder, -3.5))
{
    this->loadParametersFromOptionsVm();

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".AdvectionConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".AdvectionSolve.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
}
    
//----------------------------------------------------------------------------//
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::build()
{
    this->log("Advection","build", "start");
    
    // Mesh
    this->createMesh();
    // Function spaces
    this->createFunctionSpaces();
    // Algebraic data
    this->createAlgebraicData();
    // Bdf time scheme
    this->createTimeDiscretization();
    // Physical parameters
    this->createOthers();
    // Exporters
    this->createExporters();

    M_isBuilt = true;

    this->log("Advection","build", "finish");
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::build( mesh_ptrtype const& mesh)
{
    this->log("Advection","build(mesh)", "start");
    
    // Mesh
    M_mesh = mesh;
    // Function spaces
    this->createFunctionSpaces();
    // Algebraic data
    this->createAlgebraicData();
    // Bdf time scheme
    this->createTimeDiscretization();
    // Physical parameters
    this->createOthers();
    // Exporters
    this->createExporters();

    M_isBuilt = true;

    this->log("Advection","build(mesh)", "finish");
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::build( space_advection_ptrtype const& space)
{
    this->log("Advection","build(space)", "start");
    
    // Mesh
    M_mesh = space->mesh();
    // Function spaces
    if( this->stabilizationMethod() == AdvectionStabMethod::CIP )
    {
        CHECK(space->extendedDofTable()) << "CIP stabilization can only be used with extended dof table built.";
    }
    M_Xh = space;
    this->createFunctionSpaces();
    // Algebraic data
    this->createAlgebraicData();
    // Bdf time scheme
    this->createTimeDiscretization();
    // Physical parameters
    this->createOthers();
    // Exporters
    this->createExporters();

    M_isBuilt = true;

    this->log("Advection","build(space)", "finish");
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::build( space_advection_ptrtype const& space, space_advection_velocity_ptrtype const& spaceAdvectionVelocity )
{
    this->log("Advection","build(space)", "start");
    
    // Mesh
    M_mesh = space->mesh();
    // Function spaces
    if( this->stabilizationMethod() == AdvectionStabMethod::CIP )
    {
        CHECK(space->extendedDofTable()) << "CIP stabilization can only be used with extended dof table built.";
    }
    M_Xh = space;
    M_XhAdvectionVelocity = spaceAdvectionVelocity;
    this->createFunctionSpaces();
    // Algebraic data
    this->createAlgebraicData();
    // Bdf time scheme
    this->createTimeDiscretization();
    // Physical parameters
    this->createOthers();
    // Exporters
    this->createExporters();

    M_isBuilt = true;

    this->log("Advection","build(space)", "finish");
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::init(bool buildModelAlgebraicFactory, model_algebraic_factory_type::model_ptrtype const& app )
{
    //if ( M_isUpdatedForUse ) return;

    this->log("Advection","init", "start" );
    this->timerTool("Constructor").start();

    if( !M_isBuilt )
        this->build();


    if ( this->hasAdvection() )
    {
        if( (this->stabilizationMethod() == AdvectionStabMethod::GALS) ||
                (this->stabilizationMethod() == AdvectionStabMethod::SUPG) )
        {
            typedef StabilizationGLSParameter<mesh_type, nOrder> stab_gls_parameter_impl_type;
            M_stabilizationGLSParameter.reset( new stab_gls_parameter_impl_type( this->mesh(),prefixvm(this->prefix(),"stabilization-gls.parameter") ) );
            M_stabilizationGLSParameter->init();
        }
        if( this->stabilizationMethod() == AdvectionStabMethod::CIP )
        {
            M_stabilizationCIPCoefficient = doption( _name="stabilization-cip.coeff", _prefix=this->prefix() );
        }
    }

    // Vector solution
    this->buildVectorSolution();


    // Time step
    this->initTimeStep();
    // Post-process
    this->initPostProcess();

    // Algebraic factory
    if( buildModelAlgebraicFactory )
    {
        // matrix sparsity graph
        auto graph = this->buildMatrixGraph();
        M_algebraicFactory.reset( new model_algebraic_factory_type( app, this->backend(), graph, graph->mapRow().indexSplit() ) );
    }

    M_isUpdatedForUse = true;

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("Advection","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

//ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
//void
//ADVECTIONBASE_CLASS_TEMPLATE_TYPE::initFromMesh(
        //mesh_ptrtype const& mesh,
        //bool buildModelAlgebraicFactory, 
        //model_algebraic_factory_type::appli_ptrtype const& app )
//{
    ////if ( M_isUpdatedForUse ) return;

    //this->log("Advection","initFromMesh", "start" );
    //this->timerTool("Constructor").start();

    //if( !M_isBuilt )
        //this->build( mesh );

    //// Vector solution
    //M_blockVectorSolution.resize(1);
    //M_blockVectorSolution(0) = this->fieldSolutionPtr();
    //M_blockVectorSolution.buildVector( this->backend() );

    //// Time step
    //this->initTimeStep();
    
    //// Algebraic factory
    //if( buildModelAlgebraicFactory )
    //{
        //// matrix sparsity graph
        //auto graph = this->buildMatrixGraph();
        
        //M_algebraicFactory.reset( new model_algebraic_factory_type( app, this->backend(), graph, graph->mapRow().indexSplit() ) );
    //}

    //M_isUpdatedForUse = true;

    //double tElapsedInit = this->timerTool("Constructor").stop("init");
    //if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    //this->log("Advection","initFromMesh",(boost::format("finish in %1% s")%tElapsedInit).str() );
//}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    this->log("Advection","loadParametersFromOptionsVm", "start");
    
    M_hasSourceAdded = false;

    // Model
    std::string advection_model = this->modelProperties().model();
    if ( Environment::vm().count(prefixvm(this->prefix(),"model").c_str()) )
        advection_model = soption(_name="model",_prefix=this->prefix());
    if( !advection_model.empty() )
        this->setModelName( advection_model );
    // Solver
    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        this->setSolverName( soption(_name="solver",_prefix=this->prefix()) );

    // Stabilization method
    const std::string stabmeth = soption( _name="stabilization.method", _prefix=this->prefix() );
    CHECK(AdvectionStabMethodIdMap.count(stabmeth)) << stabmeth <<" is not in the list of possible stabilization methods\n";
    M_stabMethod = AdvectionStabMethodIdMap.at(stabmeth);

    M_doExportAll = boption(_name="export-all",_prefix=this->prefix());
    M_doExportAdvectionVelocity = boption(_name="export-advection-velocity",_prefix=this->prefix());
    M_doExportDiffusionCoefficient = boption(_name="export-diffusion", _prefix=this->prefix());
    M_doExportReactionCoefficient = boption(_name="export-reaction", _prefix=this->prefix());
    M_doExportSourceField = boption(_name="export-source", _prefix=this->prefix());
    
    this->log("Advection","loadParametersFromOptionsVm", "finish");
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::createMesh()
{
    this->log("Advection","createMesh","start");
    this->timerTool("Constructor").start();
    
    createMeshModel<mesh_type>(*this, M_mesh, this->fileNameMeshPath() );
    CHECK( M_mesh ) << "mesh generation failed";

    double tElapsed = this->timerTool("Constructor").stop("create");
    this->log("Advection","createMesh", (boost::format("finish in %1% s") %tElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::createFunctionSpaces()
{
    this->log("Advection","createFunctionSpaces","start");
    this->timerTool("Constructor").start();
    
    // Advection function space
    if( !M_Xh )
    {
        std::vector<bool> extendedDT( space_advection_type::nSpaces, false );
        if( this->stabilizationMethod() == AdvectionStabMethod::CIP )
        {
            this->log("Advection","createFunctionSpaces", "use buildDofTableMPIExtended on advection" );
            extendedDT[0] = true;
        }
        M_Xh = space_advection_type::New( 
                _mesh=M_mesh, 
                _worldscomm=this->worldsComm(), 
                _extended_doftable=extendedDT,
                _periodicity=this->periodicity()
                );
    }

    // Advection velocity function space
    if( !M_XhAdvectionVelocity )
    {
        M_XhAdvectionVelocity = space_advection_velocity_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    }

    // P0d function space
    if ( !M_spaceP0d )
    {
        M_spaceP0d = space_P0d_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    }

    double tElapsed = this->timerTool("Constructor").stop("create");
    this->log("Advection","createFunctionSpaces", (boost::format("finish in %1% s") %tElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::createAlgebraicData()
{
    this->log("Advection","createAlgebraicData","start");
    this->timerTool("Constructor").start();
    
    // Backend
    M_backend = backend_type::build( soption(_name="backend"), this->prefix(), M_Xh->worldComm() );

    double tElapsed = this->timerTool("Constructor").stop("create");
    this->log("Advection","createAlgebraicData", (boost::format("finish in %1% s") %tElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::createTimeDiscretization()
{
    this->log("Advection","createTimeDiscretization", "start" );
    this->timerTool("Constructor").start();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
         suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_bdf = bdf( _vm=Environment::vm(), _space=M_Xh,
                       _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"phi"+suffixName)),
                       _prefix=this->prefix(),
                       // don't use the advection.bdf {initial,final,step}time but the general bdf info, the order will be from advection.bdf
                       _initial_time=this->timeInitial(),
                       _final_time=this->timeFinal(),
                       _time_step=this->timeStep(),
                       _restart=this->doRestart(),
                       _restart_path=this->restartPath(),
                       _restart_at_last_save=this->restartAtLastSave(),
                       _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
    M_bdf->setfileFormat( myFileFormat );
    M_bdf->setPathSave( (fs::path(this->rootRepository()) /
                               fs::path( prefixvm(this->prefix(), (boost::format("bdf_o_%1%_dt_%2%")%M_bdf->bdfOrder()%this->timeStep() ).str() ) ) ).string() );

    double tElapsed = this->timerTool("Constructor").stop("createTimeDiscr");
    this->log("Advection","createTimeDiscretization", (boost::format("finish in %1% s") %tElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::createOthers()
{
    M_rangeMeshElements = elements(this->mesh());

    M_fieldSolution.reset( new element_advection_type(M_Xh, "phi") );

    // Advection velocity 
    M_fieldAdvectionVelocity.reset( new element_advection_velocity_type(M_XhAdvectionVelocity, "AdvectionVelocity") );
    // P0d
    M_spaceP0d = space_P0d_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    
    // Load the field velocity convection from a math expr
    if ( Environment::vm().count(prefixvm(this->prefix(),"advection-velocity").c_str()) )
    {
        std::string pathGinacExpr = this->directoryLibSymbExpr() + "/advection-velocity";
        M_exprAdvectionVelocity = expr<nDim,1>( soption(_prefix=this->prefix(),_name="advection-velocity"),
                                                                 this->modelProperties().parameters().toParameterValues(), pathGinacExpr );
        //this->updateFieldVelocityConvection();
        M_fieldAdvectionVelocity->on( _range=this->rangeMeshElements(), _expr=*M_exprAdvectionVelocity );
    }

    // Source term
    M_fieldSource.reset( new element_advection_type(M_Xh, "SourceAdded") );
    // Diffusion-reaction model
    M_diffusionReactionModel->initFromMesh( this->mesh(), this->useExtendedDofTable() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
int
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::buildBlockVectorSolution()
{
    int nBlock = this->nBlockMatrixGraph();
    M_blockVectorSolution.resize(nBlock);
    M_blockVectorSolution(0) = this->fieldSolutionPtr();
    int cptBlock=1;
    return cptBlock;
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::buildVectorSolution()
{
    this->buildBlockVectorSolution();
    M_blockVectorSolution.buildVector( this->backend() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::setPeriodicity( periodicity_type const& p )
{
    if( M_isUpdatedForUse )
       LOG(WARNING) << "Setting periodicity after initialization ! You need to init() again !";
    M_periodicity = p;
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::createExporters()
{
    this->log("Advection","createExporters", "start");
    this->timerTool("Constructor").start();

    std::string geoExportType = this->geoExportType();//change_coords_only, change, static
    M_exporter = exporter( _mesh=this->mesh(),
                           _name="Export",
                           _geo=geoExportType,
                           _path=this->exporterPath() );

    double tElapsed = this->timerTool("Constructor").stop("createExporters");
    this->log("Advection","createExporters",(boost::format("finish in %1% s")%tElapsed).str() );
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Model and solver
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::setModelName( std::string const& type )
{
    if( type != M_modelName )
        this->setNeedToRebuildCstPart(true);

    if ( type == "Advection" )
    {
        M_modelName = "Advection";
        M_solverName = "LinearSystem";
    }
    else if ( type == "Advection-Diffusion" )
    {
        M_modelName = "Advection-Diffusion";
        M_solverName = "LinearSystem";
    }
    else if ( type == "Advection-Reaction" )
    {
        M_modelName = "Advection-Reaction";
        M_solverName = "LinearSystem";
    }
    else if ( type == "Advection-Diffusion-Reaction" )
    {
        M_modelName = "Advection-Diffusion-Reaction";
        M_solverName = "LinearSystem";
    }
    else if ( type == "Diffusion" )
    {
        M_modelName = "Diffusion";
        M_solverName = "LinearSystem";
    }
    else if ( type == "Diffusion-Reaction" )
    {
        M_modelName = "Diffusion-Reaction";
        M_solverName = "LinearSystem";
    }
    else if ( type == "Reaction" )
    {
        M_modelName = "Reaction";
        M_solverName = "LinearSystem";
    }
    else
        CHECK( false ) << "invalid modelName " << type << "\n";
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
bool
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::hasAdvection() const
{
    return (this->modelName() == "Advection" || this->modelName() == "Advection-Diffusion" ||
            this->modelName() == "Advection-Reaction" || this->modelName() == "Advection-Diffusion-Reaction");
}
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
bool
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::hasDiffusion() const
{
    return (this->modelName() == "Diffusion" || this->modelName() == "Advection-Diffusion" ||
            this->modelName() == "Diffusion-Reaction" || this->modelName() == "Advection-Diffusion-Reaction");
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
bool
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::hasReaction() const
{
    return (this->modelName() == "Reaction" || this->modelName() == "Advection-Reaction" ||
            this->modelName() == "Diffusion-Reaction" || this->modelName() == "Advection-Diffusion-Reaction");
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::setSolverName( std::string const& type )
{
    if( type != M_solverName )
        this->setNeedToRebuildCstPart(true);

    if ( type == "LinearSystem" )
    {
        M_solverName = "LinearSystem";
    }
    //else if ( type == "Advection-Diffusion" )
    //{
        //M_solverName = "Advection-Diffusion";
    //}
    //else if ( type == "Advection-Diffusion-Reaction" )
    //{
        //M_solverName = "Advection-Diffusion-Reaction";
    //}
    else
        CHECK( false ) << "invalid solverName " << type << "\n";
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
bool
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::useExtendedDofTable() const
{
    if ( this->worldComm().localSize() == 1 ) return false;
    bool useExtendedDofTable=false;
    for ( bool hasExt : M_Xh->extendedDofTableComposite() )
        useExtendedDofTable = useExtendedDofTable || hasExt;
    return useExtendedDofTable;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Algebraic data
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
size_type
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::matrixPattern() const
{
    size_type pat = size_type(Pattern::COUPLED);

    if( this->stabilizationMethod() == AdvectionStabMethod::CIP )
    {
        pat = size_type(Pattern::EXTENDED);
    }

    return pat;
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
int
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = 1;
    return nBlock;
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("AdvectionBase","buildBlockMatrixGraph", "start" );
    int nBlock = this->nBlockMatrixGraph();

    BlocksBaseGraphCSR blockGraph(nBlock,nBlock);
    int indexBlock=0;

    blockGraph(indexBlock, indexBlock) = stencil( 
            _test=this->functionSpace(), _trial=this->functionSpace(),
            _pattern=this->matrixPattern(),
            _diag_is_nonzero=(nBlock==1),
            _close=(nBlock==1)
            )->graph();

    this->log("Advectionbase","buildBlockMatrixGraph", "finish" );
    return blockGraph;
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
typename ADVECTIONBASE_CLASS_TEMPLATE_TYPE::graph_ptrtype
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::buildMatrixGraph() const
{
    auto blockGraph = this->buildBlockMatrixGraph();
    blockGraph.close();

    if ( blockGraph.nRow() == 1 && blockGraph.nCol() == 1 )
        return blockGraph(0,0);
    else
        return graph_ptrtype( new graph_type( blockGraph ) );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
size_type
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    auto res = this->functionSpace()->nLocalDofWithGhost();
    return res;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Time scheme
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateTimeStepBDF() 
{
    this->log("Advection","updateTimeStepBDF", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    int previousTimeOrder = this->timeStepBDF()->timeOrder();

    M_bdf->next( *M_fieldSolution );

    int currentTimeOrder = this->timeStepBDF()->timeOrder();

    this->updateTime( M_bdf->time() );

    // maybe rebuild linear
    if ( M_algebraicFactory &&
         previousTimeOrder!=currentTimeOrder &&
         this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        if (!this->rebuildCstPartInLinearSystem())
        {
            this->log("Advection","updateTimeStepBDF", "do rebuildCstLinearPDE" );
            M_algebraicFactory->rebuildCstLinearPDE(M_fieldSolution);
        }
    }

    this->timerTool("TimeStepping").stop("updateBdf");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("Advection","updateTimeStepBDF", "finish" );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    // start or restart time step scheme
    if (!this->doRestart())
    {
        // start time step
        M_bdf->start(*M_fieldSolution);
        // up current time
        this->updateTime( M_bdf->time() );
    }
    else
    {
        // start time step
        M_bdf->restart();
        // load a previous solution as current solution
        *M_fieldSolution = M_bdf->unknown(0);
        // up initial time
        this->setTimeInitial( M_bdf->timeInitial() );
        // restart exporter
        //this->restartPostProcess();
        // up current time
        this->updateTime( M_bdf->time() );

        this->log("Advection","initTimeStep", "restart bdf/exporter done" );
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Algebraic model updates
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    using namespace Feel::vf;

    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _BuildCstPart = data.buildCstPart();
    sparse_matrix_ptrtype& A_extended = data.matrixExtended();
    bool _BuildExtendedPart = data.buildExtendedPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

    bool build_AdvectiveTerm = BuildNonCstPart;
    bool build_DiffusionTerm = BuildNonCstPart;
    bool build_ReactionTerm = BuildNonCstPart;
    //bool build_SourceTerm = BuildNonCstPart;
    //bool build_BoundaryNeumannTerm = BuildNonCstPart;

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Advection","updateLinearPDE", "start"+sc );
    this->timerTool("Solve").start();

    auto mesh = this->mesh();
    auto space = this->functionSpace();

    auto const& advection_velocity = this->fieldAdvectionVelocity();

    auto const& phi = this->fieldSolution();
    auto const& psi = this->fieldSolution();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto linearForm = form1( _test=space, _vector=F );

    // Advection
    if(this->hasAdvection() && build_AdvectiveTerm)
    {
        this->timerTool("Solve").start();
        
        bilinearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=inner((gradt(phi)*idv(advection_velocity)), id(psi)),
                _geomap=this->geomap()
                );

        double timeElapsedAdvection = this->timerTool("Solve").stop();
        this->log("Advection","updateLinearPDE","assembly advection terms in "+(boost::format("%1% s") %timeElapsedAdvection).str() );
    }

    // Diffusion
    if( this->hasDiffusion() && build_DiffusionTerm )
    {
        this->timerTool("Solve").start();

        auto const& D = this->diffusionReactionModel()->fieldDiffusionCoeff();

        bilinearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=inner(trans(gradt(phi)),idv(D)*trans(grad(psi))),
                _geomap=this->geomap()
                );

        double timeElapsedDiffusion = this->timerTool("Solve").stop();
        this->log("Advection","updateLinearPDE","assembly diffusion terms in "+(boost::format("%1% s") %timeElapsedDiffusion).str() );
    }

    // Reaction
    if( this->hasReaction() && build_ReactionTerm )
    {
        this->timerTool("Solve").start();

        auto const& R = this->diffusionReactionModel()->fieldReactionCoeff();
        
        bilinearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=inner(idv(R)*idt(phi), id(psi)),
                _geomap=this->geomap()
                );

        double timeElapsedReaction = this->timerTool("Solve").stop();
        this->log("Advection","updateLinearPDE","assembly reaction terms in "+(boost::format("%1% s") %timeElapsedReaction).str() );
    }

    // Transient terms
    if (!this->isStationary())
    {
        this->timerTool("Solve").start();
       
        this->updateLinearPDETransient( A, F, BuildCstPart ); 
        
        double timeElapsedTransient = this->timerTool("Solve").stop();
        this->log("Advection","updateLinearPDE","assembly transient terms in "+(boost::format("%1% s") %timeElapsedTransient).str() );
    }

    // Source term
    this->updateSourceTermLinearPDE( data );
    if( this->hasSourceAdded() && !BuildCstPart )
    {
        linearForm +=
            integrate( _range=this->rangeMeshElements(),
                       _expr= inner(idv(M_fieldSource),id(psi)),
                       _geomap=this->geomap() );
    }

    // User-defined additional terms
    this->updateLinearPDEAdditional( A, F, BuildCstPart );

    // Stabilization
    if ( this->hasAdvection() )
        this->updateLinearPDEStabilization( data );

    // Boundary conditions
    this->updateWeakBCLinearPDE(A, F, BuildCstPart);
    if ( BuildNonCstPart && _doBCStrongDirichlet)
        this->updateBCStrongDirichletLinearPDE(A,F);

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Advection","updateLinearPDE","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateLinearPDEStabilization( DataUpdateLinear & data ) const
{
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;
    std::string sc=(BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Advection","updateLinearPDEStabilization", "start"+sc );
    this->timerTool("Solve").start();

    uint16_type adrt = ADREnum::Advection;
    if( this->hasReaction() ) adrt |= ADREnum::Reaction;
    if( this->hasDiffusion() && nOrder >= 2 ) adrt |= ADREnum::Diffusion;
    ADREnum adrtype = static_cast<ADREnum>( adrt );
    
    switch ( this->stabilizationMethod() )
    {
        case AdvectionStabMethod::NONE : { break; } // remove -Wswitch warning

        case AdvectionStabMethod::GALS :
        {
            if( adrtype == Advection ) ADRDetails::updateLinearPDEStabilizationGLS<ADRTypes::Advection>( *this, data );
            if( adrtype == AdvectionDiffusion) ADRDetails::updateLinearPDEStabilizationGLS<ADRTypes::AdvectionDiffusion>( *this, data );
            if( adrtype == AdvectionReaction) ADRDetails::updateLinearPDEStabilizationGLS<ADRTypes::AdvectionReaction>( *this, data );
            if( adrtype == AdvectionDiffusionReaction) ADRDetails::updateLinearPDEStabilizationGLS<ADRTypes::AdvectionDiffusionReaction>( *this, data );
        } //GALS
        break ;

        case AdvectionStabMethod::SUPG :
        {
            if( adrtype == Advection ) ADRDetails::updateLinearPDEStabilizationSUPG<ADRTypes::Advection>( *this, data );
            if( adrtype == AdvectionDiffusion) ADRDetails::updateLinearPDEStabilizationSUPG<ADRTypes::AdvectionDiffusion>( *this, data );
            if( adrtype == AdvectionReaction) ADRDetails::updateLinearPDEStabilizationSUPG<ADRTypes::AdvectionReaction>( *this, data );
            if( adrtype == AdvectionDiffusionReaction) ADRDetails::updateLinearPDEStabilizationSUPG<ADRTypes::AdvectionDiffusionReaction>( *this, data );
        } //SUPG
        break;

        case AdvectionStabMethod::SGS :
        {
            if( adrtype == Advection ) ADRDetails::updateLinearPDEStabilizationSGS<ADRTypes::Advection>( *this, data );
            if( adrtype == AdvectionDiffusion) ADRDetails::updateLinearPDEStabilizationSGS<ADRTypes::AdvectionDiffusion>( *this, data );
            if( adrtype == AdvectionReaction) ADRDetails::updateLinearPDEStabilizationSGS<ADRTypes::AdvectionReaction>( *this, data );
            if( adrtype == AdvectionDiffusionReaction) ADRDetails::updateLinearPDEStabilizationSGS<ADRTypes::AdvectionDiffusionReaction>( *this, data );
        } //SGS
        break;

        case AdvectionStabMethod::CIP :
        {
            if( BuildNonCstPart )
            {
                auto mesh = this->mesh();
                auto space = this->functionSpace();
                auto const& phi = this->fieldSolution();
                sparse_matrix_ptrtype & A = data.matrix();

                auto beta = idv(this->fieldAdvectionVelocity());
                auto beta_norm = vf::sqrt(trans(beta)*beta);
                double stabCoeff = this->M_stabilizationCIPCoefficient;
                auto coeff = stabCoeff * M_gamma1 * hFace() * hFace() * beta_norm;

                auto bilinearForm_PatternExtended = form2( 
                        _test=space, _trial=space, _matrix=A, _pattern=size_type(Pattern::EXTENDED) 
                        );
                bilinearForm_PatternExtended += integrate(
                        _range=internalfaces(mesh),
                        _expr=coeff * inner(jumpt(gradt(phi)), jump(grad(phi)))
                        );
            }

        } //CIP
        break;
    } //switch

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Advection","updateLinearPDEStabilization","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateLinearPDETransient( sparse_matrix_ptrtype& A, vector_ptrtype& F, bool BuildCstPart ) const
{
    auto const& mesh = this->mesh();
    auto const& space = this->functionSpace();
    auto const& phi = this->fieldSolution();
    auto const& psi = this->fieldSolution();
    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto linearForm = form1( _test=space, _vector=F );

    bool BuildNonCstPart = !BuildCstPart;
    bool build_Form2TransientTerm = BuildNonCstPart;
    bool build_Form1TransientTerm = BuildNonCstPart;
    if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        build_Form2TransientTerm = BuildCstPart;
    }

    if (build_Form2TransientTerm)
    {
        bilinearForm += integrate( 
                _range=this->rangeMeshElements(),
                _expr=M_bdf->polyDerivCoefficient(0)*inner(idt(phi),id(psi)),
                _geomap=this->geomap() 
                );
    }
    if (build_Form1TransientTerm)
    {
        linearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=inner(idv(M_bdf->polyDeriv()),id(psi)),
                _geomap=this->geomap() 
                );
    }
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
}

//----------------------------------------------------------------------------//
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateAdvectionVelocity( element_advection_velocity_ptrtype const& u )
{
    M_exprAdvectionVelocity.reset(); // remove symbolic expr
    M_fieldAdvectionVelocity = u;
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateAdvectionVelocity( element_advection_velocity_type const& u )
{
    M_exprAdvectionVelocity.reset(); // remove symbolic expr
    *M_fieldAdvectionVelocity = u;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Solve
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Advection","solve", "start" );
    this->timerTool("Solve").start();

    std::string algebraicSolverType = "LinearSystem";
    M_algebraicFactory->solve( algebraicSolverType, M_blockVectorSolution.vectorMonolithic() ); 
    
    M_blockVectorSolution.localize();

    double tElapsed = this->timerTool("Solve").stop("solve");
    this->log("Advection","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Export results
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    if (this->doRestart() && this->restartPath().empty() )
    {
        if ( M_exporter->doExport() ) M_exporter->restart(this->timeInitial());
    }
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::exportMeasures( double time )
{
    this->exportMeasuresImpl( time );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::exportMeasuresImpl( double time )
{
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->exportResultsImpl( time );
    this->exportMeasures( time );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::exportResultsImpl( double time )
{
    if ( !M_exporter->doExport() ) return;

    this->log("Advection","exportResults", "start");
    this->timerTool("PostProcessing").start();

    M_exporter->step( time )->add( prefixvm(this->prefix(),"phi"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"phi")),
                                   this->fieldSolution() );
    if ( ( M_doExportAdvectionVelocity || M_doExportAll ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"advection_velocity"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"advection_velocity")),
                                       this->fieldAdvectionVelocity() );
    }
    if ( ( M_doExportDiffusionCoefficient || M_doExportAll ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"diffusion_coeff"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"diffusion_coeff")),
                                       this->diffusionReactionModel()->fieldDiffusionCoeff() );
    }
    if ( ( M_doExportReactionCoefficient || M_doExportAll ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"reaction_coeff"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"reaction_coeff")),
                                       this->diffusionReactionModel()->fieldReactionCoeff() );
    }
    if ( ( M_doExportSourceField || M_doExportAll ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"source"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"source")),
                                       *M_fieldSource );
    }
    M_exporter->save();

    double tElapsed = this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("Advection","exportResults", (boost::format("finish in %1% s")%tElapsed).str() );
}

} // namespace FeelModels
} // namespace Feel
