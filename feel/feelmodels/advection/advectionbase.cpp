/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
*/

#include <feel/feelmodels/advection/advectionbase.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

namespace Feel {
namespace FeelModels {

//----------------------------------------------------------------------------//
// Static member initialization
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
const std::map<std::string, AdvectionStabMethod> 
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::AdvectionStabMethodIdMap = {
    {"NONE", AdvectionStabMethod::NONE},
    {"GALS", AdvectionStabMethod::GALS},
    {"CIP", AdvectionStabMethod::CIP},
    {"SUPG", AdvectionStabMethod::SUPG},
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
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::init(bool buildModelAlgebraicFactory, model_algebraic_factory_type::appli_ptrtype const& app )
{
    //if ( M_isUpdatedForUse ) return;

    this->log("Advection","init", "start" );
    this->timerTool("Constructor").start();

    if( !M_isBuilt )
        this->build();

    // Vector solution
    M_blockVectorSolution.resize(1);
    M_blockVectorSolution(0) = this->fieldSolutionPtr();
    M_blockVectorSolution.buildVector( this->backend() );

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

    // Solver
    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        this->setSolverName( soption(_name="solver",_prefix=this->prefix()) );

    // Stabilization method
    const std::string stabmeth = soption( _name="advec-stab-method", _prefix=this->prefix() );
    CHECK(AdvectionStabMethodIdMap.count(stabmeth)) << stabmeth <<" is not in the list of possible stabilization methods\n";
    M_stabMethod = AdvectionStabMethodIdMap.at(stabmeth);
    M_stabCoeff = doption( _name="advec-stab-coeff", _prefix=this->prefix() );

    M_doExportAll = boption(_name="export-all",_prefix=this->prefix());
    M_doExportAdvectionVelocity = boption(_name="export-advection-velocity",_prefix=this->prefix());
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
    
    // Create advection function space
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
    M_fieldSolution.reset( new element_advection_type(M_Xh, "phi") );

    // Advection velocity 
    M_XhAdvectionVelocity = space_advection_velocity_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    M_fieldAdvectionVelocity.reset( new element_advection_velocity_type(M_XhAdvectionVelocity, "AdvectionVelocity") );
    
    // Load the field velocity convection from a math expr
    if ( Environment::vm().count(prefixvm(this->prefix(),"advection-velocity").c_str()) )
    {
        std::string pathGinacExpr = this->directoryLibSymbExpr() + "/advection-velocity";
        M_exprAdvectionVelocity = expr<nDim,1>( soption(_prefix=this->prefix(),_name="advection-velocity"),
                                                                 this->modelProperties().parameters().toParameterValues(), pathGinacExpr );
        //this->updateFieldVelocityConvection();
        M_fieldAdvectionVelocity->on( _range=elements(this->mesh()), _expr=*M_exprAdvectionVelocity );
    }

    // Source term
    M_fieldSource.reset( new element_advection_type(M_Xh, "SourceAdded") );
    // Diffusion-reaction model
    M_diffusionReactionModel->initFromMesh( this->mesh(), this->useExtendedDofTable() );
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

    std::string geoExportType="static";//change_coords_only, change, static
    M_exporter = exporter( _mesh=this->mesh(),
                           _name="Export",
                           _geo=geoExportType,
                           _path=this->exporterPath() );

    double tElpased = this->timerTool("Constructor").stop("createExporters");
    this->log("Advection","createExporters",(boost::format("finish in %1% s")%tElpased).str() );
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
    else if ( type == "Advection-Diffusion-Reaction" )
    {
        M_modelName = "Advection-Diffusion-Reaction";
        M_solverName = "LinearSystem";
    }
    else
        CHECK( false ) << "invalid modelName " << type << "\n";
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
bool
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::hasDiffusion() const
{
    return (this->modelName() == "Advection-Diffusion" || this->modelName() == "Advection-Diffusion-Reaction");
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
bool
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::hasReaction() const
{
    return (this->modelName() == "Advection-Diffusion-Reaction");
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
typename ADVECTIONBASE_CLASS_TEMPLATE_TYPE::graph_ptrtype
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::buildMatrixGraph() const
{
    return stencil( 
            _test=this->functionSpace(), _trial=this->functionSpace(),
            _pattern=this->matrixPattern(),
            //_diag_is_nonzero=true,
            _close=true
            )->graph();
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
    bool build_Form2TransientTerm = BuildNonCstPart;
    bool build_Form1TransientTerm = BuildNonCstPart;
    //bool build_SourceTerm = BuildNonCstPart;
    //bool build_BoundaryNeumannTerm = BuildNonCstPart;
    if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        build_Form2TransientTerm=BuildCstPart;
    }

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
    if(build_AdvectiveTerm)
    {
        this->timerTool("Solve").start();
        
        bilinearForm += integrate(
                _range=elements(mesh),
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
                _range=elements(mesh),
                _expr=idv(D)*inner(gradt(phi),grad(psi)),
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
                _range=elements(mesh),
                _expr=-idv(R)*inner(idt(phi), id(psi)),
                _geomap=this->geomap()
                );

        double timeElapsedReaction = this->timerTool("Solve").stop();
        this->log("Advection","updateLinearPDE","assembly reaction terms in "+(boost::format("%1% s") %timeElapsedReaction).str() );
    }

    // Transient terms
    if (!this->isStationary())
    {
        this->timerTool("Solve").start();
        
        if (build_Form2TransientTerm)
        {
            bilinearForm += integrate( 
                    _range=elements(mesh),
                    _expr=M_bdf->polyDerivCoefficient(0)*inner(idt(phi),id(psi)),
                    _geomap=this->geomap() 
                    );
        }
        if (build_Form1TransientTerm)
        {
            linearForm += integrate(
                    _range=elements(mesh),
                    _expr=inner(idv(M_bdf->polyDeriv()),id(psi)),
                    _geomap=this->geomap() 
                    );
        }
        
        double timeElapsedTransient = this->timerTool("Solve").stop();
        this->log("Advection","updateLinearPDE","assembly transient terms in "+(boost::format("%1% s") %timeElapsedTransient).str() );
    }

    // Source term
    if( this->hasSourceTerm() )
    {
        if( !this->hasSourceAdded() )
            M_fieldSource->zero();
        this->updateSourceTermLinearPDE(const_cast<element_advection_ptrtype&>(M_fieldSource), BuildCstPart);
    }

    if( BuildNonCstPart && (this->hasSourceAdded() || this->hasSourceTerm()) )
    {
        this->timerTool("Solve").start();

        linearForm += integrate( 
                _range=elements(mesh),
                _expr= inner(idv(*M_fieldSource),id(psi)),
                _geomap=this->geomap() 
                );

        double timeElapsedSource = this->timerTool("Solve").stop();
        this->log("Advection","updateLinearPDE","assembly source terms in "+(boost::format("%1% s") %timeElapsedSource).str() );
    }

    // Stabilization
    this->updateLinearPDEStabilization(A, F, BuildCstPart);

    // Boundary conditions
    this->updateWeakBCLinearPDE(A, F, BuildCstPart);
    if ( BuildNonCstPart && _doBCStrongDirichlet)
        this->updateBCStrongDirichletLinearPDE(A,F);

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Advection","updateLinearPDE","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateLinearPDEStabilization(sparse_matrix_ptrtype& A, vector_ptrtype& F, bool BuildCstPart) const
{
    using namespace Feel::vf;

    bool Build_StabTerm = !BuildCstPart;

    std::string sc=(BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Advection","updateLinearPDEStabilization", "start"+sc );
    this->timerTool("Solve").start();
    
    auto mesh = this->mesh();
    auto space = this->functionSpace();

    auto const& phi = this->fieldSolution();
    auto const& psi = this->fieldSolution();

    double stabCoeff = this->stabilizationCoefficient();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto bilinearForm_PatternExtended = form2( _test=space, _trial=space, _matrix=A, _pattern=size_type(Pattern::EXTENDED) );
    auto linearForm = form1( _test=space, _vector=F, _rowstart=this->rowStartInVector() );

    if(Build_StabTerm)
    {
        double sigma = this->isStationary() ? 0: M_bdf->polyDerivCoefficient(0);
        auto beta = idv(this->fieldAdvectionVelocity());
        auto const& D = this->diffusionReactionModel()->fieldDiffusionCoeff();
        auto const& R = this->diffusionReactionModel()->fieldReactionCoeff();
        auto beta_norm = vf::sqrt(trans(beta)*beta);
        auto polyDeriv = M_bdf->polyDeriv();

        auto f = (!this->isStationary()) * idv(polyDeriv) // transient rhs
            + (this->hasSourceAdded() || this->hasSourceTerm()) * idv(*M_fieldSource) // source
            ;
        
        switch ( this->stabilizationMethod() )
        {
            case AdvectionStabMethod::NONE : { break; } // remove -Wswitch warning

            case AdvectionStabMethod::GALS :
            {
                auto coeff  = val( stabCoeff /
                        ( 2*beta_norm*nOrder/h() + std::abs(sigma) ));

                //auto L_op = grad(psi)*beta + sigma*id(psi);
                //auto L_opt = gradt(phi)*beta + sigma*idt(phi);

                auto L_op = grad(psi) * beta // advection term
                    + (!this->isStationary()) * M_bdf->polyDerivCoefficient(0)*id(psi) // transient term
                    + (this->hasReaction()) * (-idv(R))*id(psi) // reaction term
                    + (this->hasDiffusion() && nOrder >= 2) * (-idv(D))*laplacian(psi) // diffusion term
                    ;
                auto L_opt = gradt(phi) * beta // advection term
                    + (!this->isStationary()) * M_bdf->polyDerivCoefficient(0)*idt(phi) // transient term
                    + (this->hasReaction()) * (-idv(R))*idt(phi) // reaction term
                    + (this->hasDiffusion() && nOrder >= 2 ) * (-idv(D))*laplaciant(phi) // diffusion term
                    ;

                bilinearForm += integrate(
                        _range=elements(mesh),
                        _expr=coeff * inner(L_op, L_opt),
                        _geomap=this->geomap() );

                linearForm += integrate(
                        _range=elements(mesh),
                        _expr=coeff * inner(L_op, f),
                        _geomap=this->geomap() );

                break ;
            } //GALS

            case AdvectionStabMethod::SUPG :
            {
                auto coeff = val(vf::h() / (2 * beta_norm + 0.001));

                auto L_op = grad(psi) * val(beta);
                //auto L_opt = gradt(phi) * val(beta) + val(sigma) * idt(phi);

                auto L_opt = gradt(phi) * beta // advection term
                    + (!this->isStationary()) * M_bdf->polyDerivCoefficient(0)*idt(phi) // transient term
                    + (this->hasReaction()) * (-idv(R))*idt(phi) // reaction term
                    + (this->hasDiffusion() && nOrder >= 2 ) * (-idv(D))*laplaciant(phi) // diffusion term
                    ;

                bilinearForm += integrate(
                        _range=elements(M_mesh),
                        _expr=coeff * inner(L_op, L_opt) 
                        );

                linearForm += integrate(
                        _range=elements(M_mesh),
                        _expr=coeff * inner(L_op, f) 
                        );

                break;
            } //SUPG

            case AdvectionStabMethod::SGS :
            {
                auto coeff = val(1. / ( (2 * beta_norm) / vf::h()  + vf::abs(sigma) ) );

                //auto L_op = grad(psi) * val(beta) - val(sigma) * id(psi);
                //auto L_opt = gradt(phi) * val(beta) + val(sigma) * idt(phi);
                
                auto L_op = grad(psi) * beta // advection term
                    - (!this->isStationary()) * M_bdf->polyDerivCoefficient(0)*id(psi) // transient term
                    - (this->hasReaction()) * (-idv(R))*id(psi) // reaction term
                    - (this->hasDiffusion() && nOrder >= 2) * (-idv(D))*laplacian(psi) // diffusion term
                    ;
                auto L_opt = gradt(phi) * beta // advection term
                    + (!this->isStationary()) * M_bdf->polyDerivCoefficient(0)*idt(phi) // transient term
                    + (this->hasReaction()) * (-idv(R))*idt(phi) // reaction term
                    + (this->hasDiffusion() && nOrder >= 2 ) * (-idv(D))*laplaciant(phi) // diffusion term
                    ;

                bilinearForm += integrate(
                        _range=elements(M_mesh),
                        _expr=coeff * inner(L_op, L_opt)
                        );

                linearForm += integrate(
                        _range=elements(M_mesh),
                        _expr=coeff * inner(L_op, f)
                        );

                break;
            } //SGS

            case AdvectionStabMethod::CIP :
            {
                auto coeff = stabCoeff * M_gamma1 * hFace() * hFace() * beta_norm;

                bilinearForm_PatternExtended += integrate(
                        _range=internalfaces(M_mesh),
                        _expr=coeff * inner(jumpt(gradt(phi)), jump(grad(phi)))
                        );

                break;
            } //CIP

        } //switch
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Advection","updateLinearPDEStabilization","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
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
    M_algebraicFactory->solve( algebraicSolverType, M_blockVectorSolution.vector() ); 
    
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
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    if ( !M_exporter->doExport() ) return;

    this->log("Advection","exportResults", "start");
    this->timerTool("PostProcessing").start();

    M_exporter->step( time )->add( prefixvm(this->prefix(),"phi"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"phi")),
                                   this->fieldSolution() );
    if ( ( M_doExportAdvectionVelocity || M_doExportAll ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"advection-velocity"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"advection_velocity")),
                                       this->fieldAdvectionVelocity() );
    }
    if ( ( M_doExportSourceField || M_doExportAll ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"source"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"source")),
                                       *M_fieldSource );
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
    this->log("Advection","exportResults", "finish");
}

} // namespace FeelModels
} // namespace Feel
