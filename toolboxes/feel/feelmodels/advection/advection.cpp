/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
*/

#include <feel/feelmodels/advection/advection.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

namespace Feel {
namespace FeelModels {

//----------------------------------------------------------------------------//
// Static member initialization
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
const std::map<std::string, AdvectionStabMethod> 
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::AdvectionStabMethodIdMap = {
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
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::AdvDiffReac( 
        std::string const& prefix,
        std::string const& keyword,
        worldcomm_ptr_t const& worldComm,
        std::string const& subPrefix,
        ModelBaseRepository const& modelRep )
:
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelPhysics<nDim>( "adr" ),
    M_isUpdatedForUse(false),
    M_diffusionReactionModel( new diffusionreaction_model_type( prefix ) ),
    M_doProjectFieldAdvectionVelocity( false ),
    M_gamma1(std::pow(nOrder, -3.5))
{
    this->loadParametersFromOptionsVm();

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".AdvectionConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".AdvectionSolve.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
typename ADVDIFFREAC_CLASS_TEMPLATE_TYPE::self_ptrtype 
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::New( 
        std::string const& prefix,
        std::string const& keyword,
        worldcomm_ptr_t const& worldComm,
        std::string const& subPrefix,
        ModelBaseRepository const& modelRep )
{
    return std::make_shared<self_type>( prefix, keyword, worldComm, subPrefix, modelRep );
}
    
//----------------------------------------------------------------------------//
// Initialization
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    if ( M_isUpdatedForUse ) return;

    this->log("AdvDiffReac","init", "start" );
    this->timerTool("Constructor").start();

    if ( !this->modelPropertiesPtr() )
        this->setModelProperties( std::make_shared<ModelProperties>( "", this->repository().expr(), this->worldCommPtr() ) );

    if( this->modelName().empty() )
    {
        std::string advection_model = this->modelProperties().models().model().equations();
        this->setModelName( advection_model );
    }

    // Boundary conditions
    this->loadConfigBCFile();

    // Mesh
    if( !M_mesh )
        this->createMesh();

    this->initMaterialProperties();

    // Function spaces
    this->initFunctionSpaces();
    // Algebraic data
    this->initAlgebraicData();
    // Bdf time scheme
    this->initTimeDiscretization();
    // Physical parameters
    this->initOthers();

    // Stabilization
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

    // Initial conditions
    this->initInitialConditions();

    // Time step
    this->initTimeStep();
    // Post-process
    this->initPostProcess();

    // Algebraic factory
    if( buildModelAlgebraicFactory )
    {
        // matrix sparsity graph
        auto graph = this->buildMatrixGraph();
        M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(), this->backend(), graph, graph->mapRow().indexSplit() ) );
    }

    M_isUpdatedForUse = true;

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("AdvDiffReac","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

//----------------------------------------------------------------------------//
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    this->log("AdvDiffReac","loadParametersFromOptionsVm", "start");

    M_hasSourceAdded = false;

    // Model
    if ( Environment::vm().count(prefixvm(this->prefix(),"model").c_str()) )
        this->setModelName( soption(_name="model",_prefix=this->prefix()) );
    // Solver
    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        this->setSolverName( soption(_name="solver",_prefix=this->prefix()) );

    // Stabilization method
    const std::string stabmeth = soption( _name="stabilization.method", _prefix=this->prefix() );
    CHECK(AdvectionStabMethodIdMap.count(stabmeth)) << stabmeth <<" is not in the list of possible stabilization methods\n";
    M_stabMethod = AdvectionStabMethodIdMap.at(stabmeth);

    this->log("AdvDiffReac","loadParametersFromOptionsVm", "finish");
}

namespace detail {

template<int Dim, bool isVectorial, 
        typename = typename std::enable_if<!isVectorial>::type
        >
map_scalar_field<2> getBCFields( 
        BoundaryConditions const& bc, 
        std::string const& field, std::string const& type
        )
{
    return bc.getScalarFields( std::string(field), std::string(type) );
}

template<int Dim, bool isVectorial, 
        typename = typename std::enable_if<isVectorial>::type
        >
map_vector_field<Dim, 1, 2> getBCFields( 
        BoundaryConditions const& bc, 
        std::string const& field, std::string const& type
        )
{
    return bc.getVectorFields<Dim>( std::string(field), std::string(type) );
}

} // namespace detail

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->M_bcInflowMarkers.clear();

    //this->M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "advection", "Dirichlet" );
    this->M_bcDirichlet = detail::getBCFields<nDim, is_vectorial>( 
            this->modelProperties().boundaryConditions(), this->prefix(), "Dirichlet");
    for( auto const& d : this->M_bcDirichlet )
        this->addMarkerDirichletBC("elimination", name(d), markers(d) );

    //this->M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "advection", "Neumann" );
    this->M_bcNeumann = detail::getBCFields<nDim, is_vectorial>(
            this->modelProperties().boundaryConditions(), this->prefix(), "Neumann"
            );
    for( auto const& d : this->M_bcNeumann )
        this->addMarkerNeumannBC(NeumannBCShape::SCALAR,name(d),markers(d));

    //this->M_bcRobin = this->modelProperties().boundaryConditions().getScalarFieldsList( "advection", "Robin" );
    //for( auto const& d : this->M_bcRobin )
        //this->addMarkerRobinBC( marker(d) );

    for( std::string const& bcMarker: this->modelProperties().boundaryConditions().markers( this->prefix(), "inflow" ) )
        if( std::find(this->M_bcInflowMarkers.begin(), this->M_bcInflowMarkers.end(), bcMarker) == this->M_bcInflowMarkers.end() )
            this->M_bcInflowMarkers.push_back( bcMarker );

    //M_sources = this->modelProperties().boundaryConditions().template getScalarFields( "advection", "Sources" );
    this->M_sources = detail::getBCFields<nDim, is_vectorial>(
            this->modelProperties().boundaryConditions(), this->prefix(), "Sources"
            );
}

//----------------------------------------------------------------------------//
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::createMesh()
{
    this->log("AdvDiffReac","createMesh","start");
    this->timerTool("Constructor").start();
    
    createMeshModel<mesh_type>(*this, M_mesh, this->fileNameMeshPath() );
    CHECK( M_mesh ) << "mesh generation failed";
    M_isUpdatedForUse = false;
    M_rangeMeshElements = elements( M_mesh );

    double tElapsed = this->timerTool("Constructor").stop("create");
    this->log("AdvDiffReac","createMesh", (boost::format("finish in %1% s") %tElapsed).str() );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::initMaterialProperties()
{
    this->log("AdvDiffReac","initMaterialProperties", "start" );
    this->timerTool("Constructor").start();

    if ( !M_materialsProperties )
    {
        auto paramValues = this->modelProperties().parameters().toParameterValues();
        this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->prefix(), this->repository().expr() ) );
        M_materialsProperties->updateForUse( M_mesh, this->modelProperties().materials(), *this );
    }

    double tElpased = this->timerTool("Constructor").stop("initMaterialProperties");
    this->log("AdvDiffReac","initMaterialProperties",(boost::format("finish in %1% s")%tElpased).str() );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    this->log("AdvDiffReac","initFunctionSpaces","start");
    this->timerTool("Constructor").start();
    
    // AdvDiffReac function space
    if( !M_Xh )
    {
        std::vector<bool> extendedDT( space_advection_type::nSpaces, false );
        if( this->stabilizationMethod() == AdvectionStabMethod::CIP )
        {
            this->log("AdvDiffReac","initFunctionSpaces", "use buildDofTableMPIExtended on advection" );
            extendedDT[0] = true;
        }
        M_Xh = space_advection_type::New( 
                _mesh=M_mesh, 
                _worldscomm=this->worldsComm(), 
                _extended_doftable=extendedDT,
                _periodicity=this->periodicity()
                );
    }
    if( this->stabilizationMethod() == AdvectionStabMethod::CIP )
    {
        CHECK(M_Xh->extendedDofTable()) << "CIP stabilization can only be used with extended dof table built.";
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
    this->log("AdvDiffReac","initFunctionSpaces", (boost::format("finish in %1% s") %tElapsed).str() );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::initAlgebraicData()
{
    this->log("AdvDiffReac","initAlgebraicData","start");
    this->timerTool("Constructor").start();
    
    // Backend
    M_backend = backend_type::build( soption(_name="backend"), this->prefix(), M_Xh->worldCommPtr() );

    double tElapsed = this->timerTool("Constructor").stop("create");
    this->log("AdvDiffReac","initAlgebraicData", (boost::format("finish in %1% s") %tElapsed).str() );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::initTimeDiscretization()
{
    this->log("AdvDiffReac","initTimeDiscretization", "start" );
    this->timerTool("Constructor").start();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
         suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_bdf = bdf( _vm=Environment::vm(), _space=M_Xh,
                       _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"phi"+suffixName)),
                       _prefix=this->prefix(),
                       // don't use the advection.bdf {initial,final,step}time but the general bdf info
                       // the order will be from [prefix].bdf.order
                       _order=this->timeOrder(),
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
    this->log("AdvDiffReac","initTimeDiscretization", (boost::format("finish in %1% s") %tElapsed).str() );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::createPostProcessExporters()
{
    if ( !this->postProcessExportsFields().empty() )
    {
        std::string geoExportType = this->geoExportType();//change_coords_only, change, static
        M_exporter = Feel::exporter( 
                _mesh=this->mesh(),
                _name="Export",
                _geo=geoExportType,
                _worldcomm=this->functionSpace()->worldComm(),
                _path=this->exporterPath() 
                );

        if ( M_exporter->doExport() && this->doRestart() && this->restartPath().empty() )
            M_exporter->restart(this->timeInitial());
    }
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::createPostProcessMeasures()
{
    // Start measures export
    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasuresIO().restart( "time", this->timeInitial() );
        else
            this->postProcessMeasuresIO().setMeasure( "time", this->timeInitial() ); //just for have time in first column
    }
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::initOthers()
{
    M_fieldSolution.reset( new element_advection_type(M_Xh, "phi") );

    // Advection velocity 
    M_fieldAdvectionVelocity.reset( new element_advection_velocity_type(M_XhAdvectionVelocity, "AdvectionVelocity") );
    // P0d
    M_spaceP0d = space_P0d_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    
    // Load the field velocity convection from a math expr
    if ( Environment::vm().count(prefixvm(this->prefix(),"advection-velocity").c_str()) )
    {
        std::string pathGinacExpr = this->repository().expr() + "/advection-velocity";
        vector_field_expression<nDim,1,2> exprAdvectionVelocity = expr<nDim,1>( soption(_prefix=this->prefix(),_name="advection-velocity"),
                                                                 this->modelProperties().parameters().toParameterValues(), pathGinacExpr );
        //this->updateFieldVelocityConvection();
        this->updateAdvectionVelocity( exprAdvectionVelocity );
        //M_fieldAdvectionVelocity->on( _range=this->rangeMeshElements(), _expr=*M_exprAdvectionVelocity );
    }

    // Source term
    M_fieldSource.reset( new element_advection_type(M_Xh, "SourceAdded") );
    // Diffusion-reaction model
    M_diffusionReactionModel->initFromMesh( this->mesh(), this->useExtendedDofTable() );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::initInitialConditions()
{
    if( !this->doRestart() )
    {
        if( !M_initialValue )
        {
            M_initialValue.reset( new element_advection_type(this->functionSpace(), "initialValue") );
            std::vector<element_advection_ptrtype> icADRFields;
            if ( this->isStationary() )
                icADRFields = { M_initialValue };
            else
                icADRFields = M_bdf->unknowns();

            auto paramValues = this->modelProperties().parameters().toParameterValues();
            this->modelProperties().initialConditions().setParameterValues( paramValues );

            this->updateInitialConditions( this->prefix(), this->rangeMeshElements(), this->symbolsExpr(), icADRFields );

            if( !this->isStationary() )
                *M_initialValue = M_bdf->unknown(0);
        }

        *this->fieldSolutionPtr() = *M_initialValue;
    }
}

//----------------------------------------------------------------------------//
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
int
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::buildBlockVectorSolution()
{
    int nBlock = this->nBlockMatrixGraph();
    M_blockVectorSolution.resize(nBlock);
    M_blockVectorSolution(0) = this->fieldSolutionPtr();
    int cptBlock=1;
    return cptBlock;
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::buildVectorSolution()
{
    this->buildBlockVectorSolution();
    M_blockVectorSolution.buildVector( this->backend() );
}

//----------------------------------------------------------------------------//
// Mesh
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::loadMesh( mesh_ptrtype m )
{
    this->log("AdvDiffReac","loadMesh", "start");
    // create or reload mesh
    if ( this->doRestart() && !m )
    {
        this->createMesh();
    }
    else
    {
        M_mesh = m;
        M_isUpdatedForUse = false;
    }
    M_rangeMeshElements = elements(M_mesh);
    this->log("AdvDiffReac","loadMesh", "finish");
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::setMesh( mesh_ptrtype const& m )
{
    this->log("AdvDiffReac","setMesh", "start");
    M_mesh = m;
    M_rangeMeshElements = elements(M_mesh);
    M_isUpdatedForUse = false;
    this->log("AdvDiffReac","setMesh", "finish");
}

//----------------------------------------------------------------------------//
// Periodicity
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::setPeriodicity( periodicity_type const& p )
{
    if( M_isUpdatedForUse )
       LOG(WARNING) << "Setting periodicity after initialization ! You need to init() again !";
    M_periodicity = p;
}
//----------------------------------------------------------------------------//
// Model and solver
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::setModelName( std::string const& type )
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

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
bool
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::hasAdvection() const
{
    return (this->modelName() == "Advection" || this->modelName() == "Advection-Diffusion" ||
            this->modelName() == "Advection-Reaction" || this->modelName() == "Advection-Diffusion-Reaction");
}
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
bool
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::hasDiffusion() const
{
    return (this->modelName() == "Diffusion" || this->modelName() == "Advection-Diffusion" ||
            this->modelName() == "Diffusion-Reaction" || this->modelName() == "Advection-Diffusion-Reaction");
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
bool
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::hasReaction() const
{
    return (this->modelName() == "Reaction" || this->modelName() == "Advection-Reaction" ||
            this->modelName() == "Diffusion-Reaction" || this->modelName() == "Advection-Diffusion-Reaction");
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::setSolverName( std::string const& type )
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
// Spaces
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::setFunctionSpace( space_advection_ptrtype space )
{
    this->log("AdvDiffReac","setFunctionSpace", "start");
    M_Xh = space;
    this->setMesh(space->mesh());
    M_isUpdatedForUse = false;
    this->log("AdvDiffReac","setFunctionSpace", "finish");
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::setFunctionSpaceAdvectionVelocity( space_advection_velocity_ptrtype space )
{
    this->log("AdvDiffReac","setFunctionSpace", "start");
    M_XhAdvectionVelocity = space;
    M_isUpdatedForUse = false;
    this->log("AdvDiffReac","setFunctionSpace", "finish");
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
bool
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::useExtendedDofTable() const
{
    if ( this->worldComm().localSize() == 1 ) return false;
    bool useExtendedDofTable=false;
    for ( bool hasExt : M_Xh->extendedDofTableComposite() )
        useExtendedDofTable = useExtendedDofTable || hasExt;
    return useExtendedDofTable;
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
typename ADVDIFFREAC_CLASS_TEMPLATE_TYPE::element_advection_velocity_ptrtype const&
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::fieldAdvectionVelocityPtr() const
{
    if( M_doProjectFieldAdvectionVelocity )
    {
       const_cast<self_type*>(this)->M_functionProjectFieldAdvectionVelocity();
    }

    return M_fieldAdvectionVelocity;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Algebraic data
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
typename ADVDIFFREAC_CLASS_TEMPLATE_TYPE::size_type
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::matrixPattern() const
{
    size_type pat = size_type(Pattern::COUPLED);

    if( this->stabilizationMethod() == AdvectionStabMethod::CIP )
    {
        pat = size_type(Pattern::EXTENDED);
    }

    return pat;
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
int
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = 1;
    return nBlock;
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("AdvDiffReac","buildBlockMatrixGraph", "start" );
    int nBlock = this->nBlockMatrixGraph();

    BlocksBaseGraphCSR blockGraph(nBlock,nBlock);
    int indexBlock=0;

    blockGraph(indexBlock, indexBlock) = stencil( 
            _test=this->functionSpace(), _trial=this->functionSpace(),
            _pattern=this->matrixPattern(),
            _diag_is_nonzero=(nBlock==1),
            _close=(nBlock==1)
            )->graph();

    this->log("AdvDiffReac","buildBlockMatrixGraph", "finish" );
    return blockGraph;
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
typename ADVDIFFREAC_CLASS_TEMPLATE_TYPE::graph_ptrtype
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::buildMatrixGraph() const
{
    auto blockGraph = this->buildBlockMatrixGraph();
    blockGraph.close();

    if ( blockGraph.nRow() == 1 && blockGraph.nCol() == 1 )
        return blockGraph(0,0);
    else
        return graph_ptrtype( new graph_type( blockGraph ) );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
typename ADVDIFFREAC_CLASS_TEMPLATE_TYPE::size_type
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    auto res = this->functionSpace()->nLocalDofWithGhost();
    return res;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Time scheme
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::updateTimeStepBDF() 
{
    this->log("AdvDiffReac","updateTimeStepBDF", "start" );
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
            this->log("AdvDiffReac","updateTimeStepBDF", "do rebuildCstLinearPDE" );
            M_algebraicFactory->rebuildCstLinearPDE(M_fieldSolution);
        }
    }

    this->timerTool("TimeStepping").stop("updateBdf");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("AdvDiffReac","updateTimeStepBDF", "finish" );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::initTimeStep()
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

        this->log("AdvDiffReac","initTimeStep", "restart bdf/exporter done" );
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Algebraic model updates

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::updateLinearPDEStabilization( DataUpdateLinear & data ) const
{
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;
    std::string sc=(BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("AdvDiffReac","updateLinearPDEStabilization", "start"+sc );
    this->timerTool("Solve").start();

    this->M_functionAssemblyLinearStabilization( data );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("AdvDiffReac","updateLinearPDEStabilization","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const
{
// TODO
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
bool
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::hasSourceTerm() const
{
    return !this->M_sources.empty(); 
}


ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
}

//----------------------------------------------------------------------------//
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::updateAdvectionVelocity( element_advection_velocity_ptrtype const& u )
{
    this->updateAdvectionVelocity( *u );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::updateAdvectionVelocity( element_advection_velocity_type const& u )
{
    *M_fieldAdvectionVelocity = u;
    M_functionAssemblyLinearAdvection = [this]( DataUpdateLinear & data ) { 
        this->updateLinearPDEAdvection( data, idv(this->M_fieldAdvectionVelocity) );
    };
    M_functionAssemblyLinearStabilization = [this]( DataUpdateLinear & data ) {
        this->updateLinearPDEStabilization( data, idv(this->M_fieldAdvectionVelocity) );
    };
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Solve
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("AdvDiffReac","solve", "start" );
    this->timerTool("Solve").start();

    CHECK( M_isUpdatedForUse ) << "The advection toolbox must be properly initialized before using it.\n";

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();

    this->M_bcDirichlet.setParameterValues( paramValues );
    this->M_bcNeumann.setParameterValues( paramValues );
    this->M_sources.setParameterValues( paramValues );

    std::string algebraicSolverType = "LinearSystem";
    M_algebraicFactory->solve( algebraicSolverType, M_blockVectorSolution.vectorMonolithic() ); 
    
    M_blockVectorSolution.localize();

    double tElapsed = this->timerTool("Solve").stop("solve");
    this->log("AdvDiffReac","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Export results
ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::initPostProcessExportsAndMeasures()
{
    // Update post-process expressions
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    std::set<std::string> postProcessAllFieldsAvailable = { "phi", "advection-velocity", "diffusion-coeff", "reaction-coeff", "source" };
    this->setPostProcessExportsAllFieldsAvailable( postProcessAllFieldsAvailable );
    this->setPostProcessExportsPidName( "pid" );
    this->setPostProcessSaveAllFieldsAvailable( postProcessAllFieldsAvailable );
    super_type::initPostProcess();

    // Point measures
    auto fieldNamesWithSpaceADR = std::make_pair( std::set<std::string>({"adr"}), this->functionSpace() );
    auto fieldNamesWithSpaces = hana::make_tuple( fieldNamesWithSpaceADR );
    M_measurePointsEvaluation = std::make_shared<measure_points_evaluation_type>( fieldNamesWithSpaces );
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( this->keyword() ) )
        M_measurePointsEvaluation->init( evalPoints );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->initPostProcessExportsAndMeasures();
    this->createPostProcessExporters();
    this->createPostProcessMeasures();
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    auto mfields = this->modelFields();
    this->exportResults( time, this->symbolsExpr( mfields ), mfields, this->allMeasuresQuantities() );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::setDoExport( bool b )
{
    if( M_exporter )
        M_exporter->setDoExport( b );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::addMarkerInflowBC( std::string const& markerName )
{
    if( std::find(this->M_bcInflowMarkers.begin(), this->M_bcInflowMarkers.end(), markerName) == this->M_bcInflowMarkers.end() )
        this->M_bcInflowMarkers.push_back( markerName );
}

} // namespace FeelModels
} // namespace Feel
