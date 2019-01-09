/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/multifluid/multifluid.hpp>

#include <feel/feelmodels/levelset/levelsetheavisideexpr.hpp>
#include <feel/feelmodels/levelset/levelsetdeltaexpr.hpp>

#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

#include <feel/feelfilters/loadmesh.hpp>

namespace Feel {
namespace FeelModels {

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
MULTIFLUID_CLASS_TEMPLATE_TYPE::MultiFluid(
        std::string const& prefix,
        worldcomm_ptr_t const& wc,
        std::string const& subPrefix,
        ModelBaseRepository const& modelRep )
: 
    super_type( prefixvm(prefix,"fluid"), false, wc, subPrefix, modelRep ),
    M_prefix( prefix ),
    M_useLagrangeP1iso( false ),
    M_doUpdateGlobalLevelset( true ),
    M_doRebuildMatrixVector( false ),
    M_usePicardIterations( false ),
    M_hasInterfaceForcesModel( false ),
    M_enableInextensibility( false ),
    M_hasInextensibilityLM( false ),
    M_doUpdateInextensibilityLM( false )
{
    //-----------------------------------------------------------------------------//
    // Load parameters
    this->loadParametersFromOptionsVm();
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
typename MULTIFLUID_CLASS_TEMPLATE_TYPE::self_ptrtype
MULTIFLUID_CLASS_TEMPLATE_TYPE::New(
        std::string const& prefix,
        worldcomm_ptr_t const& wc,
        std::string const& subPrefix,
        ModelBaseRepository const& modelRep )
{
    self_ptrtype new_multifluid( new self_type(prefix, wc, subPrefix, modelRep ) );
    return new_multifluid;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
std::string
MULTIFLUID_CLASS_TEMPLATE_TYPE::expandStringFromSpec( std::string const& s )
{
    std::string res = s;
    res = fluid_type::expandStringFromSpec( res );
    return res;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
typename MULTIFLUID_CLASS_TEMPLATE_TYPE::mesh_ptrtype
MULTIFLUID_CLASS_TEMPLATE_TYPE::createMesh()
{
    mesh_ptrtype mesh;

    this->log("MultiFluid","createMesh", "start");
    this->timerTool("Constructor").start();

    std::string fmpath = (fs::path(this->rootRepository()) / fs::path(this->fileNameMeshPath())).string();
    if (this->doRestart())
    {
        this->log("createMesh","", "restart with : "+fmpath);

        if ( !this->restartPath().empty() )
        {
            fmpath = (fs::path( this->restartPath()) / fs::path(this->fileNameMeshPath())).string();
        }
        // reload mesh path stored in file
        std::ifstream file( fmpath.c_str() );
        if ( !file )
            CHECK( false ) << "Fail to open the txt file containing path of msh file : " << fmpath << "\n";
        std::string mshfile;
        if ( ! ( file >> mshfile ) )
            CHECK( false ) << "Fail to read the msh path in file : " << fmpath << "\n";
        file.close();

        mesh = loadMesh(
                _mesh=new mesh_type( this->worldCommPtr() ),
                _filename=mshfile,
                _worldcomm=this->worldCommPtr(),
                //_prefix=this->prefix(),
                _rebuild_partitions=false,
                _savehdf5=0,
                _straighten=1,
                _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES
                );
    }
    else
    {
        if (this->hasMeshFile())
        {
            std::string mshfileRebuildPartitions = this->rootRepository() + "/" + this->prefix() + ".msh";

            this->log("createMesh","", "load mesh file : " + this->meshFile());
            std::string meshFileExt = fs::path( this->meshFile() ).extension().string();
            bool rebuildPartition = boption(_prefix=this->prefix(), _name="gmsh.partition");
            if ( rebuildPartition && meshFileExt != ".msh" )
                CHECK( false ) << "Can not rebuild at this time the mesh partitionining with other format than .msh : TODO";

            mesh = loadMesh(
                    _mesh=new mesh_type( this->worldCommPtr() ),
                    _filename=this->meshFile(),
                    _worldcomm=this->worldCommPtr(),
                    _prefix=this->prefix(),
                    _rebuild_partitions=rebuildPartition,
                    _rebuild_partitions_filename=mshfileRebuildPartitions,
                    _partitions=this->worldComm().localSize(),
                    _savehdf5=0,
                    _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES
                    );

            if (rebuildPartition) this->setMeshFile(mshfileRebuildPartitions);
        }
        else if (this->hasGeoFile() )
        {
            std::string mshfile = this->rootRepository() + "/" + this->prefix() + ".msh";
            this->setMeshFile(mshfile);

            gmsh_ptrtype geodesc = geo( 
                    _filename=this->geoFile(),
                    _prefix=this->prefix(),
                    _worldcomm=this->worldCommPtr() 
                    );
            // allow to have a geo and msh file with a filename equal to prefix
            geodesc->setPrefix(this->prefix());
            mesh = createGMSHMesh(
                    _mesh=new mesh_type,_desc=geodesc,
                    _prefix=this->prefix(),_worldcomm=this->worldCommPtr(),
                    _partitions=this->worldComm().localSize(),
                    _directory=this->rootRepository() 
                    );
        }
        this->saveMeshFile(fmpath);
    }

    CHECK( mesh ) << "mesh generation fail";

    //M_fluid->setMshfileStr( this->mshfileStr() );

    double tElapsed = this->timerTool("Constructor").stop("createMesh");
    this->log("MultiFluid","createMesh", (boost::format("finish in %1% s") %tElapsed).str() );

    return mesh;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::init()
{
    CHECK( M_nFluids >= 2 ) << "Multifluid must contain at least 2 fluids.\n";
    this->log("MultiFluid", "init", "start");

    // Read mesh info from multifluid options
    if ( Environment::vm().count( prefixvm(this->prefix(),"mesh.filename").c_str() ) )
    {
        std::string meshfile = Environment::expand( soption(_prefix=this->prefix(),_name="mesh.filename") );
        if ( fs::path( meshfile ).extension() == ".geo" )
            this->setGeoFile( meshfile );
        else
            this->setMeshFile( meshfile );
    }

    // Build inherited FluidMechanics
    this->loadMesh( this->createMesh() );

    // Init inherited FluidMechanics (to build spaces, algebraic data, ...)
    // but don't build algebraic data
    super_type::init( false, false );
    // Init LevelSets
    this->createLevelsets();
    // Init post-process
    this->initPostProcess();
    // Update inextensibility LM if needed
    if( this->M_hasInextensibilityLM )
        this->updateInextensibilityLM();
    // Finally build algebraic data
    this->buildBlockVector();
    this->initAlgebraicFactory();
    // Do not request matrix/vector update since they were just updated
    M_doRebuildMatrixVector = false;

    // "Deep" copy FluidMechanics materialProperties
    M_fluidMaterialProperties.reset( new material_properties_type( this->fluidPrefix() ) );
    // Create M_interfaceForces
    M_interfaceForces.reset( new element_levelset_vectorial_type(this->functionSpaceLevelsetVectorial(), "InterfaceForces") ); 

    // Init FluidMechanics materialProperties
    M_fluidMaterialProperties->updateForUse( this->materialProperties()->dynamicViscositySpace(), this->modelProperties().materials() );
    // Init levelsets materialProperties and interfaceForcesModels
    for( auto const& lsMaterialProperties: M_levelsetsMaterialProperties )
    {
        lsMaterialProperties.second->updateForUse(this->materialProperties()->dynamicViscositySpace(), M_levelsets[lsMaterialProperties.first]->modelProperties().materials() );
    }
    // Update density-viscosity
    this->updateFluidDensityViscosity();

    this->log("MultiFluid", "init", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    M_nFluids = ioption( _name="nfluids", _prefix=this->prefix() );

    M_usePicardIterations = boption( _name="use-picard-iterations", _prefix=this->prefix() );

    M_hasInterfaceForcesModel = false;

    M_useLagrangeP1iso = boption( _name="use-ls-P1iso-mesh", _prefix=this->prefix() );

    // Global levelset parameters
    // TODO

    uint16_type nLevelSets = M_nFluids - 1;

    M_enableInextensibility = false;
    for( uint16_type n = 0; n < nLevelSets; ++n )
    {
        std::string const lsName = levelsetName(n);
        auto levelset_prefix = prefixvm(this->prefix(), lsName);

        M_hasInextensibility[lsName] = boption( _name="enable-inextensibility", _prefix=levelset_prefix );
        if( M_hasInextensibility[lsName] ) M_enableInextensibility = true;

        M_inextensibilityMethods[lsName] = soption( _name="inextensibility-method", _prefix=levelset_prefix );
        CHECK( M_inextensibilityMethods[lsName] == "penalty" || M_inextensibilityMethods[lsName] == "lagrange-multiplier" ) 
            << "invalid inextensiblity-method " << M_inextensibilityMethods[lsName]
            << ", should be \"penalty\" or \"lagrange-multiplier\"" << std::endl;

        if( M_hasInextensibility[lsName] && M_inextensibilityMethods[lsName] == "lagrange-multiplier" )
            M_hasInextensibilityLM = true;

        M_inextensibilityGamma[lsName] = doption( _name="inextensibility-gamma", _prefix=levelset_prefix );
    }

    for( uint16_type n = 0; n < nLevelSets; ++n )
    {
        std::string const lsName = levelsetName(n);
        auto levelset_prefix = prefixvm(this->prefix(), lsName);
        M_levelsetReinitEvery[lsName] = ioption( _name="reinit-every", _prefix=levelset_prefix );
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
MULTIFLUID_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||--------------Info : MultiFluid---------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix           : " << this->prefix()
           << "\n   Root Repository  : " << this->rootRepository()
           << "\n   Number of fluids : " << M_nFluids;

    *_ostr << "\n   Fluids Parameters";
    *_ostr << "\n     -- fluid 0 (outer fluid)"
           << "\n       * rho : " << this->M_fluidMaterialProperties->cstDensity()
           << "\n       * mu  : " << this->M_fluidMaterialProperties->cstMu();
    *_ostr << this->M_fluidMaterialProperties->getInfo()->str();
    for( auto const& lsMaterialProperty: M_levelsetsMaterialProperties )
    {
    *_ostr << "\n     -- fluid " << "\"" << lsMaterialProperty.first << "\""
           << "\n       * rho : " << lsMaterialProperty.second->cstDensity()
           << "\n       * mu  : " << lsMaterialProperty.second->cstMu();
    *_ostr << lsMaterialProperty.second->getInfo()->str();
    }

    *_ostr << "\n   Level Sets Parameters";
    for( auto const& lsReinitEvery: M_levelsetReinitEvery )
    {
    *_ostr << "\n     -- level set " << "\"" << lsReinitEvery.first << "\""
           << "\n       * reinit every : " << lsReinitEvery.second;
    }

    *_ostr << "\n   Forces Parameters";
    *_ostr << "\n     -- has interface forces : " << std::boolalpha << this->M_hasInterfaceForcesModel;
    if( this->M_hasInterfaceForcesModel )
    {
        for( auto const& lsInterfaceForces: M_levelsetInterfaceForcesModels )
        {
    *_ostr << "\n     -- level set " << "\"" << lsInterfaceForces.first << "\"";
            for( auto const& force: lsInterfaceForces.second )
    *_ostr << "\n       * force model : " << force.second->getInfo()->str();
        }
    }

    *_ostr << "\n";
    *_ostr << super_type::getInfo()->str();
    for( auto const& levelset: M_levelsets )
    {
    *_ostr << levelset.second->getInfo()->str();
    }

    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n\n";

    return _ostr;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
typename MULTIFLUID_CLASS_TEMPLATE_TYPE::space_inextensibilitylm_ptrtype const&
MULTIFLUID_CLASS_TEMPLATE_TYPE::functionSpaceInextensibilityLM() const
{
    CHECK( !M_doUpdateInextensibilityLM ) << "updateInextensibilityLM() must be called before using functionSpaceInextensibilityLM()\n";
    return M_spaceInextensibilityLM;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
typename MULTIFLUID_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
MULTIFLUID_CLASS_TEMPLATE_TYPE::globalLevelsetElt() const
{
    if( !M_globalLevelsetElt )
        M_globalLevelsetElt.reset( new element_levelset_type(this->functionSpaceLevelset(), "GlobalLevelset") );
    if( M_doUpdateGlobalLevelset )
    {
        M_globalLevelsetElt->on( 
                _range=elements( M_globalLevelsetElt->mesh() ),
                _expr=this->globalLevelsetExpr()
                );
        M_doUpdateGlobalLevelset = false;
    }
    return M_globalLevelsetElt;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
int
MULTIFLUID_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    if( this->hasInextensibilityLM() )
        return super_type::nBlockMatrixGraph() + 1;
    else
        return super_type::nBlockMatrixGraph();
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
MULTIFLUID_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("MultiFluid","buildBlockMatrixGraph", "start" );
    BlocksBaseGraphCSR myBlockGraph = super_type::buildBlockMatrixGraph();

    int indexBlock = super_type::nBlockMatrixGraph();

    if( this->hasInextensibilityLM() )
    {
        this->log("MultiFluid","buildBlockMatrixGraph", "start build lagrange-multiplier" );
        // Matrix stencil
        BlocksStencilPattern patCouplingLM(1, super_type::space_fluid_type::nSpaces,size_type(Pattern::ZERO));
        patCouplingLM(0,0) = size_type(Pattern::COUPLED);

        myBlockGraph(indexBlock,0) = stencil(_test=this->functionSpaceInextensibilityLM(), _trial=this->functionSpace(),
                                             _pattern_block=patCouplingLM,
                                             _diag_is_nonzero=false,_close=false)->graph();
        myBlockGraph(0,indexBlock) = stencil(_test=this->functionSpace(), _trial=this->functionSpaceInextensibilityLM(),
                                             _pattern_block=patCouplingLM.transpose(),
                                             _diag_is_nonzero=false,_close=false)->graph();
        ++indexBlock;
    }

    this->log("MultiFluid","buildBlockMatrixGraph", "finish" );

    return myBlockGraph;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
size_type
MULTIFLUID_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    auto res = super_type::nLocalDof();
    if( this->hasInextensibilityLM() )
    {
        res += this->functionSpaceInextensibilityLM()->nLocalDofWithGhost();
    }
    return res;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
bool
MULTIFLUID_CLASS_TEMPLATE_TYPE::hasInextensibilityLM() const
{
    bool hasInextensibilityLM = false;
    if( this->M_enableInextensibility )
    {
        for( auto const& inextmethod: M_inextensibilityMethods )
        {
            if( inextmethod.second == "lagrange-multiplier" )
            {
                hasInextensibilityLM = true;
                break;
            }
        }
    }
    return hasInextensibilityLM;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateInextensibilityLM()
{
    this->log("MultiFluid", "updateInextensibilityLM", "start");
    M_inextensibleLevelsets.clear();
    for( auto const& ls: M_levelsets )
    {
        if( this->hasInextensibility(ls.first) )
        {
            if( this->inextensibilityMethod(ls.first) == "lagrange-multiplier" )
            {
                M_inextensibleLevelsets.push_back( ls.second->phi() );
            }
        }
    }
    // Compute inextensible levelsets elements range
    auto inextensibleLevelsetsExpr = Feel::FeelModels::globalLevelsetExpr( M_inextensibleLevelsets );
    auto diracls = vf::project(
            _space=this->M_levelsetSpaceManager->functionSpaceScalar(), _range=this->M_levelsetSpaceManager->rangeMeshElements(),
            _expr=Feel::FeelModels::levelsetDelta( inextensibleLevelsetsExpr, cst(M_globalLevelsetThicknessInterface) )
            );
    auto dirac = vf::project( 
            _space=this->functionSpacePressure(), _range=this->rangeMeshElements(),
            _expr=idv(diracls)
            );
    auto it_elt = this->mesh()->beginOrderedElement();
    auto en_elt = this->mesh()->endOrderedElement();

    const rank_type pid = this->mesh()->worldCommElements().localRank();
    const int ndofv = super_type::space_fluid_pressure_type::fe_type::nDof;
    elements_reference_wrapper_ptrtype diracElts( new elements_reference_wrapper_type );

    for (; it_elt!=en_elt; it_elt++)
    {
        auto const& elt = boost::unwrap_ref( *it_elt );
        if ( elt.processId() != pid )
            continue;
        bool mark_elt = true;
        for (int j=0; j<ndofv; j++)
        {
            if ( dirac.localToGlobal(elt.id(), j, 0) < 1e-6 )
            {
                mark_elt = false;
                break; //don't need to do the others dof
            }
        }
        if( mark_elt )
            diracElts->push_back( boost::cref(elt) );
    }

    M_rangeInextensibilityLM = boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
            diracElts->begin(), diracElts->end(), diracElts );

    // Lagrange-multiplier inextensibility space
    M_spaceInextensibilityLM = space_inextensibilitylm_type::New(
            _mesh=this->mesh(),
            _range=M_rangeInextensibilityLM,
            _worldscomm=this->localNonCompositeWorldsComm()
            );

    M_doUpdateInextensibilityLM = false;
    M_doRebuildMatrixVector = true;
    this->log("MultiFluid", "updateInextensibilityLM", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
bool
MULTIFLUID_CLASS_TEMPLATE_TYPE::hasInterfaceForces() const
{
    return M_hasInterfaceForcesModel || (M_additionalInterfaceForcesModel.size() > 0);
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::addInterfaceForce( interfaceforces_model_ptrtype model, std::string const& name )
{
    std::string forceName;
    if( name.empty() )
    {
        forceName = ( boost::format("force%1%" ) %(M_additionalInterfaceForcesModel.size()) ).str();
    }
    else
    {
        CHECK( M_additionalInterfaceForcesModel.find(name) == M_additionalInterfaceForcesModel.end() ) 
            << "Multifluid already has an interface force model named " << name << std::endl;
        forceName = name;
    }

    M_additionalInterfaceForcesModel[forceName] = model;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
typename MULTIFLUID_CLASS_TEMPLATE_TYPE::interfaceforces_model_ptrtype const&
MULTIFLUID_CLASS_TEMPLATE_TYPE::interfaceForce( std::string const& name ) const
{
    auto it = M_additionalInterfaceForcesModel.find(name);
    CHECK( it != M_additionalInterfaceForcesModel.end() ) << "no force named " << name << std::endl;
    return it->second;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
std::map<std::string, typename MULTIFLUID_CLASS_TEMPLATE_TYPE::interfaceforces_model_ptrtype> const&
MULTIFLUID_CLASS_TEMPLATE_TYPE::interfaceForces() const
{ 
    // TODO : get more general access to forces
    return M_additionalInterfaceForcesModel; 
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("MultiFluid", "solve", "start");
    this->timerTool("Solve").start();

    this->solveImpl();
    if( M_usePicardIterations )
    {
        double errorVelocityL2 = 0.;
        int picardIter = 0;
        auto u_old = this->fluidModel()->functionSpaceVelocity()->element();
        do
        {
            picardIter++;
            u_old = this->fieldVelocity();
            this->solveImpl();
            auto u = this->fieldVelocity();

            double uOldL2Norm = integrate(
                    _range=elements(this->mesh()),
                    _expr=trans(idv(u_old))*idv(u_old)
                    ).evaluate()(0,0);
            errorVelocityL2 = integrate(
                    _range=elements(this->mesh()),
                    _expr=trans(idv(u)-idv(u_old))*(idv(u)-idv(u_old))
                    ).evaluate()(0,0);
            errorVelocityL2 = std::sqrt(errorVelocityL2) / std::sqrt(uOldL2Norm);

            Feel::cout << "Picard iteration " << picardIter << ": errorVelocityL2 = " << errorVelocityL2 << std::endl;
        } while( errorVelocityL2 > 0.01);
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","solve","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::solveImpl()
{
    // Update density and viscosity
    this->updateFluidDensityViscosity();
    // Update inextensibility LM if needed
    if( M_doUpdateInextensibilityLM )
        this->updateInextensibilityLM();
    // Solve fluid equations (with direct assembly of interface forces)
    this->solveFluid();
    // Advect levelsets
    this->advectLevelsets();
    // Reinitialize
    for( auto const& ls: M_levelsets )
    {
        if( M_levelsetReinitEvery.at(ls.first) > 0 
                && (ls.second->iterSinceReinit()+1) % M_levelsetReinitEvery.at(ls.first) == 0 )
            ls.second->reinitialize();
    }

    // Update global levelset
    M_doUpdateGlobalLevelset = true;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateTime( double time )
{
    // Levelsets
    for( auto const& ls: M_levelsets)
    {
        ls.second->updateTime(time);
    }
    // This
    super_type::updateTime(time);
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->log("MultiFluid", "updateTimeStep", "start");
    // Fluid
    super_type::updateTimeStep();
    // Levelsets
    for( auto const& ls: M_levelsets)
    {
        ls.second->updateTimeStep();
    }

    this->log("MultiFluid", "updateTimeStep", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::exportResultsImpl( double time )
{
    this->log("MultiFluid","exportResults", "start");

    // Export forces
    for( auto const& lsInterfaceForces: M_levelsetInterfaceForcesModels )
    {
        std::string const lsName = lsInterfaceForces.first;
        auto levelset_prefix = prefixvm(this->prefix(), lsName);
        for( auto const& force: lsInterfaceForces.second )
        {
            this->M_exporter->step(time)->add( prefixvm(levelset_prefix, force.first),
                    prefixvm(levelset_prefix, force.first),
                    force.second->lastInterfaceForce() );
        }
    }
    for( auto const& force: M_additionalInterfaceForcesModel )
    {
        this->M_exporter->step(time)->add( prefixvm("additional", force.first),
                prefixvm("additional", force.first),
                force.second->lastInterfaceForce() );
    }

    // Export measures
    this->exportMeasures( time );
    // Export fluid
    super_type::exportResults(time);
    // Export global levelsets
    if( this->nLevelsets() > 1 )
    {
        if( !M_globalLevelsetExporter )
        {
            M_globalLevelsetExporter = exporter(
                    _mesh=M_levelsetSpaceManager->mesh(),
                    _name="ExportLS",
                    _geo="static",
                    _path=this->exporterPath()
                    );
        }
        M_globalLevelsetExporter->step( time )->add(
                prefixvm(this->prefix(),"GlobalLevelset.Phi"),
                prefixvm(this->prefix(),prefixvm(this->subPrefix(),"GlobalLevelset.Phi")),
                *this->globalLevelsetElt()
                );
        M_globalLevelsetExporter->save();
    }
    for( auto const& ls: M_levelsets)
    {
        ls.second->exportResults(time);
    }

    this->log("MultiFluid", "exportResults", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::exportMeasures( double time )
{
    // Levelset forces
    for( std::string const& levelsetName: M_postProcessMeasuresLevelsetForces )
    {
        auto measuredLevelsetForce = this->computeLevelsetForce( levelsetName );
        std::vector<double> vecMeasuredLevelsetForce = { measuredLevelsetForce(0,0) };
        if( nDim > 1 ) vecMeasuredLevelsetForce.push_back( measuredLevelsetForce(1,0) );
        if( nDim > 2 ) vecMeasuredLevelsetForce.push_back( measuredLevelsetForce(2,0) );
        this->postProcessMeasuresIO().setMeasureComp( levelsetName + ".force", vecMeasuredLevelsetForce );
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
typename MULTIFLUID_CLASS_TEMPLATE_TYPE::force_type
MULTIFLUID_CLASS_TEMPLATE_TYPE::computeLevelsetForce( std::string const& name ) const
{
    auto solFluid = this->fieldVelocityPressurePtr();
    auto u = solFluid->template element<0>();
    auto p = solFluid->template element<1>();
    std::string matName = this->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
    auto sigmav = Feel::FeelModels::fluidMecNewtonianStressTensor<2*super_type::nOrderVelocity>(u,p,*this->materialProperties(),matName,true);

    auto N = this->levelsetModel(name)->N();
    auto D = this->levelsetModel(name)->D();
    return integrate(_range=elements(this->mesh()),
                     _expr= - sigmav * idv(N) * idv(D),
                     _geomap=this->geomap() ).evaluate();
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::createLevelsets()
{
    this->log("MultiFluid", "createLevelsets", "start");
    // Get levelset mesh
    mesh_ptrtype mesh;
    if( this->M_useLagrangeP1iso )
    {
        // Build Lagrange P1 iso-U mesh and build levelsets with it
        M_opLagrangeP1iso = lagrangeP1( this->functionSpaceVelocity()->compSpace() );
        mesh = this->M_opLagrangeP1iso->mesh(); 
    }
    else
    {
        // Else build from common mesh
        mesh = this->mesh();
    }
    // Build levelsets space manager
    M_levelsetSpaceManager = std::make_shared<levelset_space_manager_type>( mesh, this->globalLevelsetPrefix() );
    // Temporary hack: ensures the defaults levelset spaces are built for interfaceForces
    M_levelsetSpaceManager->createFunctionSpaceDefault();
    // Build levelsets tool manager
    M_levelsetToolManager = std::make_shared<levelset_tool_manager_type>( M_levelsetSpaceManager, this->globalLevelsetPrefix() );
    // Update global levelset thickness interface
    if( Environment::vm().count( prefixvm(this->globalLevelsetPrefix(),"thickness-interface").c_str() ) )
        M_globalLevelsetThicknessInterface = doption( _name="thickness-interface", _prefix=this->globalLevelsetPrefix() );
    else
        M_globalLevelsetThicknessInterface = 1.5 * M_levelsetSpaceManager->mesh()->hAverage();

    // Build levelsets exporter
    //M_globalLevelsetExporter = exporter(
            //_mesh=M_levelsetSpaceManager->mesh(),
            //_name="ExportLS",
            //_geo="static",
            //_path=this->exporterPath()
            //);
    // Build levelsets
    uint16_type nLevelsets = M_nFluids - 1;
    for( uint16_type i = 0; i < nLevelsets; ++i )
    {
        std::string const lsName = levelsetName(i);
        auto levelset_prefix = prefixvm(this->prefix(), lsName);
        auto & levelset = M_levelsets[lsName];
        levelset.reset(
                new levelset_type( levelset_prefix, this->worldCommPtr(), "", this->repository().rootWithoutNumProc() )
                );
        levelset->setFunctionSpaceManager( M_levelsetSpaceManager );
        levelset->setToolManager( M_levelsetToolManager );
        hana::eval_if( std::is_same<space_levelset_advection_velocity_type, space_fluid_velocity_type>{},
                    [&](auto _) {
                        if( !M_useLagrangeP1iso && !_(levelset)->useSpaceIsoPN() )
                        {
                            Feel::cout << "Using fluid velocity function space\n";
                            _(levelset)->setFunctionSpaceAdvectionVelocity( this->functionSpaceVelocity() );
                        }
                    },
                    [&] {}
                );
        // Set global options if unspecified otherwise
        if( !Environment::vm().count( prefixvm(levelset_prefix,"thickness-interface").c_str() ) )
            levelset->setThicknessInterface( this->globalLevelsetThicknessInterface() );
        // Initialize LevelSets
        levelset->init();

        // Build levelsets materialProperties
        M_levelsetsMaterialProperties[lsName].reset(
                new material_properties_type( levelset_prefix )
                );
        // Build levelsets interfaceForcesModels
        if( Environment::vm().count( prefixvm(levelset_prefix, "interface-forces-model").c_str() ) )
        {
            std::vector<std::string> interfaceForcesModels = Environment::vm()[prefixvm(levelset_prefix, "interface-forces-model").c_str()].template as<std::vector<std::string>>();
            // Remove (unwanted) duplicates
            std::sort( interfaceForcesModels.begin(), interfaceForcesModels.end() );
            interfaceForcesModels.erase( std::unique(interfaceForcesModels.begin(), interfaceForcesModels.end()), interfaceForcesModels.end() );

            for( uint16_type n = 0; n < interfaceForcesModels.size(); ++n )
            {
                std::string const forceName = interfaceForcesModels[n];
                M_levelsetInterfaceForcesModels[lsName][forceName] = interfaceforces_factory_type::instance().createObject( 
                        forceName
                        );
                M_levelsetInterfaceForcesModels[lsName][forceName]->build( levelset_prefix, M_levelsets[lsName], this->fluidModel() );
            }

            M_hasInterfaceForcesModel = true;
        }
    }

    this->log("MultiFluid", "createLevelsets", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    std::string modelName = "multifluid";
    pt::ptree ptree = this->modelProperties().postProcess().pTree( modelName );
    std::string ppTypeMeasures = "Measures";
    for( auto const& ptreeLevel0 : ptree )
    {
        std::string ptreeLevel0Name = ptreeLevel0.first;
        if ( ptreeLevel0Name != ppTypeMeasures ) continue;
        for( auto const& ptreeLevel1 : ptreeLevel0.second )
        {
            std::string ptreeLevel1Name = ptreeLevel1.first;
            if ( ptreeLevel1Name == "LevelsetForces" )
            {
                // get list of marker
                std::set<std::string> levelsetNames;
                std::string levelsetNameUnique = ptreeLevel1.second.template get_value<std::string>();
                if ( levelsetNameUnique.empty() )
                {
                    for (auto const& ptreeMarker : ptreeLevel1.second )
                    {
                        std::string levelsetName = ptreeMarker.second.template get_value<std::string>();
                        levelsetNames.insert( levelsetName );
                    }
                }
                else
                {
                    levelsetNames.insert( levelsetNameUnique );
                }
                // save forces measure for each levelset
                for ( std::string const& levelsetName : levelsetNames )
                {
                    M_postProcessMeasuresLevelsetForces.push_back( levelsetName );
                }
            }
        }
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
size_type
MULTIFLUID_CLASS_TEMPLATE_TYPE::initStartBlockIndexFieldsInMatrix()
{
    size_type currentStartIndex = super_type::initStartBlockIndexFieldsInMatrix();

    if( this->hasInextensibilityLM() )
    {
        // Add inextensibility LM block index
        this->setStartSubBlockSpaceIndex( "inextensibility-lm", currentStartIndex++ );
    }

    return currentStartIndex;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
int
MULTIFLUID_CLASS_TEMPLATE_TYPE::initBlockVector()
{
    int currentBlockIndex = super_type::initBlockVector();
    if( this->hasInextensibilityLM() )
    {
        this->M_blockVectorSolution(currentBlockIndex) = this->backend()->newVector( this->functionSpaceInextensibilityLM() );
        ++currentBlockIndex;
    }

    return currentBlockIndex;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateFluidDensityViscosity()
{
    this->log("MultiFluid", "updateFluidDensityViscosity", "start");
    this->timerTool("Solve").start();

    auto globalH = Feel::FeelModels::levelsetHeaviside( 
            this->globalLevelsetExpr(),
            cst( this->globalLevelsetThicknessInterface() )
            );

    auto rho = vf::project( 
            this->materialProperties()->dynamicViscositySpace(),
            elements(this->mesh()),
            idv(M_fluidMaterialProperties->fieldRho())*globalH
            );

    auto mu = vf::project( 
            this->materialProperties()->dynamicViscositySpace(),
            elements(this->mesh()),
            idv(M_fluidMaterialProperties->fieldMu())*globalH
            );

    for( auto const& ls: M_levelsets )
    {
        auto Hi = Feel::FeelModels::levelsetHeaviside( 
                idv( ls.second->phi() ),
                cst( this->globalLevelsetThicknessInterface() )
                );
        rho += vf::project( 
                this->materialProperties()->dynamicViscositySpace(),
                elements(this->mesh()),
                idv(M_levelsetsMaterialProperties[ls.first]->fieldRho())*(1. - Hi)
                );
        mu += vf::project( 
                this->materialProperties()->dynamicViscositySpace(),
                elements(this->mesh()),
                idv(M_levelsetsMaterialProperties[ls.first]->fieldMu())*(1. - Hi)
                );
    }

    this->updateRho( idv(rho) );
    this->updateMu( idv(mu) );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log( "MultiFluid", "updateFluidDensityViscosity", 
            "fluid density/viscosity update in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateInterfaceForces()
{
    this->log("MultiFluid", "updateInterfaceForces", "start");
    this->timerTool("Solve").start();

    M_interfaceForces->zero();


    if( M_hasInterfaceForcesModel )
    {
        this->timerTool("Solve").start();
        for( auto const& lsInterfaceForces: M_levelsetInterfaceForcesModels )
        {
            for( auto const& force: lsInterfaceForces.second )
            {
                if( force.second )
                {
                    force.second->updateInterfaceForces( M_interfaceForces, false );
                }
            }
        }

        double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
        this->log("MultiFluid", "updateInterfaceForces", "update interface (model) forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
    }

    if( M_additionalInterfaceForcesModel.size() > 0 )
    {
        this->timerTool("Solve").start();
        for( auto const& f: M_additionalInterfaceForcesModel )
            f.second->updateInterfaceForces( M_interfaceForces, false );

        double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
        this->log("MultiFluid", "updateInterfaceForces", "update additional interface forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
    }

    this->updateSourceAdded( idv(M_interfaceForces) );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log( "MultiFluid", "updateInterfaceForces", 
            "interface forces updated in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::solveFluid()
{
    this->log("MultiFluid", "solveFluid", "start");
    this->timerTool("Solve").start();

    if( this->M_doRebuildMatrixVector )
    {
        // Rebuild algebraic matrix and vector
        auto graph = this->buildMatrixGraph();
        //this->algebraicFactory()->backend()->clear();
        //this->algebraicFactory()->rebuildMatrixVector( graph, graph->mapRow().indexSplit() );
        this->M_backend = backend( 
                _kind=soption( _name="backend", _prefix=this->fluidPrefix() ), 
                _name=this->fluidPrefix(),
                _rebuild=true,
                _worldcomm=this->worldCommPtr() 
                );
        //this->algebraicFactory()->reset( this->M_backend, graph, graph->mapRow().indexSplit() );
        this->initAlgebraicFactory();
        // Rebuild solution vector
        this->buildBlockVector();

        M_doRebuildMatrixVector = false;
    }

    super_type::solve();

    double timeElapsed = this->timerTool("Solve").stop();
    this->log( "MultiFluid", "solveFluid", 
            "fluid problem solved in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::advectLevelsets()
{
    this->log("MultiFluid", "advectLevelsets", "start");
    this->timerTool("Solve").start();

    auto const& u = this->fieldVelocity();
    
    for( auto const& ls: M_levelsets )
    {
        ls.second->advect( idv(u) );
    }

    if( this->hasInextensibilityLM() )
    {
        M_doUpdateInextensibilityLM = true;
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log( "MultiFluid", "advectLevelsets", 
            "level-sets advection done in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEAdditional( DataUpdateLinear & data ) const
{
    bool _BuildCstPart = data.buildCstPart();
    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("MultiFluid","updateLinearPDEAdditional", "start"+sc );
    this->timerTool("Solve").start();

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

    // Update interface forces
    if( this->hasInterfaceForces() )
    {
        this->log("MultiFluid", "updateLinearPDEAdditional", "start: update interface forces");
        this->timerTool("Solve").start();
        if( M_hasInterfaceForcesModel )
        {
            this->timerTool("Solve").start();
            for( auto const& lsInterfaceForces: M_levelsetInterfaceForcesModels )
            {
                for( auto const& force: lsInterfaceForces.second )
                {
                    if( force.second )
                    {
                        force.second->updateFluidInterfaceForcesLinearPDE( data );
                    }
                }
            }

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateLinearPDEAdditional", "update interface (model) forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        if( M_additionalInterfaceForcesModel.size() > 0 )
        {
            this->timerTool("Solve").start();
            for( auto const& f: M_additionalInterfaceForcesModel )
                f.second->updateFluidInterfaceForcesLinearPDE( data );

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateLinearPDEAdditional", "update additional interface forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        double timeElapsed = this->timerTool("Solve").stop();
        this->log( "MultiFluid", "updateLinearPDEAdditional", 
                "interface forces updated in "+(boost::format("%1% s") %timeElapsed).str() );
    }

    // Update inextensibility
    if( this->M_enableInextensibility )
    {
        auto mesh = this->mesh();
        auto Xh = this->functionSpace();

        auto const& U = this->fieldVelocityPressure();
        auto u = U.template element<0>();
        auto v = U.template element<0>();
        auto Id = vf::Id<nDim, nDim>();

        auto rowStartInMatrix = this->rowStartInMatrix();
        auto colStartInMatrix = this->colStartInMatrix();
        auto rowStartInVector = this->rowStartInVector();
        auto bilinearForm_PatternDefault = form2( 
                _test=Xh,_trial=Xh,_matrix=A,
                _pattern=size_type(Pattern::DEFAULT),
                _rowstart=rowStartInMatrix,
                _colstart=colStartInMatrix 
                );

        for( auto const& ls: M_levelsets )
        {
            std::string const& lsName = ls.first;
            if( this->hasInextensibility(lsName) )
            {
                auto N = ls.second->N();
                auto NxN = idv(N)*trans(idv(N));

                if( this->inextensibilityMethod(lsName) == "penalty" )
                {
                    this->timerTool("Solve").start();

                    if( BuildNonCstPart )
                    {
                        auto D = ls.second->D();
                        bilinearForm_PatternDefault += integrate(
                                _range=elements(mesh),
                                _expr=this->M_inextensibilityGamma.at(lsName)*trace((Id-NxN)*gradt(u))*trace((Id-NxN)*grad(v))*idv(D)/h(),
                                _geomap=this->geomap()
                                );
                    }

                    double timeElapsedInextensibility_Penalty = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateLinearPDEAdditional",
                            "assembly inextensibility (penalty) in "+(boost::format("%1% s") 
                                %timeElapsedInextensibility_Penalty).str() );
                }
                if( this->inextensibilityMethod(lsName) == "lagrange-multiplier" )
                {
                    CHECK( this->hasStartSubBlockSpaceIndex("inextensibility-lm") ) << " start dof index for inextensibility-lm is not present\n";
                    this->timerTool("Solve").start();

                    size_type startBlockIndexInextensibilityLM = this->startSubBlockSpaceIndex("inextensibility-lm");
                    auto lambda = this->functionSpaceInextensibilityLM()->element();

                    if( BuildNonCstPart )
                    {
                        auto inextensibleLevelsetsExpr = Feel::FeelModels::globalLevelsetExpr( M_inextensibleLevelsets );
                        auto inextensibleLevelsetsDeltaExpr = Feel::FeelModels::levelsetDelta( inextensibleLevelsetsExpr, cst(M_globalLevelsetThicknessInterface) );
                        auto inextensibleLevelsetsDelta = vf::project(
                                _space=this->M_levelsetSpaceManager->functionSpaceScalar(), 
                                _range=this->M_levelsetSpaceManager->rangeMeshElements(),
                                _expr=inextensibleLevelsetsDeltaExpr
                                );

                        form2( _trial=this->functionSpaceInextensibilityLM(), _test=this->functionSpace(), _matrix=A,
                               _rowstart=rowStartInMatrix,
                               _colstart=colStartInMatrix+startBlockIndexInextensibilityLM ) +=
                            integrate( _range=this->M_rangeInextensibilityLM,
                                       _expr=idt(lambda)*trace((Id-NxN)*grad(v))*idv(inextensibleLevelsetsDelta),
                                       _geomap=this->geomap()
                                       );
                        form2( _trial=this->functionSpace(), _test=this->functionSpaceInextensibilityLM(), _matrix=A,
                               _rowstart=rowStartInMatrix+startBlockIndexInextensibilityLM,
                               _colstart=colStartInMatrix ) +=
                            integrate( _range=this->M_rangeInextensibilityLM,
                                       _expr=id(lambda)*trace((Id-NxN)*gradt(u))*idv(inextensibleLevelsetsDelta),
                                       _geomap=this->geomap()
                                       );
                    }

                    double timeElapsedInextensibility_LagrangeMult = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateLinearPDEAdditional",
                            "assembly inextensibility (lagrange-multiplier) in "+(boost::format("%1% s") 
                                %timeElapsedInextensibility_LagrangeMult).str() );
                }
            }
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","updateLinearPDEAdditional","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobianAdditional( DataUpdateJacobian & data ) const
{
    bool _BuildCstPart = data.buildCstPart();
    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("MultiFluid","updateJacobianAdditional", "start"+sc );
    this->timerTool("Solve").start();

    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

    // Update interface forces
    if( this->hasInterfaceForces() )
    {
        this->log("MultiFluid", "updateJacobianAdditional", "start: update interface forces");
        this->timerTool("Solve").start();
        if( M_hasInterfaceForcesModel )
        {
            this->timerTool("Solve").start();
            for( auto const& lsInterfaceForces: M_levelsetInterfaceForcesModels )
            {
                for( auto const& force: lsInterfaceForces.second )
                {
                    if( force.second )
                    {
                        force.second->updateFluidInterfaceForcesJacobian( data );
                    }
                }
            }

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateJacobianAdditional", "update interface (model) forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        if( M_additionalInterfaceForcesModel.size() > 0 )
        {
            this->timerTool("Solve").start();
            for( auto const& f: M_additionalInterfaceForcesModel )
                f.second->updateFluidInterfaceForcesJacobian( data );

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateJacobianAdditional", "update additional interface forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        double timeElapsed = this->timerTool("Solve").stop();
        this->log( "MultiFluid", "updateJacobianAdditional", 
                "interface forces updated in "+(boost::format("%1% s") %timeElapsed).str() );
    }

    // Update inextensibility
    if( this->M_enableInextensibility )
    {
        auto mesh = this->mesh();
        auto Xh = this->functionSpace();

        auto rowStartInMatrix = this->rowStartInMatrix();
        auto colStartInMatrix = this->colStartInMatrix();
        auto rowStartInVector = this->rowStartInVector();
        auto bilinearForm_PatternCoupled = form2( 
                _test=Xh,_trial=Xh,_matrix=J,
                _pattern=size_type(Pattern::COUPLED),
                _rowstart=rowStartInMatrix,
                _colstart=colStartInMatrix 
                );

        auto U = Xh->element(XVec, rowStartInVector);
        auto u = U.template element<0>();
        auto v = U.template element<0>();
        auto Id = vf::Id<nDim, nDim>();

        for( auto const& ls: M_levelsets )
        {
            std::string const& lsName = ls.first;
            if( this->hasInextensibility(lsName) )
            {
                auto N = ls.second->N();
                auto NxN = idv(N)*trans(idv(N));
                auto D = ls.second->D();

                if( this->inextensibilityMethod(lsName) == "penalty" )
                {
                    this->timerTool("Solve").start();

                    if( BuildNonCstPart )
                    {
                        bilinearForm_PatternCoupled += integrate(
                                _range=elements(mesh),
                                _expr=this->M_inextensibilityGamma.at(lsName)*trace((Id-NxN)*gradt(u))*trace((Id-NxN)*grad(v))*idv(D)/h(),
                                _geomap=this->geomap()
                                );
                    }

                    double timeElapsedInextensibility_Penalty = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateJacobianAdditional",
                            "assembly inextensibility (penalty) in "+(boost::format("%1% s") 
                                %timeElapsedInextensibility_Penalty).str() );
                }
                if( this->inextensibilityMethod(lsName) == "lagrange-multiplier" )
                {
                    CHECK( this->hasStartSubBlockSpaceIndex("inextensibility-lm") ) << " start dof index for inextensibility-lm is not present\n";
                    this->timerTool("Solve").start();

                    size_type startBlockIndexInextensibilityLM = this->startSubBlockSpaceIndex("inextensibility-lm");
                    auto lambda = this->functionSpaceInextensibilityLM()->element();

                    if( BuildNonCstPart )
                    {
                        auto inextensibleLevelsetsExpr = Feel::FeelModels::globalLevelsetExpr( M_inextensibleLevelsets );
                        auto inextensibleLevelsetsDeltaExpr = Feel::FeelModels::levelsetDelta( inextensibleLevelsetsExpr, cst(M_globalLevelsetThicknessInterface) );
                        auto inextensibleLevelsetsDelta = vf::project(
                                _space=this->M_levelsetSpaceManager->functionSpaceScalar(), 
                                _range=this->M_levelsetSpaceManager->rangeMeshElements(),
                                _expr=inextensibleLevelsetsDeltaExpr
                                );

                        form2( _trial=this->functionSpaceInextensibilityLM(), _test=this->functionSpace(), _matrix=J,
                               _rowstart=rowStartInMatrix,
                               _colstart=colStartInMatrix+startBlockIndexInextensibilityLM ) +=
                            integrate( _range=this->M_rangeInextensibilityLM,
                                       _expr=idt(lambda)*trace((Id-NxN)*grad(v))*idv(inextensibleLevelsetsDelta),
                                       _geomap=this->geomap()
                                       );
                        form2( _trial=this->functionSpace(), _test=this->functionSpaceInextensibilityLM(), _matrix=J,
                               _rowstart=rowStartInMatrix+startBlockIndexInextensibilityLM,
                               _colstart=colStartInMatrix ) +=
                            integrate( _range=this->M_rangeInextensibilityLM,
                                       _expr=id(lambda)*trace((Id-NxN)*gradt(u))*idv(inextensibleLevelsetsDelta),
                                       _geomap=this->geomap()
                                       );
                    }

                    double timeElapsedInextensibility_LagrangeMult = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateJacobianAdditional",
                            "assembly inextensibility (lagrange-multiplier) in "+(boost::format("%1% s") 
                                %timeElapsedInextensibility_LagrangeMult).str() );
                }
            }
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","updateJacobianAdditional","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidualAdditional( DataUpdateResidual & data ) const
{
    bool _BuildCstPart = data.buildCstPart();
    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("MultiFluid","updateResidualAdditional", "start"+sc );
    this->timerTool("Solve").start();

    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool BuildCstPart = _BuildCstPart;
    bool BuildNonCstPart = !BuildCstPart;
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    // Update interface forces
    if( this->hasInterfaceForces() )
    {
        this->log("MultiFluid", "updateResidualAdditional", "start: update interface forces");
        this->timerTool("Solve").start();
        if( M_hasInterfaceForcesModel )
        {
            this->timerTool("Solve").start();
            for( auto const& lsInterfaceForces: M_levelsetInterfaceForcesModels )
            {
                for( auto const& force: lsInterfaceForces.second )
                {
                    if( force.second )
                    {
                        force.second->updateFluidInterfaceForcesResidual( data );
                    }
                }
            }

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateResidualAdditional", "update interface (model) forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        if( M_additionalInterfaceForcesModel.size() > 0 )
        {
            this->timerTool("Solve").start();
            for( auto const& f: M_additionalInterfaceForcesModel )
                f.second->updateFluidInterfaceForcesResidual( data );

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateResidualAdditional", "update additional interface forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        double timeElapsed = this->timerTool("Solve").stop();
        this->log( "MultiFluid", "updateResidualAdditional", 
                "interface forces updated in "+(boost::format("%1% s") %timeElapsed).str() );
    }

    // Update inextensibility
    if( this->M_enableInextensibility )
    {
        auto mesh = this->mesh();
        auto Xh = this->functionSpace();

        auto rowStartInMatrix = this->rowStartInMatrix();
        auto colStartInMatrix = this->colStartInMatrix();
        auto rowStartInVector = this->rowStartInVector();
        auto linearForm_PatternCoupled = form1( 
                _test=Xh,_vector=R,
                _pattern=size_type(Pattern::COUPLED),
                _rowstart=rowStartInVector
                );

        auto U = Xh->element(XVec, rowStartInVector);
        auto u = U.template element<0>();
        auto v = U.template element<0>();
        auto Id = vf::Id<nDim, nDim>();

        for( auto const& ls: M_levelsets )
        {
            std::string const& lsName = ls.first;
            if( this->hasInextensibility(lsName) )
            {
                auto N = ls.second->N();
                auto NxN = idv(N)*trans(idv(N));
                auto D = ls.second->D();

                if( this->inextensibilityMethod(lsName) == "penalty" )
                {
                    this->timerTool("Solve").start();

                    if( BuildNonCstPart && !UseJacobianLinearTerms )
                    {
                        linearForm_PatternCoupled += integrate(
                                _range=elements(mesh),
                                _expr=this->M_inextensibilityGamma.at(lsName)*trace((Id-NxN)*gradv(u))*trace((Id-NxN)*grad(v))*idv(D)/h(),
                                _geomap=this->geomap()
                                );
                    }

                    double timeElapsedInextensibility_Penalty = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateResidualAdditional",
                            "assembly inextensibility (penalty) in "+(boost::format("%1% s") 
                                %timeElapsedInextensibility_Penalty).str() );
                }
                if( this->inextensibilityMethod(lsName) == "lagrange-multiplier" )
                {
                    CHECK( this->hasStartSubBlockSpaceIndex("inextensibility-lm") ) << " start dof index for inextensibility-lm is not present\n";
                    this->timerTool("Solve").start();

                    size_type startBlockIndexInextensibilityLM = this->startSubBlockSpaceIndex("inextensibility-lm");
                    auto lambda = this->functionSpaceInextensibilityLM()->element(XVec,rowStartInVector+startBlockIndexInextensibilityLM);

                    if( BuildNonCstPart )
                    {
                        auto inextensibleLevelsetsExpr = Feel::FeelModels::globalLevelsetExpr( M_inextensibleLevelsets );
                        auto inextensibleLevelsetsDeltaExpr = Feel::FeelModels::levelsetDelta( inextensibleLevelsetsExpr, cst(M_globalLevelsetThicknessInterface) );
                        auto inextensibleLevelsetsDelta = vf::project(
                                _space=this->M_levelsetSpaceManager->functionSpaceScalar(), 
                                _range=this->M_levelsetSpaceManager->rangeMeshElements(),
                                _expr=inextensibleLevelsetsDeltaExpr
                                );

                        form1( _test=this->functionSpace(), _vector=R,
                               _rowstart=rowStartInVector ) +=
                            integrate( _range=this->M_rangeInextensibilityLM,
                                       _expr=idv(lambda)*trace((Id-NxN)*grad(v))*idv(inextensibleLevelsetsDelta),
                                       _geomap=this->geomap()
                                       );
                        form1( _test=this->functionSpaceInextensibilityLM(), _vector=R,
                               _rowstart=rowStartInVector+startBlockIndexInextensibilityLM ) += 
                            integrate( _range=this->M_rangeInextensibilityLM,
                                       _expr=id(lambda)*trace((Id-NxN)*gradv(u))*idv(inextensibleLevelsetsDelta),
                                       _geomap=this->geomap()
                                       );
                    }

                    double timeElapsedInextensibility_LagrangeMult = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateResidualAdditional",
                            "assembly inextensibility (lagrange-multiplier) in "+(boost::format("%1% s") 
                                %timeElapsedInextensibility_LagrangeMult).str() );
                }
            }
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","updateResidualAdditional","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

} // namespace FeelModels
} // namespace Feel
