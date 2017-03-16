/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/multifluid/multifluid.hpp>

#include <feel/feelfilters/loadmesh.hpp>

namespace Feel {
namespace FeelModels {

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
MULTIFLUID_CLASS_TEMPLATE_TYPE::MultiFluid(
        std::string const& prefix,
        WorldComm const& wc,
        std::string const& subPrefix,
        std::string const& rootRepository )
: super_type( prefixvm(prefix,"fluid"), false, wc, subPrefix, self_type::expandStringFromSpec( rootRepository ) )
, M_prefix( prefix )
{
    //-----------------------------------------------------------------------------//
    // Load parameters
    this->loadParametersFromOptionsVm();
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
typename MULTIFLUID_CLASS_TEMPLATE_TYPE::self_ptrtype
MULTIFLUID_CLASS_TEMPLATE_TYPE::New(
        std::string const& prefix,
        WorldComm const& wc,
        std::string const& subPrefix,
        std::string const& rootRepository )
{
    self_ptrtype new_multifluid( new self_type(prefix, wc, subPrefix, rootRepository) );
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
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::build()
{
    CHECK( M_nFluids >= 2 ) << "Multifluid must contain at least 2 fluids.\n";
    uint16_type nLevelSets = M_nFluids - 1;

    this->log("MultiFluid", "build", "start");

    // Read mesh info from multifluid options
    if (Environment::vm().count(prefixvm(this->prefix(),"mshfile").c_str()))
        this->setMshfileStr( Environment::expand( soption(_prefix=this->prefix(),_name="mshfile") ) );
    if (Environment::vm().count(prefixvm(this->prefix(),"geofile").c_str()))
        this->setGeofileStr( Environment::expand( soption(_prefix=this->prefix(),_name="geofile") ) );

    // Build inherited FluidMechanics
    this->loadMesh( this->createMesh() );
    // Deal with lagrange-multiplier inextensibility space
    if( this->M_enableInextensibility && this->inextensibilityMethod() == "lagrange-multiplier" )
    {
        M_spaceInextensibilityLM = space_inextensibilitylm_type::New( 
                _mesh=this->mesh(), 
                _worldscomm=this->localNonCompositeWorldsComm()
                );
    }

    M_globalLevelset.reset(
            new levelset_type( prefixvm(this->prefix(),"levelset"), this->worldComm(), "", this->rootRepositoryWithoutNumProc() )
            );

    M_globalLevelset->build( this->mesh() );
    if( nLevelSets < 2 )
        M_globalLevelset->getExporter()->setDoExport( false );

    // "Deep" copy
    M_fluidDensityViscosityModel.reset( new densityviscosity_model_type( this->fluidPrefix() ) );
    M_fluidDensityViscosityModel->initFromSpace( this->densityViscosityModel()->dynamicViscositySpace() );
    M_fluidDensityViscosityModel->updateFromModelMaterials( this->modelProperties().materials() );

    M_levelsets.resize( nLevelSets );
    M_levelsetDensityViscosityModels.resize( nLevelSets );
    M_levelsetInterfaceForcesModels.resize( nLevelSets );
    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
        auto levelset_prefix = prefixvm(this->prefix(), (boost::format( "levelset%1%" ) %(i+1)).str());
        M_levelsets[i].reset(
                new levelset_type( levelset_prefix, this->worldComm(), "", this->rootRepositoryWithoutNumProc() )
                );
        M_levelsets[i]->build(
                _space=M_globalLevelset->functionSpace(),
                _space_vectorial=M_globalLevelset->functionSpaceVectorial(),
                _space_markers=M_globalLevelset->functionSpaceMarkers(),
                _reinitializer=M_globalLevelset->reinitializer(),
                _projectorL2=M_globalLevelset->projectorL2(),
                _projectorL2_vectorial=M_globalLevelset->projectorL2Vectorial(),
                _smoother=M_globalLevelset->smoother(),
                _smoother_vectorial=M_globalLevelset->smootherVectorial()
                );

        M_levelsetDensityViscosityModels[i].reset(
                new densityviscosity_model_type( levelset_prefix )
                );
        M_levelsetDensityViscosityModels[i]->initFromMesh( this->mesh(), this->useExtendedDofTable() );
        M_levelsetDensityViscosityModels[i]->updateFromModelMaterials( M_levelsets[i]->modelProperties().materials() );

        if( Environment::vm().count( prefixvm(levelset_prefix, "interface-forces-model").c_str() ) )
        {
            std::vector<std::string> interfaceForcesModels = Environment::vm()[prefixvm(levelset_prefix, "interface-forces-model").c_str()].template as<std::vector<std::string>>();
            // Remove (unwanted) duplicates
            std::sort( interfaceForcesModels.begin(), interfaceForcesModels.end() );
            interfaceForcesModels.erase( std::unique(interfaceForcesModels.begin(), interfaceForcesModels.end()), interfaceForcesModels.end() );

            M_levelsetInterfaceForcesModels[i].resize( interfaceForcesModels.size() );

            for( uint16_type n = 0; n < M_levelsetInterfaceForcesModels[i].size(); ++n )
            {
                M_levelsetInterfaceForcesModels[i][n].reset( 
                        interfaceforces_factory_type::instance().createObject( 
                            interfaceForcesModels[n]
                            )
                        );
                M_levelsetInterfaceForcesModels[i][n]->build( levelset_prefix, M_levelsets[i] );
            }

            M_hasInterfaceForcesModel = true;
        }
    }

    M_interfaceForces.reset( new element_levelset_vectorial_type(this->functionSpaceLevelsetVectorial(), "InterfaceForces") ); 

    this->log("MultiFluid", "build", "finish");
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
                _mesh=new mesh_type( this->worldComm() ),
                _filename=mshfile,
                _worldcomm=this->worldComm(),
                //_prefix=this->prefix(),
                _rebuild_partitions=false,
                _savehdf5=0,
                _straighten=1,
                _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES
                );
    }
    else
    {
        if (this->hasMshfileStr())
        {
            std::string mshfileRebuildPartitions = this->rootRepository() + "/" + this->prefix() + ".msh";

            this->log("createMesh","", "load mesh file : " + this->mshfileStr());
            std::string meshFileExt = fs::path( this->mshfileStr() ).extension().string();
            bool rebuildPartition = boption(_prefix=this->prefix(), _name="gmsh.partition");
            if ( rebuildPartition && meshFileExt != ".msh" )
                CHECK( false ) << "Can not rebuild at this time the mesh partitionining with other format than .msh : TODO";

            mesh = loadMesh(
                    _mesh=new mesh_type( this->worldComm() ),
                    _filename=this->mshfileStr(),
                    _worldcomm=this->worldComm(),
                    _prefix=this->prefix(),
                    _rebuild_partitions=rebuildPartition,
                    _rebuild_partitions_filename=mshfileRebuildPartitions,
                    _partitions=this->worldComm().localSize(),
                    _savehdf5=0,
                    _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES
                    );

            if (rebuildPartition) this->setMshfileStr(mshfileRebuildPartitions);
        }
        else if (this->hasGeofileStr())
        {
            std::string mshfile = this->rootRepository() + "/" + this->prefix() + ".msh";
            this->setMshfileStr(mshfile);

            gmsh_ptrtype geodesc = geo( 
                    _filename=this->geofileStr(),
                    _prefix=this->prefix(),
                    _worldcomm=this->worldComm() 
                    );
            // allow to have a geo and msh file with a filename equal to prefix
            geodesc->setPrefix(this->prefix());
            mesh = createGMSHMesh(
                    _mesh=new mesh_type,_desc=geodesc,
                    _prefix=this->prefix(),_worldcomm=this->worldComm(),
                    _partitions=this->worldComm().localSize(),
                    _directory=this->rootRepository() 
                    );
        }
        this->saveMSHfilePath(fmpath);
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
    this->log("MultiFluid", "init", "start");

    // Initialize FluidMechanics
    super_type::init();

    // Initialize LevelSets
    M_globalLevelset->init();
    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
        M_levelsets[i]->init();
    }

    //this->updateTime( this->timeStepBase()->time() );

    this->updateGlobalLevelset();

    this->log("MultiFluid", "init", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    M_nFluids = ioption( _name="nfluids", _prefix=this->prefix() );

    M_usePicardIterations = boption( _name="use-picard-iterations", _prefix=this->prefix() );

    M_enableSurfaceTension = boption( _name="enable-surface-tension", _prefix=this->prefix() );
    M_hasInterfaceForcesModel = false;

    if( M_enableSurfaceTension )
    {
        std::vector<double> sigma = Environment::vm()[prefixvm(this->prefix(),"surface-tension-coeff").c_str()].template as<std::vector<double> >();

        CHECK( sigma.size() >= M_nFluids - 1 ) << sigma.size() << " surface tension coefficients found.\n"
                                               << "You must at least provide the surface tension coefficients between the "
                                               << M_nFluids - 1
                                               << " levelset fluids and the surrounding fluid.\n";

        M_surfaceTensionCoeff = ublas::symmetric_matrix<double, ublas::upper>(M_nFluids, M_nFluids);
        uint16_type k = 0;
        for( uint16_type i = 0; i < M_surfaceTensionCoeff.size1(); ++i )
            for( uint16_type j = i+1; j < M_surfaceTensionCoeff.size2(); ++j )
            {
                if( k < sigma.size() )
                    M_surfaceTensionCoeff(i,j) = sigma[k];
                else
                    M_surfaceTensionCoeff(i,j) = 0;
                ++k;
            }
    }

    uint16_type nLevelSets = M_nFluids - 1;

    M_hasInextensibility.resize(nLevelSets);
    M_enableInextensibility = false;
    M_inextensibilityMethod.resize(nLevelSets); 
    M_inextensibilityGamma.resize(nLevelSets); 
    for( uint16_type n = 0; n < nLevelSets; ++n )
    {
        auto levelset_prefix = prefixvm(this->prefix(), (boost::format( "levelset%1%" ) %(n+1)).str());

        M_hasInextensibility[n] = boption( _name="enable-inextensibility", _prefix=levelset_prefix );
        if( M_hasInextensibility[n] ) M_enableInextensibility = true;

        M_inextensibilityMethod[n] = soption( _name="inextensibility-method", _prefix=levelset_prefix );
        CHECK( M_inextensibilityMethod[n] == "penalty" || M_inextensibilityMethod[n] == "lagrange-multiplier" ) 
            << "invalid inextensiblity-method " << M_inextensibilityMethod[n]
            << ", should be \"penalty\" or \"lagrange-multiplier\"" << std::endl;

        M_inextensibilityGamma[n] = doption( _name="inextensibility-gamma", _prefix=levelset_prefix );
    }

    M_levelsetReinitEvery.resize(nLevelSets);
    M_levelsetReinitSmoothEvery.resize(nLevelSets);
    for( uint16_type n = 0; n < nLevelSets; ++n )
    {
        auto levelset_prefix = prefixvm(this->prefix(), (boost::format( "levelset%1%" ) %(n+1)).str());
        M_levelsetReinitEvery[n] = ioption( _name="reinit-every", _prefix=levelset_prefix );
        M_levelsetReinitSmoothEvery[n] = ioption( _name="reinit-smooth-every", _prefix=levelset_prefix );
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
boost::shared_ptr<std::ostringstream>
MULTIFLUID_CLASS_TEMPLATE_TYPE::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
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
           << "\n       * rho : " << this->M_fluidDensityViscosityModel->cstRho()
           << "\n       * mu  : " << this->M_fluidDensityViscosityModel->cstMu()
           << "\n       * nu  : " << this->M_fluidDensityViscosityModel->cstNu();
    *_ostr << this->M_fluidDensityViscosityModel->getInfo("")->str();
    for( uint16_type i = 0; i < M_levelsetDensityViscosityModels.size(); ++i )
    {
    *_ostr << "\n     -- fluid " << i+1
           << "\n       * rho : " << this->M_levelsetDensityViscosityModels[i]->cstRho()
           << "\n       * mu  : " << this->M_levelsetDensityViscosityModels[i]->cstMu()
           << "\n       * nu  : " << this->M_levelsetDensityViscosityModels[i]->cstNu();
    *_ostr << this->M_levelsetDensityViscosityModels[i]->getInfo("")->str();
    }

    *_ostr << "\n   Level Sets Parameters";
    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
    *_ostr << "\n     -- level set " << i
           << "\n       * reinit every : " << this->M_levelsetReinitEvery[i];
    }

    *_ostr << "\n   Forces Parameters"
           << "\n     -- has surface tension  : " << std::boolalpha << this->M_enableSurfaceTension;
    if( this->M_enableSurfaceTension )
    {
        for( uint16_type i = 0; i < M_surfaceTensionCoeff.size1(); ++i )
            for( uint16_type j = i+1; j < M_surfaceTensionCoeff.size2(); ++j )
    *_ostr << "\n       * surface tension (" << i << "," << j << ") : " << this->M_surfaceTensionCoeff(i,j);
    }
    *_ostr << "\n     -- has interface forces : " << std::boolalpha << this->M_hasInterfaceForcesModel;
    if( this->M_hasInterfaceForcesModel )
    {
        for( uint16_type i = 0; i < M_levelsets.size(); ++i )
        {
    *_ostr << "\n     -- level set " << i;
            for( uint16_type n = 0; n < M_levelsetInterfaceForcesModels[i].size(); ++n )
    *_ostr << "\n       * force model : " << this->M_levelsetInterfaceForcesModels[i][n]->getInfo()->str();
        }
    }

    *_ostr << "\n";
    *_ostr << super_type::getInfo()->str();
    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
    *_ostr << M_levelsets[i]->getInfo()->str();
    }

    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n\n";

    return _ostr;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
int
MULTIFLUID_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    if( this->M_enableInextensibility && this->inextensibilityMethod() == "lagrange-multiplier" )
        return super_type::nBlockMatrixGraph() + 1;
    else
        return super_type::nBlockMatrixGraph();
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
MULTIFLUID_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    BlocksBaseGraphCSR myBlockGraph = super_type::buildBlockMatrixGraph();

    int indexBlock = super_type::nBlockMatrixGraph();
    if( this->M_enableInextensibility && this->inextensibilityMethod() == "lagrange-multiplier" )
    {
        BlocksStencilPattern patCouplingLM(1, super_type::space_fluid_type::nSpaces,size_type(Pattern::ZERO));
        patCouplingLM(0,1) = size_type(Pattern::COUPLED);

        myBlockGraph(indexBlock,0) = stencil(_test=this->functionSpaceInextensibilityLM(), _trial=this->functionSpace(),
                                             _pattern_block=patCouplingLM,
                                             _diag_is_nonzero=false,_close=false)->graph();
        myBlockGraph(0,indexBlock) = stencil(_test=this->functionSpace(), _trial=this->functionSpaceInextensibilityLM(),
                                             _pattern_block=patCouplingLM.transpose(),
                                             _diag_is_nonzero=false,_close=false)->graph();
        ++indexBlock;
    }

    return myBlockGraph;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
size_type
MULTIFLUID_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    auto res = super_type::nLocalDof();
    if( this->M_enableInextensibility && this->inextensibilityMethod() == "lagrange-multiplier" )
    {
        res += this->functionSpaceInextensibilityLM()->nLocalDofWithGhost();
    }
    return res;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
bool
MULTIFLUID_CLASS_TEMPLATE_TYPE::hasInterfaceForces() const
{
    return this->hasSurfaceTension() || M_hasInterfaceForcesModel;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("MultiFluid", "solve", "start");
    this->timerTool("Solve").start();

    double errorVelocityL2 = 0.;
    int picardIter = 0;
    auto u_old = this->fluidModel()->functionSpaceVelocity()->element();
    do
    {
        u_old = this->fieldVelocity();
        // Update density and viscosity
        this->updateFluidDensityViscosity();
        // Update interface forces
        if( this->hasInterfaceForces() )
        {
            this->updateInterfaceForces();
        }
        // Solve fluid equations
        this->solveFluid();
        // Advect levelsets
        this->advectLevelsets();
        // Reinitialize
        for( uint16_type n = 0; n < M_levelsets.size(); ++n )
        {
	    if( M_levelsetReinitEvery[n] > 0 
			    && (M_levelsets[n]->iterSinceReinit()+1) % M_levelsetReinitEvery[n] == 0 )
		    M_levelsets[n]->reinitialize();
	    else if( M_levelsetReinitSmoothEvery[n] > 0 
			    && (M_levelsets[n]->iterSinceReinit()+1) % M_levelsetReinitSmoothEvery[n] == 0 )
		    M_levelsets[n]->reinitialize( true );
        }

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
        picardIter++;
    } while( M_usePicardIterations && errorVelocityL2 > 0.01);

    // Update global levelset
    this->updateGlobalLevelset();

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","solve","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateTime( double time )
{
    // Levelsets
    for( uint16_type i = 0; i < M_levelsets.size(); ++i)
    {
        M_levelsets[i]->updateTime(time);
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
    for( uint16_type i = 0; i < M_levelsets.size(); ++i)
    {
        M_levelsets[i]->updateTimeStep();
    }

    this->log("MultiFluid", "updateTimeStep", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("MultiFluid","exportResults", "start");

    super_type::exportResults(time);

    if( this->nLevelsets() > 1 )
        M_globalLevelset->exportResults(time);

    for( uint16_type i = 0; i < M_levelsets.size(); ++i)
    {
        M_levelsets[i]->exportResults(time);
    }

    this->log("MultiFluid", "exportResults", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
size_type
MULTIFLUID_CLASS_TEMPLATE_TYPE::initStartBlockIndexFieldsInMatrix()
{
    size_type currentStartIndex = super_type::initStartBlockIndexFieldsInMatrix();

    if( this->M_enableInextensibility && this->inextensibilityMethod() == "lagrange-multiplier" )
    {
        // Add inextensibility LM block index
        this->M_startBlockIndexFieldsInMatrix["inextensibility-lm"] = currentStartIndex++;
    }

    return currentStartIndex;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
int
MULTIFLUID_CLASS_TEMPLATE_TYPE::initBlockVector()
{
    int currentBlockIndex = super_type::initBlockVector();
    if( this->M_enableInextensibility && this->inextensibilityMethod() == "lagrange-multiplier" )
    {
        this->M_blockVectorSolution(currentBlockIndex) = this->backend()->newVector( this->functionSpaceInextensibilityLM() );
        ++currentBlockIndex;
    }

    return currentBlockIndex;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateGlobalLevelset()
{
    this->log("MultiFluid", "updateGlobalLevelset", "start");

    auto minPhi = M_globalLevelset->phi();

    *minPhi = *(M_levelsets[0]->phi());
    for( uint16_type i = 1; i < M_levelsets.size(); ++i )
    {
        *minPhi = vf::project( 
                M_globalLevelset->functionSpace(), 
                elements(M_globalLevelset->mesh()),
                vf::min( idv(minPhi), idv(M_levelsets[i]->phi()) )
                );
    }

    M_globalLevelset->updateInterfaceQuantities();

    this->log("MultiFluid", "updateGlobalLevelset", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateFluidDensityViscosity()
{
    this->log("MultiFluid", "updateFluidDensityViscosity", "start");
    this->timerTool("Solve").start();

    auto globalH = M_globalLevelset->H();

    auto rho = vf::project( 
            this->densityViscosityModel()->dynamicViscositySpace(),
            elements(this->mesh()),
            idv(M_fluidDensityViscosityModel->fieldRho())*idv(globalH)
            );

    auto mu = vf::project( 
            this->densityViscosityModel()->dynamicViscositySpace(),
            elements(this->mesh()),
            idv(M_fluidDensityViscosityModel->fieldMu())*idv(globalH)
            );

    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
        rho += vf::project( 
                this->densityViscosityModel()->dynamicViscositySpace(),
                elements(this->mesh()),
                idv(M_levelsetDensityViscosityModels[i]->fieldRho())*(1. - idv(M_levelsets[i]->H()))
                );
        mu += vf::project( 
                this->densityViscosityModel()->dynamicViscositySpace(),
                elements(this->mesh()),
                idv(M_levelsetDensityViscosityModels[i]->fieldMu())*(1. - idv(M_levelsets[i]->H()))
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

    if( this->hasSurfaceTension() )
    {
        this->timerTool("Solve").start();
        for( uint16_type n = 0; n < M_levelsets.size(); ++n )
        {
            *M_interfaceForces += vf::project( 
                    this->functionSpaceLevelsetVectorial(),
                    elements(this->mesh()),
                    - M_surfaceTensionCoeff(0,n+1)*idv(M_levelsets[n]->K())*idv(M_levelsets[n]->N())*idv(M_levelsets[n]->D())
                    );
        }
        double timeElapsedSurfaceTension = this->timerTool("Solve").stop();
        this->log("MultiFluid", "updateInterfaceForces", "update surface tension forces in "+(boost::format("%1% s")%timeElapsedSurfaceTension).str() );
    }

    if( M_hasInterfaceForcesModel )
    {
        this->timerTool("Solve").start();
        for( uint16_type i = 0; i < M_levelsets.size(); ++i )
        {
            for( uint16_type n = 0; n < M_levelsetInterfaceForcesModels[i].size(); ++n )
            {
                if( M_levelsetInterfaceForcesModels[i][n] )
                {
                    M_levelsetInterfaceForcesModels[i][n]->updateInterfaceForces( M_interfaceForces, false );
                }
            }
        }

        double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
        this->log("MultiFluid", "updateInterfaceForces", "update interface (model) forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
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

    auto u = this->fieldVelocity();
    
    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
        M_levelsets[i]->advect( idv(u) );
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log( "MultiFluid", "advectLevelsets", 
            "level-sets advection done in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEAdditional( 
        sparse_matrix_ptrtype & A, vector_ptrtype & F, bool _BuildCstPart ) const
{
    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("MultiFluid","updateLinearPDEAdditional", "start"+sc );
    this->timerTool("Solve").start();

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

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

        for( uint16_type n = 0; n < M_levelsets.size(); ++n )
        {
            if( this->hasInextensibility(n) )
            {
                auto N = this->M_levelsets[n]->N();
                auto NxN = idv(N)*trans(idv(N));
                auto D = this->M_levelsets[n]->D();

                if( this->inextensibilityMethod(n) == "penalty" )
                {
                    Feel::cout << "Assembling inextensibility (penalty)\n";
                    if( BuildNonCstPart )
                    {
                        bilinearForm_PatternDefault += integrate(
                                _range=elements(mesh),
                                _expr=this->M_inextensibilityGamma[n]*trace((Id-NxN)*gradt(u))*trace((Id-NxN)*grad(v))*idv(D)/h(),
                                _geomap=this->geomap()
                                );
                    }
                }
                if( this->inextensibilityMethod(n) == "lagrange-method" )
                {
                    Feel::cout << "Assembling inextensibility (lagrange-multiplier)\n";
                    size_type startBlockIndexInextensibilityLM = this->startBlockIndexFieldsInMatrix().find("inextensibility-lm")->second;
                    auto lambda = this->functionSpaceInextensibilityLM()->element();
                    if( BuildNonCstPart )
                    {
                        form2( _trial=this->functionSpaceInextensibilityLM(), _test=this->functionSpace(), _matrix=A,
                               _rowstart=rowStartInMatrix,
                               _colstart=colStartInMatrix+startBlockIndexInextensibilityLM ) +=
                            integrate( _range=elements(this->mesh()),
                                       _expr=idt(lambda)*trace((Id-NxN)*grad(v))*idv(D),
                                       _geomap=this->geomap()
                                       );
                        form2( _trial=this->functionSpace(), _test=this->functionSpaceInextensibilityLM(), _matrix=A,
                               _rowstart=rowStartInMatrix+startBlockIndexInextensibilityLM,
                               _colstart=colStartInMatrix ) +=
                            integrate( _range=elements(this->mesh()),
                                       _expr=id(lambda)*trace((Id-NxN)*gradt(u))*idv(D),
                                       _geomap=this->geomap()
                                       );
                    }
                }
            }
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","updateLinearPDEAdditional","assembly additional terms in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobianAdditional( sparse_matrix_ptrtype & J, bool BuildCstPart ) const
{
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidualAdditional( vector_ptrtype & R, bool BuildCstPart ) const
{
}

} // namespace FeelModels
} // namespace Feel
