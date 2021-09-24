/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */
#include <numerics>
#include <fmt/core.h>

#include <feel/feelmodels/multifluid/multifluid.hpp>

#include <feel/feelmodels/levelset/levelsetheavisideexpr.hpp>
#include <feel/feelmodels/levelset/levelsetdeltaexpr.hpp>

#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

#include <boost/iterator/transform_iterator.hpp>

namespace Feel {
namespace FeelModels {

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
MULTIFLUID_CLASS_TEMPLATE_TYPE::MultiFluid(
        std::string const& prefix,
        std::string const& keyword,
        worldcomm_ptr_t const& wc,
        std::string const& subPrefix,
        ModelBaseRepository const& modelRep )
: 
    super_type( prefix, keyword, wc, subPrefix, modelRep ),
    ModelBase( prefix, keyword, wc, subPrefix, modelRep ),
    M_useLagrangeP1iso( false ),
    M_usePicardIterations( false ),
    M_hasInterfaceForcesModel( false ),
    M_enableInextensibility( false ),
    M_hasInextensibilityLM( false ),
    M_doUpdateInextensibilityLM( false ),
    M_globalLevelset( [this]( element_levelset_scalar_ptrtype & globalLs ) { this->updateGlobalLevelset( globalLs ); } )
{
    this->log("MultiFluid","constructor", "start" );

    // timers
    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MultiFluidConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MultiFluidSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MultiFluidPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MultiFluidTimeStepping.data";
    this->addTimerTool("Constructor", nameFileConstructor);
    this->addTimerTool("Solve", nameFileSolve);
    this->addTimerTool("PostProcessing", nameFilePostProcessing);
    this->addTimerTool("TimeStepping", nameFileTimeStepping);

    // option in cfg files
    this->loadParameterFromOptionsVm();

    this->log("MultiFluid","constructor", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("MultiFluid","initMesh", "start");
    this->timerTool("Constructor").start();

    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( M_mesh ) << "mesh generation fail";

    double tElapsed = this->timerTool("Constructor").stop("initMesh");
    this->log("MultiFluid","initMesh", (boost::format("finish in %1% s") %tElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    CHECK( M_nLevelsets >= 1 ) << "Multifluid must contain at least 1 levelset.\n";
    this->log("MultiFluid", "init", "start");
    this->timerTool("Constructor").start();

    M_fluidModel = std::make_shared<fluid_model_type>(
            prefixvm(this->prefix(),"fluid"), "fluid", this->worldCommPtr(),
            this->subPrefix(), this->repository() );
    for( index_type i = 0; i < M_nLevelsets; ++i )
    {
        std::string const lsKeyword = fmt::format( "levelset{}", i );
        std::string const lsPrefix = prefixvm( this->prefix(), lsKeyword );
        M_levelsetModels.push_back( std::make_shared<levelset_model_type>(
                lsPrefix, lsKeyword, this->worldCommPtr(),
                this->subPrefix(), this->repository() )
                );
    }

    //if ( this->physics().empty() )
    //{
        //typename ModelPhysics<mesh_type::nDim>::subphysic_description_type subPhyicsDesc;
        //subPhyicsDesc[M_fluidModel->physicType()] = std::make_tuple( M_fluidModel->keyword(), M_fluidModel );
        ////subPhyicsDesc[M_levelsetModel[0]->physicType()] = std::make_tuple( M_levelsetModel[0]->keyword(), M_levelsetModel[0] );
        //this->initPhysics( this->keyword(), this->modelProperties().models(), subPhyicsDesc );
    //}

    // Physical properties
    if ( !M_materialsProperties )
    {
        auto paramValues = this->modelProperties().parameters().toParameterValues();
        this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    // Init mesh
    if ( !this->mesh() )
        this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    // Init fluid toolbox
    //M_fluidModel->setPhysics( this->physics( M_fluidModel->physicType() ), M_fluidModel->keyword() );
    M_fluidModel->setManageParameterValues( false );
    if ( !M_fluidModel->modelPropertiesPtr() )
    {
        M_fluidModel->setModelProperties( this->modelPropertiesPtr() );
        M_fluidModel->setManageParameterValuesOfModelProperties( false );
    }
    M_fluidModel->setMesh( this->mesh() );
    M_fluidModel->setMaterialsProperties( M_materialsProperties );
    M_fluidModel->init( false );

    // Init levelset toolboxes
    this->initLevelsets();

    // Update current time
    if ( !this->isStationary() )
    {
        this->updateTime( this->timeStepBase()->time() );
        this->setTimeInitial( this->timeStepBase()->timeInitial() );
    }

    // Set advection velocity in levelset toolboxes
    std::string fluidVelocitySymbolUsed = "U"; // TODO : get this info from fluid toolbox
    std::string velConvExprStr = fmt::format( 
            nDim==2 ? 
            "{{ {0}_{1}_0, {0}_{1}_1 }}:{0}_{1}_0:{0}_{1}_1" : 
            "{{ {0}_{1}_0, {0}_{1}_1, {0}_{1}_2 }}:{0}_{1}_0:{0}_{1}_1:{0}_{1}_2",
            M_fluidModel->keyword(), fluidVelocitySymbolUsed 
            );
    auto velConvExpr = expr<nDim,1>( velConvExprStr, "",this->worldComm(),this->repository().expr() );
    for( levelset_model_ptrtype const& lsModel: M_levelsetModels )
        lsModel->setVelocityConvectionExpr( velConvExpr );

    // Update constant parameters
    this->updateParameterValues();

    // Update initial conditions
    this->updateInitialConditions( this->symbolsExpr() );

    // Backend
    this->initAlgebraicBackend();

    // Block vector solution
    this->buildBlockVectorSolution();

    // Update inextensibility LM if needed
    if( this->M_hasInextensibilityLM )
        this->updateInextensibilityLM();
    // Build algebraic data
    if( this->useImplicitCoupling() || this->hasInextensibilityLM() )
    {
        if( buildModelAlgebraicFactory )
            this->initAlgebraicFactory();
    }
    else
    {
        M_fluidModel->initAlgebraicFactory();
        for( levelset_model_ptrtype const& lsModel: M_levelsetModels )
            lsModel->initAlgebraicFactory();
        // Add specific multifluid terms into the matrix assembly
        M_fluidModel->algebraicFactory()->setFunctionLinearAssembly(
                [this]( DataUpdateLinear & data ) { updateLinear_Fluid( data ); }
                );
        M_fluidModel->algebraicFactory()->setFunctionJacobianAssembly(
                [this]( DataUpdateLinear & data ) { updateJacobian_Fluid( data ); }
                );
        M_fluidModel->algebraicFactory()->setFunctionResidualAssembly(
                [this]( DataUpdateLinear & data ) { updateResidual_Fluid( data ); }
                );
        for( size_type i = 0; i < M_levelsetModels.size(); ++i )
        {
            M_levelsetModels[i]->algebraicFactory()->setFunctionLinearAssembly(
                    [this, i]( DataUpdateLinear & data ) { updateLinear_Levelset( i, data ); }
                    );
            M_levelsetModels[i]->algebraicFactory()->setFunctionJacobianAssembly(
                    [this, i]( DataUpdateLinear & data ) { updateJacobian_Levelset( i, data ); }
                    );
            M_levelsetModels[i]->algebraicFactory()->setFunctionResidualAssembly(
                    [this, i]( DataUpdateLinear & data ) { updateResidual_Levelset( i, data ); }
                    );
        }
    }

    //// "Deep" copy FluidMechanics materialProperties
    //M_fluidMaterialProperties.reset( new materials_properties_type( this->fluidModel()->prefix() ) );
    //// Create M_interfaceForces
    //M_interfaceForces.reset( new element_levelset_vectorial_type(this->functionSpaceLevelsetVectorial(), "InterfaceForces") ); 

    //// Init FluidMechanics materialProperties
    //M_fluidMaterialProperties->updateForUse( this->fluidModel()->materialProperties()->dynamicViscositySpace(), this->fluidModel()->modelProperties().materials() );
    //// Init levelsets materialProperties and interfaceForcesModels
    //for( auto const& lsMaterialProperties: M_levelsetsMaterialProperties )
    //{
        //lsMaterialProperties.second->updateForUse(this->fluidModel()->materialProperties()->dynamicViscositySpace(), M_levelsetModels[lsMaterialProperties.first]->modelProperties().materials() );
    //}
    // Update density-viscosity
    //this->updateFluidDensityViscosity();

    // Init post-process
    this->initPostProcess();

    // Mark as updated for use
    this->setIsUpdatedForUse( true );

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("MultiFluid","init", fmt::format("finish in {} s", tElapsedInit ) );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(), this->backend() );
    this->setAlgebraicFactory( algebraicFactory );
    if ( M_fluidModel->hasOperatorPCD() )
        algebraicFactory->preconditionerTool()->attachOperatorPCD( "pcd", M_fluidModel->operatorPCD() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    M_nLevelsets = ioption( _name="nlevelsets", _prefix=this->prefix() );

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
        M_levelsetRedistEvery[lsName] = ioption( _name="redist-every", _prefix=levelset_prefix );
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
           << "\n   Number of levelsets : " << M_nLevelsets;

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
    for( auto const& lsRedistEvery: M_levelsetRedistEvery )
    {
    *_ostr << "\n     -- level set " << "\"" << lsRedistEvery.first << "\""
           << "\n       * redist every : " << lsRedistEvery.second;
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
    *_ostr << this->fluidModel()->getInfo()->str();
    for( auto const& levelset: M_levelsetModels )
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
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateGlobalLevelset( element_levelset_scalar_ptrtype & globalLevelset ) const
{
    this->log("MultiFluid", "updateGlobalLevelset", "start");
    this->timerTool("UpdateLevelsetData").start();

    if( !globalLevelset )
        globalLevelset.reset( new element_levelset_scalar_type( this->functionSpaceLevelset(), "GlobalLevelset" ) );

    globalLevelset->on(
            _range=elements( globalLevelset->mesh() ),
            _expr=this->globalLevelsetExpr()
            );

    double timeElapsed = this->timerTool("UpdateLevelsetData").stop();
    this->log("MultiFluid", "updateGlobalLevelset", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
int
MULTIFLUID_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlocks = this->fluidModel()->nBlockMatrixGraph();
    if( this->useImplicitCoupling() )
    {
        //for( auto const& ls: this->levelsetModels() )
            //nBlocks += ls.second->nBlockMatrixGraph();
        CHECK( false ) << "TODO: implicit coupling\n";
    }
    if( this->hasInextensibilityLM() )
        nBlocks += 1;

    return nBlocks;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
MULTIFLUID_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("MultiFluid","buildBlockMatrixGraph", "start" );

    int nBlocks = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myBlockGraph(nBlocks, nBlocks);

    int nBlocksFluid = this->fluidModel()->nBlockMatrixGraph();

    int startIndexBlockFluid = 0;
    int indexBlock = startIndexBlockFluid;

    auto blockMatFluid = this->fluidModel()->buildBlockMatrixGraph();
    for (int tk1=0; tk1<nBlocksFluid; ++tk1 )
        for (int tk2=0; tk2<nBlocksFluid; ++tk2 )
            myBlockGraph(startIndexBlockFluid+tk1,startIndexBlockFluid+tk2) = blockMatFluid(tk1,tk2);
    indexBlock += nBlocksFluid;

    if( this->useImplicitCoupling() )
    {
        //int nBlocksLevelsets = std::accumulate( 
                //this->levelsetModels().begin(), this->levelsetModels().end(), 0,
                //[]( auto res, auto const& rhs ) {
                //return res + rhs.second->nBlockMatrixGraph();
                //}
            //);
        CHECK( false ) << "TODO: implicit coupling\n";
    }

    if( this->hasInextensibilityLM() )
    {
        this->log("MultiFluid","buildBlockMatrixGraph", "start build inextensibility lagrange-multiplier" );
        int startIndexBlockInextensibility = indexBlock;
        // Matrix stencil
        //BlocksStencilPattern patCouplingLM(1, fluid_model_type::space_fluid_type::nSpaces,size_type(Pattern::ZERO));
        //patCouplingLM(0,0) = size_type(Pattern::COUPLED);

        myBlockGraph(startIndexBlockInextensibility,startIndexBlockFluid) = stencil(
                _test=this->functionSpaceInextensibilityLM(), _trial=this->functionSpaceVelocity(),
                //_pattern_block=patCouplingLM,
                _diag_is_nonzero=false,_close=false)->graph();
        myBlockGraph(startIndexBlockFluid,startIndexBlockInextensibility) = stencil(
                _test=this->functionSpaceVelocity(), _trial=this->functionSpaceInextensibilityLM(),
                //_pattern_block=patCouplingLM.transpose(),
                _diag_is_nonzero=false,_close=false)->graph();
        myBlockGraph(startIndexBlockInextensibility,startIndexBlockInextensibility) = stencil(
                _test=this->functionSpaceInextensibilityLM(), _trial=this->functionSpaceInextensibilityLM(),
                _pattern=size_type(Pattern::ZERO),
                _diag_is_nonzero=false,_close=false)->graph();
        ++indexBlock;
    }

    myBlockGraph.close();

    this->log("MultiFluid","buildBlockMatrixGraph", "finish" );

    return myBlockGraph;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
typename MULTIFLUID_CLASS_TEMPLATE_TYPE::size_type
MULTIFLUID_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    auto res = this->fluidModel()->nLocalDof();
    if( this->useImplicitCoupling() )
    {
        //res += std::accumulate( 
            //this->levelsetModels().begin(), this->levelsetModels().end(), 0,
            //[]( auto res, auto const& rhs ) {
                //return res + rhs.second->nLocalDof();
            //}
            //);
            CHECK( false ) << "TODO: implicit coupling\n";

    }
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
    for( auto const& ls: M_levelsetModels )
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
    auto inextensibleLevelsets = vf::project(
            _space=this->M_levelsetSpaceManager->functionSpaceScalar(), 
            _range=this->M_levelsetSpaceManager->rangeMeshElements(),
            _expr=inextensibleLevelsetsExpr
            );
    auto dirac = vf::project( 
            _space=this->fluidModel()->functionSpacePressure(), _range=this->fluidModel()->rangeMeshElements(),
            _expr=Feel::FeelModels::levelsetDelta( inextensibleLevelsets, M_globalLevelsetThicknessInterface )
            );
    auto it_elt = this->mesh()->beginOrderedElement();
    auto en_elt = this->mesh()->endOrderedElement();

    const rank_type pid = this->mesh()->worldCommElements().localRank();
    const int ndofv = fluid_model_type::space_pressure_type::fe_type::nDof;
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

    // Solve
    if( this->useImplicitCoupling() )
    {
        this->solveImplicitCoupling();
    }
    else /* explicit coupling */
    {
        if( M_usePicardIterations )
            this->solvePicard();
        else
            this->solveSemiImplicitCoupling();
    }

    // Redistantiate
    for( auto const& ls: M_levelsetModels )
    {
        if( M_levelsetRedistEvery.at(ls.first) > 0 
                && (ls.second->iterSinceRedistanciation()+1) % M_levelsetRedistEvery.at(ls.first) == 0 )
            ls.second->redistanciate();
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","solve","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::solveSemiImplicitCoupling()
{
    if ( this->hasInextensibilityLM() ) /* inextensibility-lm block in fluid */
    {
        this->updateParameterValues();

        // Update inextensibility LM if requested
        if( M_doUpdateInextensibilityLM )
            this->updateInextensibilityLM();
        // Rebuild matrix and vector
        if( this->M_doRebuildMatrixVector )
        {
            // Rebuild algebraic matrix and vector
            auto graph = this->buildMatrixGraph();
            M_algebraicFactory->rebuildMatrixVector( graph, graph->mapRow().indexSplit() );
            // Rebuild solution vector
            this->buildBlockVectorSolution();
            // Rebuild backend (required since PETSc stores the size of the matrix, 
            // which can change here because of the LM space support)
            this->backend()->clear();

            M_doRebuildMatrixVector = false;
        }

        this->fluidModel()->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("fluid") );

        M_blockVectorSolution.updateVectorFromSubVectors();

        // Solve fluid
        M_algebraicFactory->solve( this->fluidModel()->solverName(), M_blockVectorSolution.vectorMonolithic() );
        M_blockVectorSolution.localize();
    }
    else
    {
        // Solve fluid equations (with direct assembly of interface forces)
        M_fluidModel->solve();
    }
    // Advect levelsets
    this->advectLevelsets();
    // Update density and viscosity
    this->updateFluidDensityViscosity();
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::solveImplicitCoupling()
{
    this->updateParameterValues();

    this->fluidModel()->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("fluid") );
    //for( auto const& ls: M_levelsetModels )
    //{
        //ls.second->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex( ls.first ) );
    //}
    CHECK(false) << "TODO: implicit coupling\n";

    M_blockVectorSolution.updateVectorFromSubVectors();
    //TODO
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::solvePicard()
{
    this->solveSemiImplicitCoupling();

    double errorVelocityL2 = 0.;
    int picardIter = 0;
    auto u_old = this->fluidModel()->functionSpaceVelocity()->element();
    do
    {
        picardIter++;
        u_old = this->fieldVelocity();
        // Sub-solves
        if( M_doUpdateInextensibilityLM )
            this->updateInextensibilityLM();
        this->solveSemiImplicitCoupling();
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

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    bool _BuildCstPart = data.buildCstPart();
    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("MultiFluid","updateLinearPDE", "start"+sc );
    this->timerTool("Solve").start();

    //sparse_matrix_ptrtype& A = data.matrix();
    //vector_ptrtype& F = data.rhs();
    //bool BuildNonCstPart = !_BuildCstPart;
    //bool BuildCstPart = _BuildCstPart;

    // Update fluid
    this->fluidModel()->updateLinearPDE( data );

    // Update interface forces
    this->updateLinearPDEInterfaceForces( data );

    // Update inextensibility
    this->updateLinearPDEInextensibility( data );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","updateLinearPDE","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEInterfaceForces( DataUpdateLinear & data ) const
{
    // Update interface forces
    if( this->hasInterfaceForces() )
    {
        this->log("MultiFluid", "updateLinearPDEInterfaceForces", "start: update interface forces");
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
            this->log("MultiFluid", "updateLinearPDEInterfaceForces", "update interface (model) forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        if( M_additionalInterfaceForcesModel.size() > 0 )
        {
            this->timerTool("Solve").start();
            for( auto const& f: M_additionalInterfaceForcesModel )
                f.second->updateFluidInterfaceForcesLinearPDE( data );

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateLinearPDEInterfaceForces", "update additional interface forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        double timeElapsed = this->timerTool("Solve").stop();
        this->log( "MultiFluid", "updateLinearPDEInterfaceForces", 
                "interface forces updated in "+(boost::format("%1% s") %timeElapsed).str() );
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEInextensibility( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;

    if( this->M_enableInextensibility )
    {
        if( BuildNonCstPart )
        {
            auto mesh = this->mesh();
            auto XhV = this->functionSpaceVelocity();

            auto const& u = this->fieldVelocity();
            auto const& v = u;
            auto Id = vf::Id<nDim, nDim>();

            auto myBfV = form2( 
                    _test=XhV, _trial=XhV, _matrix=A,
                    _rowstart=this->fluidModel()->rowStartInMatrix(),
                    _colstart=this->fluidModel()->colStartInMatrix()
                    );

            for( auto const& ls: M_levelsetModels )
            {
                std::string const& lsName = ls.first;
                if( this->hasInextensibility(lsName) && this->inextensibilityMethod(lsName) == "penalty" )
                {
                    auto N = ls.second->N();
                    auto NxN = idv(N)*trans(idv(N));
                    auto D = ls.second->D();

                    this->timerTool("Solve").start();

                    myBfV += integrate(
                            _range=elements(mesh),
                            _expr=this->M_inextensibilityGamma.at(lsName)*trace((Id-NxN)*gradt(u))*trace((Id-NxN)*grad(v))*idv(D)/h(),
                            _geomap=this->geomap()
                            );

                    double timeElapsedInextensibility_Penalty = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateLinearPDEInextensibility",
                            "assembly inextensibility (penalty) in "+(boost::format("%1% s") 
                                %timeElapsedInextensibility_Penalty).str() );
                }
            }

            if( this->hasInextensibilityLM() )
            {
                CHECK( this->hasStartSubBlockSpaceIndex("inextensibility-lm") ) << " start dof index for inextensibility-lm is not present\n";
                this->timerTool("Solve").start();

                size_type startBlockIndexInextensibilityLM = this->startSubBlockSpaceIndex("inextensibility-lm");
                auto lambda = this->functionSpaceInextensibilityLM()->element();

                auto inextensibleLevelsetsExpr = Feel::FeelModels::globalLevelsetExpr( M_inextensibleLevelsets );
                auto inextensibleLevelsets = vf::project(
                        _space=this->M_levelsetSpaceManager->functionSpaceScalar(), 
                        _range=this->M_levelsetSpaceManager->rangeMeshElements(),
                        _expr=inextensibleLevelsetsExpr
                        );
                auto inextensibleLevelsetsDeltaExpr = Feel::FeelModels::levelsetDelta( inextensibleLevelsets, M_globalLevelsetThicknessInterface );
                auto N = trans(gradv(inextensibleLevelsets)) / sqrt( gradv(inextensibleLevelsets)*trans(gradv(inextensibleLevelsets)) );
                auto NxN = N*trans(N);

                form2( _trial=this->functionSpaceInextensibilityLM(), _test=XhV, 
                        _matrix=A,
                        _rowstart=this->fluidModel()->rowStartInMatrix(),
                        _colstart=startBlockIndexInextensibilityLM ) +=
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=idt(lambda)*trace((Id-NxN)*grad(v))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );
                form2( _trial=XhV, _test=this->functionSpaceInextensibilityLM(), 
                        _matrix=A,
                        _rowstart=startBlockIndexInextensibilityLM,
                        _colstart=this->fluidModel()->colStartInMatrix() ) +=
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=id(lambda)*trace((Id-NxN)*gradt(u))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );

                double timeElapsedInextensibility_LagrangeMult = this->timerTool("Solve").stop();
                this->log("MultiFluid","updateLinearPDEInextensibility",
                        "assembly inextensibility (lagrange-multiplier) in "+(boost::format("%1% s") 
                            %timeElapsedInextensibility_LagrangeMult).str() );
            }
        }
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    this->fluidModel()->updateLinearPDEDofElimination( data );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    bool _BuildCstPart = data.buildCstPart();
    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("MultiFluid","updateJacobian", "start"+sc );
    this->timerTool("Solve").start();

    //const vector_ptrtype& XVec = data.currentSolution();
    //sparse_matrix_ptrtype& J = data.jacobian();

    //bool BuildNonCstPart = !_BuildCstPart;
    //bool BuildCstPart = _BuildCstPart;

    // Update interface forces
    this->updateJacobianInterfaceForces( data );
    // Update inextensibility
    this->updateJacobianInextensibility( data );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","updateJacobian","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobianInterfaceForces( DataUpdateJacobian & data ) const
{
    if( this->hasInterfaceForces() )
    {
        this->log("MultiFluid", "updateJacobianInterfaceForces", "start: update interface forces");
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
            this->log("MultiFluid", "updateJacobianInterfaceForces", "update interface (model) forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        if( M_additionalInterfaceForcesModel.size() > 0 )
        {
            this->timerTool("Solve").start();
            for( auto const& f: M_additionalInterfaceForcesModel )
                f.second->updateFluidInterfaceForcesJacobian( data );

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateJacobianInterfaceForces", "update additional interface forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        double timeElapsed = this->timerTool("Solve").stop();
        this->log( "MultiFluid", "updateJacobianInterfaceForces", 
                "interface forces updated in "+(boost::format("%1% s") %timeElapsed).str() );
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobianInextensibility( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();

    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;

    if( this->M_enableInextensibility )
    {
        if( BuildNonCstPart )
        {
            auto mesh = this->mesh();
            auto XhV = this->functionSpaceVelocity();

            auto myBfV = form2(
                    _test=XhV,_trial=XhV,_matrix=J,
                    //_pattern=size_type(Pattern::COUPLED),
                    _rowstart=this->fluidModel()->rowStartInMatrix(),
                    _colstart=this->fluidModel()->colStartInMatrix()
                    );

            auto u = XhV->element(XVec, this->fluidModel()->rowStartInVector());
            auto const& v = this->fieldVelocity();
            auto Id = vf::Id<nDim, nDim>();

            for( auto const& ls: M_levelsetModels )
            {
                std::string const& lsName = ls.first;
                if( this->hasInextensibility(lsName) && this->inextensibilityMethod(lsName) == "penalty" )
                {
                    auto N = ls.second->N();
                    auto NxN = idv(N)*trans(idv(N));
                    auto D = ls.second->D();

                    this->timerTool("Solve").start();

                    if( BuildNonCstPart )
                    {
                        myBfV += integrate(
                                _range=elements(mesh),
                                _expr=this->M_inextensibilityGamma.at(lsName)*trace((Id-NxN)*gradt(u))*trace((Id-NxN)*grad(v))*idv(D)/h(),
                                _geomap=this->geomap()
                                );
                    }

                    double timeElapsedInextensibility_Penalty = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateJacobianInextensibility",
                            "assembly inextensibility (penalty) in "+(boost::format("%1% s") 
                                %timeElapsedInextensibility_Penalty).str() );
                }
            }

            if( this->hasInextensibilityLM() )
            {
                CHECK( this->hasStartSubBlockSpaceIndex("inextensibility-lm") ) << " start dof index for inextensibility-lm is not present\n";
                this->timerTool("Solve").start();

                size_type startBlockIndexInextensibilityLM = this->startSubBlockSpaceIndex("inextensibility-lm");
                auto lambda = this->functionSpaceInextensibilityLM()->element();

                auto inextensibleLevelsetsExpr = Feel::FeelModels::globalLevelsetExpr( M_inextensibleLevelsets );
                auto inextensibleLevelsets = vf::project(
                        _space=this->M_levelsetSpaceManager->functionSpaceScalar(), 
                        _range=this->M_levelsetSpaceManager->rangeMeshElements(),
                        _expr=inextensibleLevelsetsExpr
                        );
                auto inextensibleLevelsetsDeltaExpr = Feel::FeelModels::levelsetDelta( inextensibleLevelsets, M_globalLevelsetThicknessInterface );
                auto N = trans(gradv(inextensibleLevelsets)) / sqrt( gradv(inextensibleLevelsets)*trans(gradv(inextensibleLevelsets)) );
                auto NxN = N*trans(N);

                form2( _trial=this->functionSpaceInextensibilityLM(), _test=this->functionSpaceVelocity(), 
                        _matrix=J,
                        _rowstart=this->fluidModel()->rowStartInMatrix(),
                        _colstart=startBlockIndexInextensibilityLM ) +=
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=idt(lambda)*trace((Id-NxN)*grad(v))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );
                form2( _trial=this->functionSpaceVelocity(), _test=this->functionSpaceInextensibilityLM(), 
                        _matrix=J,
                        _rowstart=startBlockIndexInextensibilityLM,
                        _colstart=this->fluidModel()->colStartInMatrix() ) +=
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=id(lambda)*trace((Id-NxN)*gradt(u))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );

                double timeElapsedInextensibility_LagrangeMult = this->timerTool("Solve").stop();
                this->log("MultiFluid","updateJacobianInextensibility",
                        "assembly inextensibility (lagrange-multiplier) in "+(boost::format("%1% s") 
                            %timeElapsedInextensibility_LagrangeMult).str() );
            }
        }
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    this->fluidModel()->updateJacobianDofElimination( data );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    bool _BuildCstPart = data.buildCstPart();
    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("MultiFluid","updateResidual", "start"+sc );
    this->timerTool("Solve").start();

    //const vector_ptrtype& XVec = data.currentSolution();
    //vector_ptrtype& R = data.residual();
    //bool BuildCstPart = _BuildCstPart;
    //bool BuildNonCstPart = !BuildCstPart;
    //bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    // Update interface forces
    this->updateResidualInterfaceForces( data );

    // Update inextensibility
    this->updateResidualInextensibility( data );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","updateResidual","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidualInterfaceForces( DataUpdateResidual & data ) const
{
    if( this->hasInterfaceForces() )
    {
        this->log("MultiFluid", "updateResidualInterfaceForces", "start: update interface forces");
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
            this->log("MultiFluid", "updateResidualInterfaceForces", "update interface (model) forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        if( M_additionalInterfaceForcesModel.size() > 0 )
        {
            this->timerTool("Solve").start();
            for( auto const& f: M_additionalInterfaceForcesModel )
                f.second->updateFluidInterfaceForcesResidual( data );

            double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
            this->log("MultiFluid", "updateResidualInterfaceForces", "update additional interface forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
        }

        double timeElapsed = this->timerTool("Solve").stop();
        this->log( "MultiFluid", "updateResidualInterfaceForces", 
                "interface forces updated in "+(boost::format("%1% s") %timeElapsed).str() );
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidualInextensibility( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    if( this->M_enableInextensibility )
    {
        if( BuildNonCstPart )
        {
            auto mesh = this->mesh();
            auto XhV = this->functionSpaceVelocity();

            auto myLV = form1( 
                    _test=XhV,_vector=R,
                    //_pattern=size_type(Pattern::COUPLED),
                    _rowstart=this->fluidModel()->rowStartInVector()
                    );

            auto u = XhV->element(XVec, this->fluidModel()->rowStartInVector());
            auto const& v = this->fieldVelocity();
            auto Id = vf::Id<nDim, nDim>();

            for( auto const& ls: M_levelsetModels )
            {
                std::string const& lsName = ls.first;
                if( this->hasInextensibility(lsName) && this->inextensibilityMethod(lsName) == "penalty" )
                {
                    auto N = ls.second->N();
                    auto NxN = idv(N)*trans(idv(N));
                    auto D = ls.second->D();

                    this->timerTool("Solve").start();

                    if( !UseJacobianLinearTerms )
                    {
                        myLV += integrate(
                                _range=elements(mesh),
                                _expr=this->M_inextensibilityGamma.at(lsName)*trace((Id-NxN)*gradv(u))*trace((Id-NxN)*grad(v))*idv(D)/h(),
                                _geomap=this->geomap()
                                );
                    }

                    double timeElapsedInextensibility_Penalty = this->timerTool("Solve").stop();
                    this->log("MultiFluid","updateResidualInextensibility",
                            "assembly inextensibility (penalty) in "+(boost::format("%1% s") 
                                %timeElapsedInextensibility_Penalty).str() );
                }
            }

            if( this->hasInextensibilityLM() )
            {
                CHECK( this->hasStartSubBlockSpaceIndex("inextensibility-lm") ) << " start dof index for inextensibility-lm is not present\n";
                this->timerTool("Solve").start();

                size_type startBlockIndexInextensibilityLM = this->startSubBlockSpaceIndex("inextensibility-lm");
                auto lambda = this->functionSpaceInextensibilityLM()->element(XVec,startBlockIndexInextensibilityLM);

                auto inextensibleLevelsetsExpr = Feel::FeelModels::globalLevelsetExpr( M_inextensibleLevelsets );
                auto inextensibleLevelsets = vf::project(
                        _space=this->M_levelsetSpaceManager->functionSpaceScalar(), 
                        _range=this->M_levelsetSpaceManager->rangeMeshElements(),
                        _expr=inextensibleLevelsetsExpr
                        );
                auto inextensibleLevelsetsDeltaExpr = Feel::FeelModels::levelsetDelta( inextensibleLevelsets, M_globalLevelsetThicknessInterface );
                auto N = trans(gradv(inextensibleLevelsets)) / sqrt( gradv(inextensibleLevelsets)*trans(gradv(inextensibleLevelsets)) );
                auto NxN = N*trans(N);

                form1( _test=this->functionSpaceVelocity(), _vector=R,
                        _rowstart=this->fluidModel()->rowStartInVector() ) +=
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=idv(lambda)*trace((Id-NxN)*grad(v))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );
                form1( _test=this->functionSpaceInextensibilityLM(), _vector=R,
                        _rowstart=startBlockIndexInextensibilityLM ) += 
                    integrate( _range=this->M_rangeInextensibilityLM,
                            _expr=id(lambda)*trace((Id-NxN)*gradv(u))*inextensibleLevelsetsDeltaExpr,
                            _geomap=this->geomap()
                            );

                double timeElapsedInextensibility_LagrangeMult = this->timerTool("Solve").stop();
                this->log("MultiFluid","updateResidualInextensibility",
                        "assembly inextensibility (lagrange-multiplier) in "+(boost::format("%1% s") 
                            %timeElapsedInextensibility_LagrangeMult).str() );
            }
        }
    }
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    this->fluidModel()->updateResidualDofElimination( data );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    this->fluidModel()->updateParameterValues();
}

//MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
//void
//MULTIFLUID_CLASS_TEMPLATE_TYPE::updateTime( double time )
//{
    //// Levelsets
    //for( auto const& ls: M_levelsetModels)
    //{
        //ls.second->updateTime(time);
    //}
    //// This
    //super_type::updateTime(time);
//}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    this->log("MultiFluid", "startTimeStep", "start");

    this->fluidModel()->startTimeStep();
    this->updateTime( this->timeStepBase()->time() );

    this->log("MultiFluid", "startTimeStep", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->log("MultiFluid", "updateTimeStep", "start");
    // Fluid
    this->fluidModel()->updateTimeStep();
    // Levelsets
    for( auto const& ls: M_levelsetModels)
    {
        ls.second->updateTimeStep();
    }
    // Self
    this->updateTime( this->timeStepBase()->time() );

    this->log("MultiFluid", "updateTimeStep", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("MultiFluid","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    auto const symbolsExpr = this->symbolsExpr();
    M_fluidModel->exportResults( time, symbolsExpr );
    for( auto const& ls: M_levelsetModels)
        ls.second->exportResults( time, symbolsExpr );

    std::map<std::string, typename interfaceforces_model_type::element_ptrtype> interfaceForcesFields;
    for( auto const& lsInterfaceForces: M_levelsetInterfaceForcesModels )
    {
        std::string const lsName = lsInterfaceForces.first;
        for( auto const& force: lsInterfaceForces.second )
        {
            interfaceForcesFields[prefixvm(lsName, force.first)] = force.second->lastInterfaceForce();
        }
    }
    for( auto const& force: M_additionalInterfaceForcesModel )
    {
        interfaceForcesFields[prefixvm("additional", force.first)] = force.second->lastInterfaceForce();
    }
#if 0
    //auto const ls_fields = []( std::pair<std::string, levelset_model_ptrtype> const& p ) 
        //-> decltype( p.second->allFields() )
    //{
        //return p.second->allFields( p.second->keyword() );
    //};
    struct ls_fields {
        typedef decltype( std::declval<levelset_model_type>().allFields() ) result_type;
        result_type operator()( std::pair<std::string, levelset_model_ptrtype> const& p ) const {
            return p.second->allFields( p.second->keyword() );
        }
    };
    auto levelsetsFields = Feel::zip_with( 
            []( auto it1, auto it2 ) {
                if constexpr( Feel::is_iterable_v<decltype(*it1)> )
                {
                    using key_type = decltype( it1->begin()->first );
                    using value_type = decltype( it1->begin()->second );
                    std::map<key_type, value_type> m;
                    for( auto it = it1; it != it2; ++it )
                    {
                        auto const& itmap = *it;
                        m.insert( itmap.begin(), itmap.end() );
                    }
                    return m;
                }
                else
                {
                    using key_type = decltype( it1->first );
                    using value_type = decltype( it1->second );
                    return std::map<key_type, value_type>( it1, it2 );
                }
            },
            boost::iterators::transform_iterator( M_levelsetModels.begin(), ls_fields{} ),
            boost::iterators::transform_iterator( M_levelsetModels.end(), ls_fields{} )
            );

    auto fields = hana::flatten( hana::make_tuple( 
            M_fluidModel->allFields( M_fluidModel->keyword() ), 
            //M_levelsetModels.begin()->second->allFields( M_levelsetModels.begin()->second->keyword() ), 
            levelsetsFields,
            hana::make_tuple( 
                std::make_pair( prefixvm(this->prefix(),"global-levelset.phi"), this->globalLevelsetElt() ),
                interfaceForcesFields
                ) 
            ) );
    //TODO: support fields for ALL levelsets
    auto exprExport = M_fluidModel->exprPostProcessExports( symbolsExpr, M_fluidModel->keyword() );
    // TODO: add exprPostProcessExports support in LevelSet
    this->executePostProcessExports( M_exporter, time, fields, symbolsExpr, exprExport );
#else
    auto mfields = this->modelFields();
    auto symbolExpr = this->symbolsExpr( mfields );
    this->executePostProcessExports( M_exporter, time, mfields, symbolExpr );
#endif
    // Export measures
    this->exportMeasures( time );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
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
    auto const& u = this->fieldVelocity();
    auto const& p = this->fieldPressure();
    std::string matName = this->fluidModel()->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
    auto sigmav = Feel::FeelModels::fluidMecStressTensor(gradv(u),idv(p),*this->fluidModel()->materialProperties(),matName,true,2*fluid_model_type::nOrderVelocity,true);

    auto const& phi = this->levelsetModel(name)->phi();
    auto N_expr = trans(gradv(phi)) / sqrt( gradv(phi) * trans(gradv(phi)) );
    auto const D_expr = this->levelsetModel(name)->diracExpr();
    return integrate(_range=elements(this->mesh()),
                     _expr= - sigmav * N_expr * D_expr,
                     _geomap=this->geomap() ).evaluate();
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::initLevelsets()
{
    this->log("MultiFluid", "initLevelsets", "start");
    // Get levelset mesh
    mesh_ptrtype mesh;
    if( this->M_useLagrangeP1iso )
    {
        // Build Lagrange P1 iso-U mesh and build levelsets with it
        M_opLagrangeP1iso = lagrangeP1( this->fluidModel()->functionSpaceVelocity()->compSpace() );
        mesh = this->M_opLagrangeP1iso->mesh(); 
    }
    else
    {
        // Else build from common mesh
        mesh = this->mesh();
    }
    // Build levelsets space manager
    M_levelsetSpaceManager = std::make_shared<levelset_space_manager_type>( mesh, this->globalLevelsetPrefix() );
    // TODO: Temporary hack: ensures the defaults levelset spaces are built for interfaceForces
    M_levelsetSpaceManager->createFunctionSpaceDefault();
    // Build levelsets tool manager
    M_levelsetToolManager = std::make_shared<levelset_tool_manager_type>( M_levelsetSpaceManager, this->globalLevelsetPrefix() );
    // Update global levelset thickness interface
    if( Environment::vm().count( prefixvm(this->globalLevelsetPrefix(),"thickness-interface").c_str() ) )
        M_globalLevelsetThicknessInterface = doption( _name="thickness-interface", _prefix=this->globalLevelsetPrefix() );
    else
        M_globalLevelsetThicknessInterface = 1.5 * M_levelsetSpaceManager->mesh()->hAverage();

    // Init levelsets
    for( levelset_model_ptrtype const& lsModel: M_levelsetModels )
    {
        //// Set global options if unspecified otherwise
        //if( !Environment::vm().count( prefixvm(levelset_prefix,"thickness-interface").c_str() ) )
            //levelset->setThicknessInterface( this->globalLevelsetThicknessInterface() );

        //lsModel->setPhysics( this->physics( lsModel->physicType() ), lsModel->keyword() );
        lsModel->setManageParameterValues( false );
        if ( !lsModel->modelPropertiesPtr() )
        {
            lsModel->setModelProperties( this->modelPropertiesPtr() );
            lsModel->setManageParameterValuesOfModelProperties( false );
        }
        lsModel->setFunctionSpaceManager( M_levelsetSpaceManager );
        lsModel->setToolManager( M_levelsetToolManager );
        //lsModel->setMesh( this->mesh() );
        //// Build levelsets materialProperties
        //M_levelsetsMaterialProperties[lsName].reset(
                //new material_properties_type( levelset_prefix )
                //);
        lsModel->setMaterialsProperties( M_materialsProperties );
        lsModel->init( false );

        // TODO:Build levelsets interfaceForcesModels
        //if( Environment::vm().count( prefixvm(levelset_prefix, "interface-forces-model").c_str() ) )
        //{
            //std::vector<std::string> interfaceForcesModels = Environment::vm()[prefixvm(levelset_prefix, "interface-forces-model").c_str()].template as<std::vector<std::string>>();
            //// Remove (unwanted) duplicates
            //std::sort( interfaceForcesModels.begin(), interfaceForcesModels.end() );
            //interfaceForcesModels.erase( std::unique(interfaceForcesModels.begin(), interfaceForcesModels.end()), interfaceForcesModels.end() );

            //for( uint16_type n = 0; n < interfaceForcesModels.size(); ++n )
            //{
                //std::string const forceName = interfaceForcesModels[n];
                //M_levelsetInterfaceForcesModels[lsName][forceName] = interfaceforces_factory_type::instance().createObject( 
                        //forceName
                        //);
                //M_levelsetInterfaceForcesModels[lsName][forceName]->build( levelset_prefix, M_levelsetModels[lsName], this->fluidModel() );
            //}

            //M_hasInterfaceForcesModel = true;
        //}
    }

    this->log("MultiFluid", "initLevelsets", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("MultiFluid","initPostProcess", "start");
    this->timerTool("Constructor").start();

    std::set<std::string> ppExportsAllFieldsAvailable;
    for ( auto const& s : M_fluidModel->postProcessExportsAllFieldsAvailable() )
        ppExportsAllFieldsAvailable.insert( prefixvm( M_fluidModel->keyword(), s ) );
    for ( auto const& ls : M_levelsetModels )
        for ( auto const& s : ls.second->postProcessExportsAllFieldsAvailable() )
            ppExportsAllFieldsAvailable.insert( prefixvm( ls.second->keyword(), s ) );
    for( auto const& lsForces : M_levelsetInterfaceForcesModels )
        for( auto const& f : lsForces.second )
            ppExportsAllFieldsAvailable.insert( prefixvm( lsForces.first, f.first ) );
    ppExportsAllFieldsAvailable.insert( "global-levelset.phi" );
    this->setPostProcessExportsAllFieldsAvailable( ppExportsAllFieldsAvailable );
    this->setPostProcessExportsPidName( "pid" );
    super_type::initPostProcess();

    if ( !this->postProcessExportsFields().empty() )
    {
        std::string geoExportType="static";//change_coords_only, change, static
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export",
                               _geo=geoExportType,
                               _path=this->exporterPath() );

        if ( this->doRestart() && this->restartPath().empty() )
        {
            if ( M_exporter->doExport() )
                M_exporter->restart(this->timeInitial());
        }
    }

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
    // start or restart the export of measures
    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasuresIO().restart( "time", this->timeInitial() );
        else
            this->postProcessMeasuresIO().setMeasure( "time", this->timeInitial() ); //just to have time in the first column
    }

    double tElpased = this->timerTool("Constructor").stop("initPostProcess");
    this->log("MultiFluid","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::buildBlockVectorSolution()
{
    this->initBlockVectorSolution();
    this->algebraicBlockVectorSolution()->buildVector( this->backend() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
int
MULTIFLUID_CLASS_TEMPLATE_TYPE::initBlockVectorSolution()
{
    auto const& blockVectorSolutionFluid = *(this->fluidModel()->algebraicBlockVectorSolution());
    int nBlockFluid = blockVectorSolutionFluid.size();
    int nBlockLevelsets = std::reduce( 
            this->levelsetModels().begin(), this->levelsetModels().end(), 
                []( int res, auto const& lsModel ) { return res + lsModel->algebraicBlockVectorSolution().size(); }
            );
    int nBlock = nBlockFluid + nBlockLevelsets;

    auto bvs = this->initAlgebraicBlockVectorSolution( nBlock );

    int indexBlock=0;
    size_type currentStartBlockSpaceIndex = 0;

    this->setStartSubBlockSpaceIndex( "fluid", currentStartBlockSpaceIndex );
    int numberOfBlockSpaceFluid = 0;
    for ( int k=0; k<nBlockFluid; ++k )
    {
        bvs->operator()(indexBlock+k) = blockVectorSolutionFluid(k);
        numberOfBlockSpaceFluid += blockVectorSolutionFluid(k)->map().numberOfDofIdToContainerId();
    }
    indexBlock += nBlockFluid;
    currentStartBlockSpaceIndex += numberOfBlockSpaceFluid;

    for( levelset_model_ptrtype const& lsModel: this->levelsetModels() )
    {
        auto const& blockVectorSolutionLevelset = *(lsModel->algebraicBlockVectorSolution());
        int nBlockLevelset = blockVectorSolutionLevelset.size();
        this->setStartSubBlockSpaceIndex( lsModel->keyword(), currentStartBlockSpaceIndex );
        int numberOfBlockSpaceLevelset = 0;
        for ( int k=0; k<nBlockLevelset; ++k )
        {
            bvs->operator()(indexBlock+k) = blockVectorSolutionLevelset(k);
            numberOfBlockSpaceLevelset += blockVectorSolutionLevelset(k)->map().numberOfDofIdToContainerId();
        }
        indexBlock += nBlockLevelset;
        currentStartBlockSpaceIndex += numberOfBlockSpaceLevelset;
    }


    if( this->hasInextensibilityLM() )
    {
        this->setStartSubBlockSpaceIndex( "inextensibility-lm", currentStartBlockSpaceIndex );
        bvs->operator()(indexBlock) = this->backend()->newVector( this->functionSpaceInextensibilityLM() );
        int numberOfBlockSpaceInextensibilityLM = bvs->operator()(indexBlock)->map().numberOfDofIdToContainerId();
        ++indexBlock;
        currentStartBlockSpaceIndex += numberOfBlockSpaceInextensibilityLM;
    }

    return indexBlock;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
bool
MULTIFLUID_CLASS_TEMPLATE_TYPE::useImplicitCoupling() const
{
    return false;
}

//MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
//void
//MULTIFLUID_CLASS_TEMPLATE_TYPE::updateFluidDensityViscosity()
//{
    //this->log("MultiFluid", "updateFluidDensityViscosity", "start");
    //this->timerTool("Solve").start();

    //auto globalH = Feel::FeelModels::levelsetHeaviside( 
            //this->globalLevelsetExpr(),
            //cst( this->globalLevelsetThicknessInterface() )
            //);

    //auto rho = vf::project( 
        //_space=this->fluidModel()->materialProperties()->dynamicViscositySpace(),
        //_range=elements(this->mesh()),
        //_expr=idv(M_fluidMaterialProperties->fieldRho())*globalH
            //);

    //auto mu = vf::project( 
        //_space=this->fluidModel()->materialProperties()->dynamicViscositySpace(),
        //_range=elements(this->mesh()),
        //_expr=idv(M_fluidMaterialProperties->fieldMu())*globalH
            //);

    //for( auto const& ls: M_levelsetModels )
    //{
        //auto Hi = Feel::FeelModels::levelsetHeaviside( 
                //idv( ls.second->phi() ),
                //cst( this->globalLevelsetThicknessInterface() )
                //);
        //rho += vf::project( 
            //_space=this->fluidModel()->materialProperties()->dynamicViscositySpace(),
            //_range=elements(this->mesh()),
            //_expr=idv(M_levelsetsMaterialProperties[ls.first]->fieldRho())*(1. - Hi)
                //);
        //mu += vf::project( 
            //_space=this->fluidModel()->materialProperties()->dynamicViscositySpace(),
            //_range=elements(this->mesh()),
            //_expr=idv(M_levelsetsMaterialProperties[ls.first]->fieldMu())*(1. - Hi)
                //);
    //}

    //this->fluidModel()->updateRho( idv(rho) );
    //this->fluidModel()->updateMu( idv(mu) );

    //double timeElapsed = this->timerTool("Solve").stop();
    //this->log( "MultiFluid", "updateFluidDensityViscosity", 
            //"fluid density/viscosity update in "+(boost::format("%1% s") %timeElapsed).str() );
//}

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

    this->fluidModel()->updateSourceAdded( idv(M_interfaceForces) );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log( "MultiFluid", "updateInterfaceForces", 
            "interface forces updated in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::advectLevelsets()
{
    this->log("MultiFluid", "advectLevelsets", "start");
    this->timerTool("Solve").start();

    auto const& u = this->fieldVelocity();
    
    for( auto const& ls: M_levelsetModels )
    {
        ls.second->advect( idv(u) );
    }
    // Request global levelset update
    M_globalLevelset.setDoUpdate( true );

    if( this->hasInextensibilityLM() )
    {
        M_doUpdateInextensibilityLM = true;
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log( "MultiFluid", "advectLevelsets", 
            "level-sets advection done in "+(boost::format("%1% s") %timeElapsed).str() );
}

} // namespace FeelModels
} // namespace Feel
