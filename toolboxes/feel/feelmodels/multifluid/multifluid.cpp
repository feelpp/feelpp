/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <numeric>
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
    ModelPhysics<mesh_type::nDim>( "multifluid" ),
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
    this->loadParametersFromOptionsVm();

    this->log("MultiFluid","constructor", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("MultiFluid","initMesh", "start");
    this->timerTool("Constructor").start();

    if ( this->modelProperties().jsonData().contains("Meshes") )
        super_type::super_model_meshes_type::setup( this->modelProperties().jsonData().at("Meshes"), {this->keyword()} );
    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( this->mesh() ) << "mesh generation failed";

    double tElapsed = this->timerTool("Constructor").stop("initMesh");
    this->log("MultiFluid","initMesh", fmt::format( "finish in {} s", tElapsed ) );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    CHECK( M_nLevelsets >= 1 ) << "Multifluid must contain at least 1 levelset.\n";
    this->log("MultiFluid", "init", "start");
    this->timerTool("Constructor").start();

    this->initModelProperties();

    M_fluidModel = std::make_shared<fluid_model_type>(
            prefixvm(this->prefix(),"fluid"), "fluid", this->worldCommPtr(),
            this->subPrefix(), this->repository() );
    //for( index_type i = 0; i < M_nLevelsets; ++i )
    //{
        //std::string const lsKeyword = fmt::format( "levelset{}", i );
        //std::string const lsPrefix = prefixvm( this->prefix(), lsKeyword );
        //M_levelsetModels.push_back( std::make_shared<levelset_model_type>(
                //lsPrefix, lsKeyword, this->worldCommPtr(),
                //this->subPrefix(), this->repository() )
                //);
    //}

    this->initPhysics(
        this->shared_from_this(),
        [this]( typename super_physics_type::PhysicsTree & physicsTree ) {
            physicsTree.updatePhysics( this->shared_from_this(), this->modelProperties().models() );
            physicsTree.updatePhysics( M_fluidModel, this->modelProperties().models() );
            for( auto const& lsModel: M_levelsetModels )
                physicsTree.updatePhysics( lsModel, this->modelProperties().models() );
        } );
    //using physic_id_type = typename model_physics_type::physic_id_type;
    //std::map< physic_id_type, std::shared_ptr<ModelPhysic<nDim>> > innerFluidsPhysics;
    //std::string fluidType = M_fluidModel->keyword();
    //for ( auto const& [physicId, physic]: this->physicsFromCurrentType() )
    //{
        //for( uint32_t n = 0; n < this->nFluids() - 1; ++n )
        //{
            //std::string fluidName = fmt::format( "fluid{}", n );
            //innerFluidsPhysics[physic_id_type( { "fluid", fluidName } )] = ModelPhysic<nDim>::New( M_fluidModel, "fluid", fluidType, fluidName, this->modelProperties().models() );
            //innerFluidsPhysics[fluidKeyword] = std::make_shared< ModelPhysics< nDim > >( "fluid", fluidKeyword );
            //physicsTree.addChild( M_modelPhysicsFluids[fluidKeyword], models );
        //}
        //M_fluidModel->addPhysics( );
    //}
    // Set multifluid density and viscosity for inner and outer fluids
    for ( auto const& [physicId, physic]: this->physicsFromCurrentType() )
    {
        for ( auto subphysic: physic->subphysicsFromType( M_fluidModel->physicType() ) )
        {
            auto subphysicFluid = std::dynamic_pointer_cast<ModelPhysicFluid<nDim>>(subphysic);
            CHECK( subphysicFluid ) << "invalid fluid physic ptr cast";
            subphysicFluid->dynamicViscosity()->setMultifluid( true );
        }
    }

    // Physical properties
    if ( !M_materialsProperties )
    {
        //auto paramValues = this->modelProperties().parameters().toParameterValues();
        //this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    // Mesh
    if ( !this->mesh() )
        this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    // Init fluid toolbox
    this->initFluidToolbox();
    // Init levelset toolboxes
    this->initLevelsetToolboxes();

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
    for ( levelset_model_ptrtype const& lsModel: M_levelsetModels )
        lsModel->setAdvectionVelocityExpr( velConvExprStr );

    // Update constant parameters
    this->updateParameterValues();

    // Update initial conditions
    this->updateInitialConditions( this->symbolsExpr() );

    // Backend
    this->initAlgebraicBackend();

    // Block vector solution
    this->buildBlockVectorSolution();

    // Cache levelset models block vector solutions
    M_algebraicBlockVectorSolutionLevelsets.clear();
    std::transform(
            M_levelsetModels.begin(), M_levelsetModels.end(), 
            std::back_inserter( M_algebraicBlockVectorSolutionLevelsets ),
            []( levelset_model_ptrtype const& lsModel ) { return lsModel->algebraicBlockVectorSolution()->vectorMonolithic(); } 
            );

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
                [this]( DataUpdateJacobian & data ) { updateJacobian_Fluid( data ); }
                );
        M_fluidModel->algebraicFactory()->setFunctionResidualAssembly(
                [this]( DataUpdateResidual & data ) { updateResidual_Fluid( data ); }
                );
        for( size_type i = 0; i < M_levelsetModels.size(); ++i )
        {
            M_levelsetModels[i]->algebraicFactory()->setFunctionLinearAssembly(
                    [this, i]( DataUpdateLinear & data ) { updateLinear_Levelset( i, data ); }
                    );
            M_levelsetModels[i]->algebraicFactory()->setFunctionJacobianAssembly(
                    [this, i]( DataUpdateJacobian & data ) { updateJacobian_Levelset( i, data ); }
                    );
            M_levelsetModels[i]->algebraicFactory()->setFunctionResidualAssembly(
                    [this, i]( DataUpdateResidual & data ) { updateResidual_Levelset( i, data ); }
                    );
        }
    }

    //// "Deep" copy FluidMechanics materialProperties
    //M_fluidMaterialProperties.reset( new materialsproperties_type( this->fluidModel()->prefix() ) );
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

    M_enableInextensibility = false;
    for( index_type i = 0; i < M_nLevelsets; ++i )
    {
        std::string const lsKeyword = fmt::format( "levelset{}", i );
        std::string const lsPrefix = prefixvm( this->prefix(), lsKeyword );

        // Inextensibility
        bool hasInextensibility = boption( _name="enable-inextensibility", _prefix=lsPrefix );
        M_hasInextensibility.push_back( hasInextensibility );

        if( hasInextensibility ) M_enableInextensibility = true;

        std::string inextensibilityMethod = soption( _name="inextensibility-method", _prefix=lsPrefix );
        CHECK( inextensibilityMethod == "penalty" || inextensibilityMethod == "lagrange-multiplier" ) 
            << "invalid inextensiblity-method " << inextensibilityMethod
            << ", should be \"penalty\" or \"lagrange-multiplier\"" << std::endl;
        M_inextensibilityMethod.push_back( inextensibilityMethod );

        if( hasInextensibility && inextensibilityMethod == "lagrange-multiplier" )
            M_hasInextensibilityLM = true;

        double inextensibilityGamma = doption( _name="inextensibility-gamma", _prefix=lsPrefix );
        M_inextensibilityGamma.push_back( inextensibilityGamma );

        // Redistanciation
        int redistEvery = ioption( _name="redist-every", _prefix=lsPrefix );
        M_levelsetRedistEvery.push_back( redistEvery );
    }
}

#if 0
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

    //*_ostr << "\n   Fluids Parameters";
    //*_ostr << "\n     -- fluid 0 (outer fluid)"
           //<< "\n       * rho : " << this->M_fluidMaterialProperties->cstDensity()
           //<< "\n       * mu  : " << this->M_fluidMaterialProperties->cstMu();
    //*_ostr << this->M_fluidMaterialProperties->getInfo()->str();
    //for( auto const& lsMaterialProperty: M_levelsetsMaterialProperties )
    //{
    //*_ostr << "\n     -- fluid " << "\"" << lsMaterialProperty.first << "\""
           //<< "\n       * rho : " << lsMaterialProperty.second->cstDensity()
           //<< "\n       * mu  : " << lsMaterialProperty.second->cstMu();
    //*_ostr << lsMaterialProperty.second->getInfo()->str();
    //}

    //*_ostr << "\n   Level Sets Parameters";
    //for( index_type i = 0; i < M_nLevelsets; ++i )
    //{
    //*_ostr << "\n     -- level set " << "\"" << i << "\""
           //<< "\n       * redist every : " << M_levelsetRedistEvery[i];
    //}

    //*_ostr << "\n   Forces Parameters";
    //*_ostr << "\n     -- has interface forces : " << std::boolalpha << this->M_hasInterfaceForcesModel;
    //if( this->M_hasInterfaceForcesModel )
    //{
        //for( auto const& lsInterfaceForces: M_levelsetInterfaceForcesModels )
        //{
    //*_ostr << "\n     -- level set " << "\"" << lsInterfaceForces.first << "\"";
            //for( auto const& force: lsInterfaceForces.second )
    //*_ostr << "\n       * force model : " << force.second->getInfo()->str();
        //}
    //}

    //*_ostr << "\n";
    //*_ostr << this->fluidModel()->getInfo()->str();
    //for( auto const& levelset: M_levelsetModels )
    //{
    //*_ostr << levelset->getInfo()->str();
    //}

    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n\n";

    return _ostr;
}
#endif

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    //TODO
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
MULTIFLUID_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    //TODO
    return tabInfo;
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
    this->log("MultiFluid", "updateGlobalLevelset", fmt::format( "finish in {} s", timeElapsed ) );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
int
MULTIFLUID_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = this->fluidModel()->nBlockMatrixGraph();
    for( levelset_model_ptrtype const& lsModel: this->levelsetModels() )
        nBlock += lsModel->nBlockMatrixGraph();
    if( this->hasInextensibilityLM() )
        nBlock += 1;

    return nBlock;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
MULTIFLUID_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("MultiFluid","buildBlockMatrixGraph", "start" );

    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myBlockGraph(nBlock, nBlock);

    int nBlockFluid = this->fluidModel()->nBlockMatrixGraph();

    int startIndexBlockFluid = 0;
    int indexBlock = startIndexBlockFluid;

    auto blockMatFluid = this->fluidModel()->buildBlockMatrixGraph();
    for (int tk1=0; tk1<nBlockFluid; ++tk1 )
        for (int tk2=0; tk2<nBlockFluid; ++tk2 )
            myBlockGraph(startIndexBlockFluid+tk1,startIndexBlockFluid+tk2) = blockMatFluid(tk1,tk2);
    indexBlock += nBlockFluid;

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
                _test=this->functionSpaceInextensibilityLM(), _trial=this->fluidModel()->functionSpaceVelocity(),
                //_pattern_block=patCouplingLM,
                _diag_is_nonzero=false,_close=false)->graph();
        myBlockGraph(startIndexBlockFluid,startIndexBlockInextensibility) = stencil(
                _test=this->fluidModel()->functionSpaceVelocity(), _trial=this->functionSpaceInextensibilityLM(),
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
bool
MULTIFLUID_CLASS_TEMPLATE_TYPE::hasInextensibilityLM() const
{
    return M_hasInextensibilityLM;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateInextensibilityLM()
{
    this->log("MultiFluid", "updateInextensibilityLM", "start");
    M_inextensibleLevelsets.clear();
    for( size_type i = 0; i < M_nLevelsets; ++i )
    {
        if( this->hasInextensibility(i) )
        {
            if( this->inextensibilityMethod(i) == "lagrange-multiplier" )
            {
                M_inextensibleLevelsets.push_back( this->levelsetModel(i)->phi() );
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
            _expr=Feel::FeelModels::levelsetDelta( _element=inextensibleLevelsets, _thickness=M_globalLevelsetThicknessInterface )
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
    //TODO:M_doRebuildMatrixVector = true;
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
        forceName = fmt::format( "force{}", M_additionalInterfaceForcesModel.size() );
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
    for( size_type i = 0; i < this->levelsetModels().size(); ++i )
    {
        if( M_levelsetRedistEvery[i] > 0 
                && (this->levelsetModel(i)->iterSinceRedistanciation()+1) % M_levelsetRedistEvery[i] == 0 )
            this->levelsetModel(i)->redistanciate();
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid", "solve", fmt::format( "finish in {} s", timeElapsed ) );
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
        //// Rebuild matrix and vector
        //if( this->M_doRebuildMatrixVector )
        //{
            //// Rebuild algebraic matrix and vector
            //auto graph = this->buildMatrixGraph();
            //M_algebraicFactory->rebuildMatrixVector( graph, graph->mapRow().indexSplit() );
            //// Rebuild solution vector
            //this->buildBlockVectorSolution();
            //// Rebuild backend (required since PETSc stores the size of the matrix, 
            //// which can change here because of the LM space support)
            //this->backend()->clear();

            //M_doRebuildMatrixVector = false;
        //}

        this->fluidModel()->setStartBlockSpaceIndex( this->startSubBlockSpaceIndex("fluid") );

        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();

        // Solve fluid
        this->algebraicFactory()->solve( this->fluidModel()->solverName(), this->algebraicBlockVectorSolution()->vectorMonolithic() );
        this->algebraicBlockVectorSolution()->localize();
    }
    else
    {
        // Solve fluid equations (with direct assembly of interface forces)
        M_fluidModel->solve();
        for( auto const& lsModel: M_levelsetModels )
            lsModel->solve();

        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
    }
    // Request global levelset update
    M_globalLevelset.setDoUpdate( true );
    // Possibly request inextensiblity LM update
    if( this->hasInextensibilityLM() )
        M_doUpdateInextensibilityLM = true;
    //// Update density and viscosity
    //this->updateFluidDensityViscosity();
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

    this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
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
        u_old = this->fluidModel()->fieldVelocity();
        // Sub-solves
        if( M_doUpdateInextensibilityLM )
            this->updateInextensibilityLM();
        this->solveSemiImplicitCoupling();
        auto u = this->fluidModel()->fieldVelocity();

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
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateParameterValues()
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

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
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

    M_fluidModel->setParameterValues( paramValues );
    for( auto const& lsModel: M_levelsetModels )
        lsModel->setParameterValues( paramValues );
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
    for( auto const& lsModel: this->levelsetModels() )
        lsModel->startTimeStep();

    this->updateTime( this->timeStepBase()->time() );

    this->updateParameterValues();

    this->log("MultiFluid", "startTimeStep", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->log("MultiFluid", "updateTimeStep", "start");

    this->fluidModel()->updateTimeStep();
    for( auto const& lsModel: this->levelsetModels() )
        lsModel->updateTimeStep();

    this->updateTime( this->timeStepBase()->time() );

    this->updateParameterValues();

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
    for( auto const& lsModel: M_levelsetModels)
        lsModel->exportResults( time, symbolsExpr );

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
    //// Levelset forces
    //for( std::string const& levelsetName: M_postProcessMeasuresLevelsetForces )
    //{
        //auto measuredLevelsetForce = this->computeLevelsetForce( levelsetName );
        //std::vector<double> vecMeasuredLevelsetForce = { measuredLevelsetForce(0,0) };
        //if( nDim > 1 ) vecMeasuredLevelsetForce.push_back( measuredLevelsetForce(1,0) );
        //if( nDim > 2 ) vecMeasuredLevelsetForce.push_back( measuredLevelsetForce(2,0) );
        //this->postProcessMeasuresIO().setMeasureComp( levelsetName + ".force", vecMeasuredLevelsetForce );
    //}
    //TODO
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
typename MULTIFLUID_CLASS_TEMPLATE_TYPE::force_type
MULTIFLUID_CLASS_TEMPLATE_TYPE::computeLevelsetForce( std::string const& name ) const
{
#if 0
    auto const& u = this->fluidModel()->fieldVelocity();
    auto const& p = this->fluidModel()->fieldPressure();
    std::string matName = *(this->fluidModel()->materialsProperties()->physicToMaterials( this->fluidModel()->physicType() ).begin());
    auto sigmav = Feel::FeelModels::fluidMecStressTensor(gradv(u),idv(p),*this->fluidModel()->materialsProperties(),matName,true,2*fluid_model_type::nOrderVelocity,true);

    auto const& phi = this->levelsetModel(name)->phi();
    auto N_expr = trans(gradv(phi)) / sqrt( gradv(phi) * trans(gradv(phi)) );
    auto const D_expr = this->levelsetModel(name)->diracExpr();
    return integrate(_range=elements(this->mesh()),
                     _expr= - sigmav * N_expr * D_expr,
                     _geomap=this->geomap() ).evaluate();
#endif
    CHECK( false ) << "TODO\n";
    return force_type();
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::initFluidToolbox()
{
    this->log("MultiFluid", "initFluidToolbox", "start");

    M_fluidModel->setManageParameterValues( false );
    if ( !M_fluidModel->modelPropertiesPtr() )
    {
        M_fluidModel->setModelProperties( this->modelPropertiesPtr() );
        M_fluidModel->setManageParameterValuesOfModelProperties( false );
    }
    M_fluidModel->setModelMeshAsShared( this->modelMesh() );
    M_fluidModel->setMaterialsProperties( M_materialsProperties );
    M_fluidModel->init( false );

    this->log("MultiFluid", "initLevelsets", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::initLevelsetToolboxes()
{
    this->log("MultiFluid", "initLevelsetToolboxes", "start");
    // Get levelset mesh
    mesh_ptrtype mesh;
    if( this->M_useLagrangeP1iso )
    {
        // Build Lagrange P1 iso-U mesh and build levelsets with it
        M_opLagrangeP1iso = lagrangeP1( _space=this->fluidModel()->functionSpaceVelocity()->compSpace() );
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
    M_levelsetSpaceManager->initFunctionSpaceDefault();
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

    //TODO
    std::set<std::string> ppExportsAllFieldsAvailable;
    for ( auto const& s : M_fluidModel->postProcessExportsAllFieldsAvailable() )
        ppExportsAllFieldsAvailable.insert( prefixvm( M_fluidModel->keyword(), s ) );
    for ( auto const& lsModel : M_levelsetModels )
        for ( auto const& s : lsModel->postProcessExportsAllFieldsAvailable() )
            ppExportsAllFieldsAvailable.insert( prefixvm( lsModel->keyword(), s ) );
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

#if 0
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
#endif

    auto se = this->symbolsExpr();
    this->template initPostProcessMeshes<mesh_type>( se );

    // start or restart the export of measures
    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasures().restart( this->timeInitial() );
    }

    double tElapsed = this->timerTool("Constructor").stop("initPostProcess");
    this->log("MultiFluid","initPostProcess",fmt::format("finish in {} s", tElapsed ) );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updatePhysics( typename super_physics_type::PhysicsTreeNode & physicsTree, ModelModels const& models )
{
    auto currentPhysic = std::dynamic_pointer_cast<ModelPhysicMultifluid<nDim>>( physicsTree.physic() );
    CHECK( currentPhysic ) << "wrong physic: expected multifluid, got " << physicsTree.physic()->name();

    physicsTree.addChild( M_fluidModel, models );

    //for( uint32_t n = 0; n < currentPhysic->nFluids() - 1; ++n )
    //{
        //std::string fluidKeyword = fmt::format( "fluid{}", n );
        //M_modelPhysicsFluids[fluidKeyword] = std::make_shared< ModelPhysics< nDim > >( "fluid", fluidKeyword );
        //physicsTree.addChild( M_modelPhysicsFluids[fluidKeyword], models );
    //}

    for( uint32_t n = 0; n < currentPhysic->nFluids() - 1; ++n )
    {
        std::string const levelsetKeyword = fmt::format( "levelset{}", n );
        std::string const levelsetPrefix = prefixvm( this->prefix(), levelsetKeyword );
        levelset_model_ptrtype levelsetModel = std::make_shared<levelset_model_type>(
                levelsetPrefix, levelsetKeyword, this->worldCommPtr(),
                this->subPrefix(), this->repository() );
        M_levelsetModels.push_back( levelsetModel );
        physicsTree.addChild( levelsetModel, models );
    }

    //physicsTree.updateMaterialSupportFromChildren( "intersect" );
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
            this->levelsetModels().begin(), this->levelsetModels().end(), 0,
                []( int res, auto const& lsModel ) { return res + lsModel->algebraicBlockVectorSolution()->size(); }
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
        this->log("MultiFluid", "updateInterfaceForces", fmt::format( "update interface (model) forces in {} s", timeElapsedInterfaceForces ) );
    }

    if( M_additionalInterfaceForcesModel.size() > 0 )
    {
        this->timerTool("Solve").start();
        for( auto const& f: M_additionalInterfaceForcesModel )
            f.second->updateInterfaceForces( M_interfaceForces, false );

        double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
        this->log("MultiFluid", "updateInterfaceForces", fmt::format( "update additional interface forces in {} s", timeElapsedInterfaceForces ) );
    }

    this->fluidModel()->updateSourceAdded( idv(M_interfaceForces) );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log( "MultiFluid", "updateInterfaceForces", 
            fmt::format( "interface forces updated in {} s", timeElapsed ) );
}

} // namespace FeelModels
} // namespace Feel
