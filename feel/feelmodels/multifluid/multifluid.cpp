/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/multifluid/multifluid.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

namespace Feel {
namespace FeelModels {

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
MULTIFLUID_CLASS_TEMPLATE_TYPE::MultiFluid(
        std::string const& prefix,
        WorldComm const& wc,
        std::string const& subPrefix,
        std::string const& rootRepository )
: super_type( prefix, wc, subPrefix, self_type::expandStringFromSpec( rootRepository ) )
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

    M_fluid.reset( 
            new fluid_type( prefixvm(this->prefix(),"fluid"), false, this->worldComm(), "", this->rootRepositoryWithoutNumProc() ) 
            ); 
    M_globalLevelset.reset(
            new levelset_type( prefixvm(this->prefix(),"levelset"), this->worldComm(), "", this->rootRepositoryWithoutNumProc() )
            );

    this->createMesh();

    M_fluid->loadMesh( M_mesh );
    M_globalLevelset->build( this->mesh() );

    // "Deep" copy
    M_fluidDensityViscosityModel.reset( new densityviscosity_model_type( M_fluid->prefix() ) );
    M_fluidDensityViscosityModel->initFromSpace( M_fluid->densityViscosityModel()->dynamicViscositySpace() );
    M_fluidDensityViscosityModel->updateFromModelMaterials( M_fluid->modelProperties().materials() );

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
        M_levelsetDensityViscosityModels[i]->initFromMesh( this->mesh(), M_fluid->useExtendedDofTable() );
        M_levelsetDensityViscosityModels[i]->updateFromModelMaterials( M_levelsets[i]->modelProperties().materials() );

        if( Environment::vm().count( prefixvm(levelset_prefix, "interface-forces-model").c_str() ) )
        {
            M_levelsetInterfaceForcesModels[i].reset( 
                    interfaceforces_factory_type::instance().createObject( 
                        soption( _name="interface-forces-model", _prefix=levelset_prefix ) 
                        )
                    );
            M_levelsetInterfaceForcesModels[i]->build( levelset_prefix, M_levelsets[i] );

            M_hasInterfaceForcesModel = true;
        }
    }

    M_interfaceForces.reset( new element_levelset_vectorial_type(this->functionSpaceLevelsetVectorial(), "InterfaceForces") ); 

    this->log("MultiFluid", "build", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::createMesh()
{
    this->log("MultiFluid","createMesh", "start");
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    M_fluid->setMshfileStr( this->mshfileStr() );

    double tElapsed = this->timerTool("Constructor").stop("createMesh");
    this->log("MultiFluid","createMesh", (boost::format("finish in %1% s") %tElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::init()
{
    this->log("MultiFluid", "init", "start");

    M_fluid->init();
    M_globalLevelset->init();
    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
        M_levelsets[i]->init();
    }

    this->updateTime( this->timeStepBase()->time() );

    this->updateGlobalLevelset();

    this->log("MultiFluid", "init", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    M_nFluids = ioption( _name="nfluids", _prefix=this->prefix() );

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
    M_levelsetReinitEvery.resize(nLevelSets);
    for( uint16_type n = 0; n < nLevelSets; ++n )
    {
        auto levelset_prefix = prefixvm(this->prefix(), (boost::format( "levelset%1%" ) %(n+1)).str());
        M_levelsetReinitEvery[n] = ioption( _name="reinit-every", _prefix=levelset_prefix );
    }
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
    }
    // Update global levelset
    this->updateGlobalLevelset();

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","solve","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateTime( double time )
{
    // Fluid mechanics
    M_fluid->updateTime(time);
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
    // Fluid mechanics
    M_fluid->updateTimeStep();
    // Levelsets
    for( uint16_type i = 0; i < M_levelsets.size(); ++i)
    {
        M_levelsets[i]->updateTimeStep();
    }

    this->updateTime( this->timeStepBase()->time() );

    this->log("MultiFluid", "updateTimeStep", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("MultiFluid","exportResults", "start");

    M_fluid->exportResults(time);

    if( this->nLevelsets() > 1 )
        M_globalLevelset->exportResults(time);

    for( uint16_type i = 0; i < M_levelsets.size(); ++i)
    {
        M_levelsets[i]->exportResults(time);
    }

    this->log("MultiFluid", "exportResults", "finish");
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
            M_fluid->densityViscosityModel()->dynamicViscositySpace(),
            elements(M_fluid->mesh()),
            idv(M_fluidDensityViscosityModel->fieldRho())*idv(globalH)
            );

    auto mu = vf::project( 
            M_fluid->densityViscosityModel()->dynamicViscositySpace(),
            elements(M_fluid->mesh()),
            idv(M_fluidDensityViscosityModel->fieldMu())*idv(globalH)
            );

    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
        rho += vf::project( 
                M_fluid->densityViscosityModel()->dynamicViscositySpace(),
                elements(M_fluid->mesh()),
                idv(M_levelsetDensityViscosityModels[i]->fieldRho())*(1. - idv(M_levelsets[i]->H()))
                );
        mu += vf::project( 
                M_fluid->densityViscosityModel()->dynamicViscositySpace(),
                elements(M_fluid->mesh()),
                idv(M_levelsetDensityViscosityModels[i]->fieldMu())*(1. - idv(M_levelsets[i]->H()))
                );
    }

    M_fluid->updateRho( idv(rho) );
    M_fluid->updateMu( idv(mu) );

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
                    elements(M_fluid->mesh()),
                    - M_surfaceTensionCoeff(0,n+1)*idv(M_levelsets[n]->K())*idv(M_levelsets[n]->N())*idv(M_levelsets[n]->D())
                    );
        }
        double timeElapsedSurfaceTension = this->timerTool("Solve").stop();
        this->log("MultiFluid", "updateInterfaceForces", "update surface tension forces in "+(boost::format("%1% s")%timeElapsedSurfaceTension).str() );
    }

    if( M_hasInterfaceForcesModel )
    {
        this->timerTool("Solve").start();
        for( uint16_type n = 0; n < M_levelsets.size(); ++n )
        {
            if( M_levelsetInterfaceForcesModels[n] )
            {
                M_levelsetInterfaceForcesModels[n]->updateInterfaceForces( M_interfaceForces, false );
            }
        }
        double timeElapsedInterfaceForces = this->timerTool("Solve").stop();
        this->log("MultiFluid", "updateInterfaceForces", "update interface (model) forces in "+(boost::format("%1% s")%timeElapsedInterfaceForces).str() );
    }

    M_fluid->updateSourceAdded( idv(M_interfaceForces) );

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

    M_fluid->solve();

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

    auto u = M_fluid->fieldVelocity();
    
    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
        M_levelsets[i]->advect( idv(u) );
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log( "MultiFluid", "advectLevelsets", 
            "level-sets advection done in "+(boost::format("%1% s") %timeElapsed).str() );
}

} // namespace FeelModels
} // namespace Feel
