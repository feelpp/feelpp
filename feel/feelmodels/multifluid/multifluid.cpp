/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/multifluid/multifluid.hpp>

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
    CHECK( M_nFluids > 2 ) << "Multifluid must contain at least 2 fluids.\n";
    uint16_type nLevelSets = M_nFluids - 1;

    this->log("MultiFluid", "build", "start");

    M_fluid = fluid_ptrtype( 
            new fluid_type("fluid", false, this->worldComm(), "", this->rootRepositoryWithoutNumProc() ) 
            ); 
    M_fluid->build();

    M_fluidDensityViscosityModel.reset( new densityviscosity_model_type(*M_fluid->densityViscosityModel()) );

    M_globalLevelset = levelset_ptrtype(
            new levelset_type( "levelset", this->worldComm(), "", this->rootRepositoryWithoutNumProc() )
            );
    M_globalLevelset->build( M_fluid->mesh() );

    M_levelsets.resize( nLevelSets );
    M_levelsetDensityViscosityModels.resize( nLevelSets );
    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
        auto levelset_prefix = (boost::format( "levelset%1%" ) %i).str();
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
                _smoother_curvature=M_globalLevelset->smootherCurvature()
                );

        M_levelsetDensityViscosityModels[i].reset(
                new densityviscosity_model_type( levelset_prefix )
                );
        M_levelsetDensityViscosityModels[i]->initFromMesh( M_fluid->mesh(), M_fluid->useExtendedDofTable() );
        M_levelsetDensityViscosityModels[i]->updateFromModelMaterials( M_levelsets[i]->modelProperties().materials() );
    }

    M_interfaceForces.reset( new element_levelset_vectorial_type(this->functionSpaceLevelsetVectorial(), "InterfaceForces") ); 

    this->log("MultiFluid", "build", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    M_nFluids = ioption( _name="nfluids", _prefix=this->prefix() );

    M_enableSurfaceTension = boption( _name="enable-surface-tension", _prefix=this->prefix() );

    if( M_enableSurfaceTension )
    {
        std::vector<double> sigma = Environment::vm()[prefixvm(this->prefix(),"surface-tension-coeff").c_str()].template as<std::vector<double> >();

        CHECK( sigma.size() >= M_nFluids - 1 ) << sigma.size() << " surface tension coefficients found.\n"
                                               << "You must at least provide the surface tension coefficients between the "
                                               << M_nFluids - 1
                                               << " levelset fluids and the surrounding fluid.\n";
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
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
bool
MULTIFLUID_CLASS_TEMPLATE_TYPE::hasInterfaceForces() const
{
    return this->hasSurfaceTension();
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
        this->updateInterfaceForces();
    // Solve fluid equations
    M_fluid->solve();
    // Advect levelsets
    this->advectLevelsets();

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("MultiFluid","solve","finish in "+(boost::format("%1% s") %timeElapsed).str() );
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
    this->log("MultiFluid", "updateTimeStep", "finish");
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::updateGlobalLevelset()
{
    this->log("MultiFluid", "updateGlobalLevelset", "start");

    auto minPhi = M_globalLevelset->phi();

    *minPhi = *M_levelsets[0]->phi();
    for( uint16_type i = 1; i < M_levelsets.size(); ++i )
    {
        *minPhi = vf::project( 
                M_globalLevelset->functionSpace(), 
                elements(M_globalLevelset->mesh()),
                vf::min( idv(minPhi), idv(M_levelsets[i]->phi()) )
                );
    }

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
                idv(M_levelsetDensityViscosityModels[i]->fieldRho())*(1 - idv(M_levelsets[i]->H()))
                );
        mu += vf::project( 
                M_fluid->densityViscosityModel()->dynamicViscositySpace(),
                elements(M_fluid->mesh()),
                idv(M_levelsetDensityViscosityModels[i]->fieldMu())*(1 - idv(M_levelsets[i]->H()))
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
    M_interfaceForces->zero();

    if( this->hasSurfaceTension() )
    {
        for( uint16_type n = 0; n < M_levelsets.size(); ++n )
        {
            *M_interfaceForces += vf::project( 
                    this->functionSpaceLevelsetVectorial(),
                    elements(this->mesh()),
                    - M_surfaceTensionCoeff(0,n+1)*idv(M_levelsets[n]->K())*idv(M_levelsets[n]->N())*idv(M_levelsets[n]->D())
                    );
        }
    }

    M_fluid->updateSourceAdded( idv(M_interfaceForces) );
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
     this->updateGlobalLevelset();

    double timeElapsed = this->timerTool("Solve").stop();
    this->log( "MultiFluid", "advectLevelsets", 
            "level-sets advection done in "+(boost::format("%1% s") %timeElapsed).str() );
}

} // namespace FeelModels
} // namespace Feel
