#include <feel/feelmodels/levelset/levelset.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

#include <feel/feelmodels/levelset/cauchygreentensorexpr.hpp>
#include <feel/feelmodels/levelset/cauchygreeninvariantsexpr.hpp>
#include <feel/feelmodels/levelset/levelsetdeltaexpr.hpp>
#include <feel/feelmodels/levelset/levelsetheavisideexpr.hpp>
//#include <feel/feelmodels/levelset/levelsetcurvatureexpr.hpp>

#include <boost/assign/list_of.hpp>
#include <boost/assign/list_inserter.hpp>

namespace Feel {
namespace FeelModels {

//----------------------------------------------------------------------------//
// Static member initialization

//----------------------------------------------------------------------------//

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
LEVELSET_CLASS_TEMPLATE_TYPE::LevelSet( 
        std::string const& prefix,
        std::string const& keyword,
        worldcomm_ptr_t const& worldComm,
        std::string const& subPrefix,
        ModelBaseRepository const& modelRep ) 
:
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep ),
    M_advectionToolbox( new cfpde_toolbox_type( 
                typename FeelModels::ModelGenericPDE<nDim>::infos_type( "levelset_phi", "phi", "phi", (boost::format("Pch%1%")%Order).str() ),
                prefix, keyword, worldComm, subPrefix, modelRep ) ),
    M_doUpdateCauchyGreenTensor(true),
    M_doUpdateCauchyGreenInvariant1(true),
    M_doUpdateCauchyGreenInvariant2(true),
    M_iterSinceRedistanciation(0),
    M_useOrder1AfterRedist( false )
{
    //-----------------------------------------------------------------------------//
    // Set advection model
#if 0
    M_advectionToolbox->setModelName( "Advection" );
#endif
    //-----------------------------------------------------------------------------//
    // Load parameters
    this->loadParametersFromOptionsVm();
    // Get periodicity from options (if needed)
    //this->loadPeriodicityFromOptionsVm();

    /*// --------------- mesh adaptation -----------------
#if defined (MESH_ADAPTATION)
    auto backend_mesh_adapt = backend_type::build(Environment::vm(), "mesh-adapt-backend");
    mesh_adapt.reset( new mesh_adaptation_type ( backend_mesh_adapt ));
#endif*/
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::self_ptrtype
LEVELSET_CLASS_TEMPLATE_TYPE::New(
        std::string const& prefix,
        std::string const& keyword,
        worldcomm_ptr_t const& worldComm,
        std::string const& subPrefix,
        ModelBaseRepository const& modelRep )
{
    self_ptrtype new_ls( new self_type(prefix, keyword, worldComm, subPrefix, modelRep ) );
    return new_ls;
}

//----------------------------------------------------------------------------//
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("LevelSet", "init", "start");

    // Build extended dof table if CIP stabilisation
#if 0
    if( this->M_advectionToolbox->stabilizationMethod() == AdvectionStabMethod::CIP )
        this->setBuildExtendedDofSpace( true );
#endif
    // Init LevelSetBase
    super_type::init();

    // Function spaces
    this->initFunctionSpaces();
    // Tools
    this->initInterfaceQuantities();
    this->initTools();

    // Init advection toolbox
    M_advectionToolbox->setSpaceUnknown( this->functionSpace() );
    M_advectionToolbox->init( buildModelAlgebraicFactory );
    // Set transient coefficient d=1
    for( std::string const& matName : M_advectionToolbox->materialsProperties()->physicToMaterials( M_advectionToolbox->physicDefault() ) )
    {
        auto & matProp = M_advectionToolbox->materialsProperties()->materialProperties( matName );
        ModelExpression firstTimeDerivativeCoefficientExpr;
        firstTimeDerivativeCoefficientExpr.setExpr( "1", this->worldComm(), this->repository().expr() );
        M_advectionToolbox->materialsProperties()->addProperty( matProp, M_advectionToolbox->firstTimeDerivativeCoefficientName(), firstTimeDerivativeCoefficientExpr, true );
    }

    if ( !this->isStationary() )
    {
        this->setTimeInitial( M_advectionToolbox->timeInitial() );
        this->setTimeOrder( M_advectionToolbox->timeOrder() );
        this->updateTime( M_advectionToolbox->currentTime() );
    }

    // Set advection initial value
    if( !this->doRestart() )
    {
        // Set other initial values
        this->updateInitialConditions();
    }
    // Init boundary conditions
    this->initBoundaryConditions();

    // We init the exporters only after having set the correct initial time
    // to prevent restart mishmash
    this->initPostProcessExporters();
    this->initPostProcessMeasures();

    // Init iterSinceRedistanciation
    if( this->doRestart() )
    {
        // Reload saved iterSinceRedistanciation data
        auto iterSinceRedistanciationPath = fs::path(this->rootRepository()) / fs::path( prefixvm(this->prefix(), "itersincereinit") );
        if( fs::exists( iterSinceRedistanciationPath ) )
        {
            fs::ifstream ifs( iterSinceRedistanciationPath );
            boost::archive::text_iarchive ia( ifs );
            ia >> BOOST_SERIALIZATION_NVP( M_vecIterSinceRedistanciation );
            M_iterSinceRedistanciation = M_vecIterSinceRedistanciation.back();
        }
        else
        {
            // If iterSinceRedistanciation not found, we assume that last step redistanciated by default
            M_iterSinceRedistanciation = 0;
        }
    }
    else
    {
            M_vecIterSinceRedistanciation.push_back( M_iterSinceRedistanciation );
    }
    // Adjust BDF order with iterSinceRedistanciation
    if( this->useOrder1AfterRedist() && M_iterSinceRedistanciation < this->timeOrder() )
    {
        this->timeStepBDF()->setTimeOrder( M_iterSinceRedistanciation + 1 );
        if( M_useGradientAugmented )
            M_modGradPhiAdvection->timeStepBdfUnknown()->setTimeOrder( M_iterSinceRedistanciation + 1 );
        if( M_useStretchAugmented )
            M_stretchAdvection->timeStepBdfUnknown()->setTimeOrder( M_iterSinceRedistanciation + 1 );
        if( M_useCauchyAugmented )
            M_backwardCharacteristicsAdvection->timeStepBdfUnknown()->setTimeOrder( M_iterSinceRedistanciation + 1 );
    }

    this->log("LevelSet", "init", "finish");
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateInitialConditions()
{
    this->log("LevelSet", "updateInitialConditions", "start");

    if( !this->doRestart() )
    {
        if( this->initialValue() )
        {
            // We need to hack the InitialCondition mechanism since it does not allow for proper management
            if( !this->isStationary() )
            {
                for( auto const& unknown : M_advectionToolbox->timeStepBdfUnknown()->unknowns() )
                    *unknown = *(this->initialValue());
            }
        }
    }

    if( M_useGradientAugmented )
    {
        // Initialize modGradPhi
        if( this->redistInitialValue() )
        {
            for( auto const& unknown : M_modGradPhiAdvection->timeStepBdfUnknown()->unknowns() )
                unknown->setConstant(1.);
        }
        else
        {
            *(M_modGradPhiAdvection->fieldUnknownPtr()) = this->modGrad( this->phiElt() );
            for( auto const& unknown : M_modGradPhiAdvection->timeStepBdfUnknown()->unknowns() )
                *unknown = M_modGradPhiAdvection->fieldUnknown();
        }
    }
    if( M_useStretchAugmented )
    {
        // Initialize stretch modGradPhi
        M_stretchAdvection->fieldUnknownPtr()->setConstant(1.);
        for( auto const& unknown : M_stretchAdvection->timeStepBdfUnknown()->unknowns() )
            unknown->setConstant(1.);
    }
    if( M_useCauchyAugmented )
    {
        // Initialize backward characteristics
        if( M_hasInitialBackwardCharacteristics )
        {
            M_backwardCharacteristicsAdvection->fieldUnknownPtr()->on(
                    _range=M_backwardCharacteristicsAdvection->rangeMeshElements(),
                    _expr=M_initialBackwardCharacteristics
                    );
        }
        else
        {
            M_backwardCharacteristicsAdvection->fieldUnknownPtr()->on(
                    _range=M_backwardCharacteristicsAdvection->rangeMeshElements(),
                    _expr=vf::P()
                    );
        }
        for( auto const& unknown : M_backwardCharacteristicsAdvection->timeStepBdfUnknown()->unknowns() )
            *unknown = M_backwardCharacteristicsAdvection->fieldUnknown();
    }

    this->log("LevelSet", "updateInitialConditions", "finish");
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    super_type::initPostProcessExportsAndMeasures();
    // Do not call the initPostProcess* functions here, this is done after time initial
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    if( !M_spaceAdvectionVelocity )
        M_spaceAdvectionVelocity = space_vectorial_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );

    if( M_useCauchyAugmented )
    {
        this->functionSpaceManager()->initFunctionSpaceTensor2Symm();
        M_spaceTensor2Symm = this->functionSpaceManager()->functionSpaceTensor2Symm();
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::initInterfaceQuantities()
{
    if( M_useGradientAugmented )
    {
        M_modGradPhiAdvection.reset( new modgradphi_advection_type( 
                    typename FeelModels::ModelGenericPDE<nDim>::infos_type( prefixvm( this->keyword(), "modgradphi", "_" ), "modgradphi", "modgradphi", (boost::format("Pch%1%")%Order).str() ),
                    this->prefix(), this->keyword(), this->worldCommPtr(), this->subPrefix(), this->rootRepository() ) );
#if 0
        M_modGradPhiAdvection->setModelName( "Advection-Reaction" );
        M_modGradPhiAdvection->setTimeOrder( this->timeOrder() );
#endif
        M_modGradPhiAdvection->setSpaceUnknown( this->functionSpace() );
        M_modGradPhiAdvection->init();
    }
    if( M_useStretchAugmented )
    {
        M_stretchAdvection.reset( new stretch_advection_type( 
                    typename FeelModels::ModelGenericPDE<nDim>::infos_type( prefixvm( this->keyword(), "stretch", "_" ), "stretch", "stretch", (boost::format("Pch%1%")%Order).str() ),
                    this->prefix(), this->keyword(), this->worldCommPtr(), this->subPrefix(), this->rootRepository() ) );
#if 0
        M_stretchAdvection->setModelName( "Advection-Reaction" );
        M_stretchAdvection->setTimeOrder( this->timeOrder() );
#endif
        M_stretchAdvection->setSpaceUnknown( this->functionSpace() );
        M_stretchAdvection->init();
    }
    if( M_useCauchyAugmented )
    {
        M_backwardCharacteristicsAdvection.reset( new backwardcharacteristics_advection_type( 
                    typename FeelModels::ModelGenericPDE<nDim>::infos_type( prefixvm( this->keyword(), "backward_characteristics", "_" ), "backward_characteristics", "backward_characteristics", (boost::format("Pchv%1%")%Order).str() ),
                    this->prefix(), this->keyword(), this->worldCommPtr(), this->subPrefix(), this->rootRepository() ) );
#if 0
        M_backwardCharacteristicsAdvection->setModelName( "Advection" );
        M_backwardCharacteristicsAdvection->setTimeOrder( this->timeOrder() );
#endif
        M_backwardCharacteristicsAdvection->setSpaceUnknown( this->functionSpaceVectorial() );
        M_backwardCharacteristicsAdvection->init();

        M_leftCauchyGreenTensor.reset( new element_tensor2symm_type(this->functionSpaceTensor2Symm(), "LeftCauchyGreenTensor") );
        M_cauchyGreenInvariant1.reset( new element_cauchygreen_invariant_type(this->functionSpace(), "CauchyGreenI1(TrC)") );
        M_cauchyGreenInvariant2.reset( new element_cauchygreen_invariant_type(this->functionSpace(), "CauchyGreenI2(TrCofC)") );
    }
#if 0
    for( std::string const& bcMarker: this->M_bcMarkersInflow )
    {
        if( M_useGradientAugmented )
            M_modGradPhiAdvection->addMarkerInflowBC( bcMarker );
        if( M_useStretchAugmented )
            M_stretchAdvection->addMarkerInflowBC( bcMarker );
        if( M_useCauchyAugmented )
            M_backwardCharacteristicsAdvection->addMarkerInflowBC( bcMarker );
    }
#endif
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::initTools()
{
    if( M_useCauchyAugmented )
    {
        this->toolManager()->initProjectorL2Tensor2Symm();
        M_projectorL2Tensor2Symm = this->toolManager()->projectorL2Tensor2Symm();
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_stretch_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::stretch() const
{
    if( !M_levelsetStretch )
        M_levelsetStretch.reset( new element_stretch_type(this->M_stretchAdvection->spaceUnknown(), "Stretch") );

    if( M_useStretchAugmented )
    {
        *M_levelsetStretch = *(M_stretchAdvection->fieldUnknownPtr());
    }
    else
    {
        *M_levelsetStretch = *(this->modGradPhi());
    }

    M_levelsetStretch->add(-1.);

    return M_levelsetStretch;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_backwardcharacteristics_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::backwardCharacteristics() const
{
    if( M_useCauchyAugmented )
    {
        return M_backwardCharacteristicsAdvection->fieldUnknownPtr();
    }
    else
    {
        throw std::logic_error( this->prefix()+".use-cauchy-augmented option must be true to use backward characteristics" );
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    M_useOrder1AfterRedist = boption( _name="use-order1-after-redist", _prefix=this->prefix() );

    M_useGradientAugmented = boption( _name="use-gradient-augmented", _prefix=this->prefix() );
    M_reinitGradientAugmented = boption( _name="reinit-gradient-augmented", _prefix=this->prefix() );

    M_useStretchAugmented = boption( _name="use-stretch-augmented", _prefix=this->prefix() );
    M_reinitStretchAugmented = boption( _name="reinit-stretch-augmented", _prefix=this->prefix() );

    M_useCauchyAugmented = boption( _name="use-cauchy-augmented", _prefix=this->prefix() );
    if ( Environment::vm().count(prefixvm(this->prefix(),"initial-backward-characteristics").c_str()) )
    {
        M_initialBackwardCharacteristics = expr<nDim,1>( soption(_name="initial-backward-characteristics", _prefix=this->prefix()) );
        M_hasInitialBackwardCharacteristics = true;
    }
    else
    {
        M_hasInitialBackwardCharacteristics = false;
    }

    M_doExportAdvection = boption(_name="do_export_advection", _prefix=this->prefix());

    M_fixVolume = boption( _name="fix-volume", _prefix=this->prefix() );
    M_fixArea = boption( _name="fix-area", _prefix=this->prefix() );

    M_useExtensionVelocity = boption( _name="use-extension-velocity", _prefix=this->prefix() );
    M_extensionVelocityNitscheGamma = doption( _name="extension-velocity.gamma", _prefix=this->prefix() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    M_bcMarkersInflow.clear();
    M_levelsetParticleInjectors.clear();

    // Inflow BC
    for( std::string const& bcMarker: this->modelProperties().boundaryConditions().markers( this->prefix(), "inflow" ) )
    {
        if( std::find(M_bcMarkersInflow.begin(), M_bcMarkersInflow.end(), bcMarker) == M_bcMarkersInflow.end() )
            M_bcMarkersInflow.push_back( bcMarker );
    }

    // Particle injectors
    for( std::string const& injectorMarker: this->modelProperties().boundaryConditions().markers( this->prefix(), "ParticleInjector" ) )
    {
        auto particleInjector = std::make_shared<levelsetparticleinjector_type>( this->shared_from_this(), markedelements( this->mesh(), injectorMarker) );

        // Injector method
        std::string particleInjectorMethod = "fixed_position"; // default value
        std::pair<bool, std::string> particleInjectorMethodRead = this->modelProperties().boundaryConditions().sparam( this->prefix(), "ParticleInjector", injectorMarker, "method" );
        if( particleInjectorMethodRead.first )
        {
            particleInjectorMethod = particleInjectorMethodRead.second;
        }
        particleInjector->setParticleInjectionMethod( particleInjectorMethod );

        // Injector particles
        auto const& injectorPTree = this->modelProperties().pTree().get_child( pt::ptree::path_type( "BoundaryConditions/"+this->prefix()+"/ParticleInjector/"+injectorMarker, '/' ) );
        if( auto const& injectorParticlesPTree = injectorPTree.get_child_optional( "particles" ) )
        {
            // read all particles with form "id": { "shape":[shape], [params]... }
            for( auto const& part: *injectorParticlesPTree )
            {
                std::string partId = part.first;
                std::string partShape = part.second.template get<std::string>( "shape" );
                parameter_map partParams = particleInjector->levelsetParticleShapes()->readShapeParams( partShape, part.second );
                particleInjector->addParticle( partId, partShape, partParams );
            }
        }

        // Injector particle number
        int nParticles = 1; // default
        std::pair<bool, int> nParticlesRead = this->modelProperties().boundaryConditions().iparam( this->prefix(), "ParticleInjector", injectorMarker, "nparticles" );
        if( nParticlesRead.first )
        {
            nParticles = nParticlesRead.second;
        }
        particleInjector->setNParticles( nParticles );

        M_levelsetParticleInjectors.push_back( particleInjector );
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Update levelset-dependent functions

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Markers accessors

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerCrossedElements() const
{
    if( !M_markerCrossedElements )
        M_markerCrossedElements.reset( new element_markers_type(this->functionSpaceMarkers(), "MarkerCrossedElements") );

    if( this->M_doUpdateMarkers )
       const_cast<self_type*>(this)->updateMarkerCrossedElements(); 

    return M_markerCrossedElements;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Advection
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletLinearPDE(
        sparse_matrix_ptrtype& A, vector_ptrtype& F ) const 
{
#if 0
    if( this->M_bcMarkersInflow.empty() ) return;

    this->log("LevelSet","updateBCStrongDirichletLinearPDE","start" );

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto const& phi = *this->phi();//this->fieldSolution();
    auto gradPhi = this->gradPhi();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for( auto const& bcMarker: this->M_bcMarkersInflow )
    {
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, bcMarker),
                    _element=phi,
                    _rhs=F,
                    _expr=(idv(this->timeStepBDF()->polyDeriv())-trans(idv(gradPhi))*idv(M_advectionToolbox->fieldAdvectionVelocity()))/(this->timeStepBDF()->polyDerivCoefficient(0))
              );
    }

    this->log("LevelSet","updateBCStrongDirichletLinearPDE","finish" );
#endif
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("LevelSet", "solve", "start");
    this->timerTool("Solve").start();

#if 0
    if( this->M_useExtensionVelocity )
    {
        auto uExt = this->extensionVelocity( idv(M_advectionToolbox->fieldAdvectionVelocity()) );
        this->updateAdvectionVelocity( uExt );
    }
#endif

    // Solve phi
    M_advectionToolbox->solve();
    this->setPhi( M_advectionToolbox->fieldUnknownPtr(), false );

#if 0
    if( M_useGradientAugmented )
    {
        // Update modGraPhi advection and reaction coefficients
        for( std::string const& matName : M_modGradPhiAdvection->materialsProperties()->physicToMaterials( M_modGradPhiAdvection->physicDefault() ) )
        {
            auto const& advectionVelocityExpr = M_advectionToolbox->materialsProperties()->materialProperty( matName, //TODO );
        }
        // Solve modGradPhi
        auto u = M_advectionToolbox->fieldAdvectionVelocityPtr();
        auto NxN = idv(this->N()) * trans(idv(this->N()));
        auto Du = sym( gradv(u) );
        M_modGradPhiAdvection->updateAdvectionVelocity( idv(u) );
        //M_modGradPhiAdvection->updateReactionCoeff( inner(NxN, Du) );
        //auto NxNDu = this->smoother()->project( inner(NxN, Du) );
        auto NxNDu = this->projectorL2()->project( inner(NxN, Du) );
        M_modGradPhiAdvection->updateReactionCoeff( idv(NxNDu) );
        //M_modGradPhiAdvection->updateSourceAdded(
                //- idv(modGradPhi) * inner( NxN, Du)
                //);
        M_modGradPhiAdvection->solve();
    }
    if( M_useStretchAugmented )
    {
        // Solve stretch modGradPhi
        //auto modGradPhi = M_stretchAdvection->fieldSolutionPtr();
        auto u = M_advectionToolbox->fieldAdvectionVelocityPtr();
        auto NxN = idv(this->N()) * trans(idv(this->N()));
        auto Du = sym( gradv(u) );
        auto NxNDu = this->projectorL2()->project( inner(NxN, Du) );
        M_stretchAdvection->updateAdvectionVelocity( idv(u) );
        //M_stretchAdvection->updateReactionCoeff( inner(NxN, Du) );
        M_stretchAdvection->updateReactionCoeff( idv(NxNDu) );
        //auto NxNDu = this->smoother()->project( inner(NxN, Du) );
        //auto NxNDu = vf::project( 
                //_space=M_stretchAdvection->functionSpace(),
                //_range=elements(M_stretchAdvection->mesh()),
                //_expr=inner(NxN, Du) 
                //);
        //M_stretchAdvection->updateReactionCoeff( idv(NxNDu) );
        //M_stretchAdvection->updateSourceAdded(
                //- idv(modGradPhi) * inner( NxN, Du)
                //);
        M_stretchAdvection->solve();
    }
    if( M_useCauchyAugmented )
    {
        auto u = M_advectionToolbox->fieldAdvectionVelocityPtr();
        M_backwardCharacteristicsAdvection->updateAdvectionVelocity( idv(u) );
        M_backwardCharacteristicsAdvection->solve();
    }
#endif

    // Update interface-related quantities
    this->updateInterfaceQuantities();

    // Correct volume if requested
    if( this->M_fixVolume )
    {
        auto const& phi = this->phi();
        double lambda = ( this->volume() - this->initialVolume() ) / this->perimeter();
        phi->add( lambda );
        // Request update interface-related quantities again since phi has changed
        // Note that updateInterfaceQuantities has lazy evaluation
        this->updateInterfaceQuantities();
    }
    // Correct area if requested
    if( this->M_fixArea )
    {
        auto const& phi = this->phi();
        auto const& K = this->K();
        auto const& D = this->D();
        double L = this->perimeter();
        double L0 = this->initialPerimeter();
        double Kbar = integrate( _range=this->rangeMeshElements(), _expr=idv(K)*idv(D)*idv(D) ).evaluate()(0,0) 
            / integrate( _range=this->rangeMeshElements(), _expr=idv(D)*idv(D) ).evaluate()(0,0);
        double mu = (L0-L) / integrate( _range=this->rangeMeshElements(), _expr=idv(K)*(idv(K)-Kbar)*idv(D)*idv(D) ).evaluate()(0,0);
        *phi = vf::project(
            _space=this->functionSpace(),
            _range=this->rangeMeshElements(),
            _expr=idv(phi) + mu*(idv(K)-Kbar)*idv(D)
            );
        // Request update interface-related quantities again since phi has changed
        // Note that updateInterfaceQuantities has lazy evaluation
        this->updateInterfaceQuantities();
    }

    if( !M_levelsetParticleInjectors.empty() )
    {
        auto const& phi = this->phi();
        // Inject particles
        for( auto const& particleInjector: M_levelsetParticleInjectors )
        {
            *phi = particleInjector->inject( *phi );
            this->updateInterfaceQuantities();
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("LevelSet","solve","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    M_advectionToolbox->startTimeStep();
    this->updateTime( M_advectionToolbox->time() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    double current_time = this->timeStepBDF()->time();

    if( this->hasRedistanciated() )
        M_iterSinceRedistanciation = 0;
    else
        ++M_iterSinceRedistanciation;

    M_vecIterSinceRedistanciation.push_back( M_iterSinceRedistanciation );

    this->saveCurrent();

    // Up advection field in case phi has been changed
    *(M_advectionToolbox->fieldUnknownPtr()) = this->phiElt();
    // Update time step
    M_advectionToolbox->updateTimeStep();
    if( M_useGradientAugmented )
        M_modGradPhiAdvection->updateTimeStep();
    if( M_useStretchAugmented )
        M_stretchAdvection->updateTimeStep();
    if( M_useCauchyAugmented )
        M_backwardCharacteristicsAdvection->updateTimeStep();

    this->updateTime( M_advectionToolbox->currentTime() );

    if( this->useOrder1AfterRedist() && M_iterSinceRedistanciation < this->timeOrder() )
    {
        this->timeStepBDF()->setTimeOrder( M_iterSinceRedistanciation + 1 );
        if( M_useGradientAugmented )
            M_modGradPhiAdvection->timeStepBdfUnknown()->setTimeOrder( M_iterSinceRedistanciation + 1 );
        if( M_useStretchAugmented )
            M_stretchAdvection->timeStepBdfUnknown()->setTimeOrder( M_iterSinceRedistanciation + 1 );
        if( M_useCauchyAugmented )
            M_backwardCharacteristicsAdvection->timeStepBdfUnknown()->setTimeOrder( M_iterSinceRedistanciation + 1 );
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_tensor2symm_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::leftCauchyGreenTensor() const
{
    if( M_useCauchyAugmented )
    {
        if( M_doUpdateCauchyGreenTensor )
        {
            const_cast<self_type*>(this)->updateLeftCauchyGreenTensor();
        }
    }
    else
    {
        throw std::logic_error( this->prefix()+".use-cauchy-augmented option must be true to use Cauchy-Green tensor" );
    }
    
    return M_leftCauchyGreenTensor;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
auto
LEVELSET_CLASS_TEMPLATE_TYPE::leftCauchyGreenTensorExpr() const
{
    CHECK( this->M_useCauchyAugmented ) << this->prefix()+".use-cauchy-augmented option must be true to use Cauchy-Green tensor";

    auto Y = M_backwardCharacteristicsAdvection->fieldUnknownPtr();
    auto const& N = this->N();

    return Feel::FeelModels::leftCauchyGreenTensorExpr( *Y, *N );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateLeftCauchyGreenTensor()
{
    CHECK( this->M_useCauchyAugmented ) << this->prefix()+".use-cauchy-augmented option must be true to use Cauchy-Green tensor";

    this->log("LevelSet", "updateLeftCauchyGreenTensor", "start");
    this->timerTool("UpdateInterfaceData").start();

    M_leftCauchyGreenTensor->on(
            _range=this->rangeMeshElements(),
            _expr=this->leftCauchyGreenTensorExpr()
            );

    M_doUpdateCauchyGreenTensor = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSet", "updateLeftCauchyGreenTensor", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
auto
LEVELSET_CLASS_TEMPLATE_TYPE::cauchyGreenInvariant1Expr() const
{
    return Feel::FeelModels::cauchyGreenInvariant1Expr( this->leftCauchyGreenTensorExpr() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_cauchygreen_invariant_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::cauchyGreenInvariant1() const
{
    if( !M_useCauchyAugmented )
        throw std::logic_error( "use-cauchy-augmented option must be true to use Cauchy-Green invariants" );

    if( this->M_doUpdateCauchyGreenInvariant1 )
    {
        this->log("LevelSet", "cauchyGreenInvariant1", "start");
        this->timerTool("UpdateInterfaceData").start();
#if 0
        M_cauchyGreenInvariant1->zero();
        M_cauchyGreenInvariant1->on(
                _range=this->rangeDiracElements(),
                _expr=trace(idv(this->leftCauchyGreenTensor()))
                );
#elif 0
        *M_cauchyGreenInvariant1 = vf::project(
                _space=M_cauchyGreenInvariant1->functionSpace(),
                _expr=trace(idv(this->leftCauchyGreenTensor()))
                );
#elif 1 // New implementation sqrt(tr(cof A))
        auto A = idv(this->leftCauchyGreenTensor());
        M_cauchyGreenInvariant1->zero();
        M_cauchyGreenInvariant1->on(
                _range=this->rangeDiracElements(),
                _expr=this->cauchyGreenInvariant1Expr()
                );
#endif
        M_doUpdateCauchyGreenInvariant1 = false;

        double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
        this->log("LevelSet", "cauchyGreenInvariant1", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
    }

    return M_cauchyGreenInvariant1;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
auto
LEVELSET_CLASS_TEMPLATE_TYPE::cauchyGreenInvariant2Expr() const
{
    return 0.5*trace( this->leftCauchyGreenTensorExpr() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_cauchygreen_invariant_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::cauchyGreenInvariant2() const
{
    if( !M_useCauchyAugmented )
        throw std::logic_error( "use-cauchy-augmented option must be true to use Cauchy-Green invariants" );

    if( this->M_doUpdateCauchyGreenInvariant2 )
    {
        this->log("LevelSet", "cauchyGreenInvariant2", "start");
        this->timerTool("UpdateInterfaceData").start();
#if 0
        auto A = idv(this->leftCauchyGreenTensor());
        auto trA = trace(A);
        // 3D: TrCofA = 1/2 (Tr2(A)-TrA2)
        *M_cauchyGreenInvariant2 = vf::project(
                _space=M_cauchyGreenInvariant2->functionSpace(),
                _expr=0.5*( trA*trA-trace(A*A) )
                );
#elif 0
        auto const& N = this->N();
        auto const& K = this->M_leftCauchyGreenTensor_K;
        auto const& KN = this->M_leftCauchyGreenTensor_KN;
        M_cauchyGreenInvariant2->zero();
        M_cauchyGreenInvariant2->on(
                _range=this->rangeDiracElements(),
                _expr=det(idv(K))/(trans(idv(N))*idv(KN))
                );
//#elif 1 // New implementation TrA / (2 sqrt(cofA))
#elif 1 // Old implementation TrA/2
        auto A = idv(this->leftCauchyGreenTensor());
        auto trA = trace(A);
        M_cauchyGreenInvariant2->zero();
        M_cauchyGreenInvariant2->on(
                _range=this->rangeDiracElements(),
                _expr=this->cauchyGreenInvariant2Expr()
                );
#endif
        M_doUpdateCauchyGreenInvariant2 = false;
        double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
        this->log("LevelSet", "cauchyGreenInvariant2", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
    }

    return M_cauchyGreenInvariant2;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateInterfaceQuantities()
{
    super_type::updateInterfaceQuantities();
    M_doUpdateCauchyGreenTensor = true;
    M_doUpdateCauchyGreenInvariant1 = true;
    M_doUpdateCauchyGreenInvariant2 = true;
}

//----------------------------------------------------------------------------//
// Redistanciation
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::redistanciate()
{ 
    this->log("LevelSet", "redistanciate", "start");
    this->timerTool("Redist").start();

    super_type::redistanciate();

    if( M_useGradientAugmented && M_reinitGradientAugmented )
    {
        CHECK( false ) << "TODO: ensure correct BDF restart";
        //auto sol = M_modGradPhiAdvection->fieldSolutionPtr();
        //sol->setConstant(1.);
    }
    if( M_useStretchAugmented && M_reinitStretchAugmented )
    {
        CHECK( false ) << "TODO: ensure correct BDF restart";
        //auto R = this->interfaceRectangularFunction();
        //auto sol = M_stretchAdvection->fieldSolutionPtr();
        //*sol = vf::project(
                //_space=M_stretchAdvection->functionSpace(),
                //_range=elements(M_stretchAdvection->mesh()),
                //_expr = 1. + (idv(sol)-1.)*idv(R)
                //);
    }

    double timeElapsed = this->timerTool("Redist").stop();
    this->log("LevelSet","redistanciate","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//LEVELSET_CLASS_TEMPLATE_DECLARATIONS
//std::shared_ptr<std::ostringstream>
//LEVELSET_CLASS_TEMPLATE_TYPE::getInfo() const
//{
    //std::string reinitMethod;
    //std::string reinitmethod = soption( _name="reinit-method", _prefix=this->prefix() );
    //if( reinitmethod == "fm" )
    //{
        //reinitMethod = "Fast-Marching";
        //std::string fmInitMethod = super_type::FastMarchingInitializationMethodIdMap.right.at( this->fastMarchingInitializationMethod() );
        //reinitMethod += " (" + fmInitMethod + ")";
    //}
    //else if( reinitmethod == "hj" )
        //reinitMethod = "Hamilton-Jacobi";

    //std::string scalarSmootherParameters;
    //if ( this->projectorSMScalar() )
    //{
        //double scalarSmootherCoeff = this->smoother()->epsilon() * Order / this->mesh()->hAverage();
        //scalarSmootherParameters = "coeff (h*c/order) = " 
            //+ std::to_string(this->smoother()->epsilon())
            //+ " (" + std::to_string(this->mesh()->hAverage()) + " * " + std::to_string(scalarSmootherCoeff) + " / " + std::to_string(Order) + ")"
            //;
    //}
    //std::string vectorialSmootherParameters;

    //if ( this->projectorSMVectorial() )
    //{
        //double vectorialSmootherCoeff = this->smootherVectorial()->epsilon() * Order / this->mesh()->hAverage();
        //vectorialSmootherParameters = "coeff (h*c/order) = " 
            //+ std::to_string(this->smootherVectorial()->epsilon())
            //+ " (" + std::to_string(this->mesh()->hAverage()) + " * " + std::to_string(vectorialSmootherCoeff) + " / " + std::to_string(Order) + ")"
            //;
    //}

    //std::string restartMode = (this->doRestart())? "ON": "OFF";

    //std::string exporterType = this->exporter()->type();
    //std::string hovisuMode = "OFF";
    //int exporterFreq = this->M_exporter->freq();
    //std::string exportedFields;
    //if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::GradPhi ) )
        //exportedFields = (exportedFields.empty())? "GradPhi": exportedFields+" - GradPhi";
    //if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::ModGradPhi ) )
        //exportedFields = (exportedFields.empty())? "ModGradPhi": exportedFields+" - ModGradPhi";
    //if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::Distance ) )
        //exportedFields = (exportedFields.empty())? "Distance": exportedFields+" - Distance";
    //if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::DistanceNormal ) )
        //exportedFields = (exportedFields.empty())? "DistanceNormal": exportedFields+" - DistanceNormal";
    //if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::DistanceCurvature ) )
        //exportedFields = (exportedFields.empty())? "DistanceCurvature": exportedFields+" - DistanceCurvature";
    //if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::AdvectionVelocity ) )
        //exportedFields = (exportedFields.empty())? "AdvectionVelocity": exportedFields+" - AdvectionVelocity";
    //if ( this->M_useStretchAugmented )
        //exportedFields = (exportedFields.empty())? "Stretch": exportedFields+" - Stretch";

    //std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    //*_ostr << "\n||==============================================||"
           //<< "\n||---------------Info : LevelSet----------------||"
           //<< "\n||==============================================||"
           //<< "\n   Prefix          : " << this->prefix()
           //<< "\n   Root Repository : " << this->rootRepository()
           //<< "\n   Dim             : " << nDim
           //<< "\n   Order           : " << Order
           //<< "\n   Periodicity     : " << this->periodicity().isPeriodic()

           //<< "\n   Level Set Parameters"
           //<< "\n     -- thickness interface (use adaptive)  : " << this->thicknessInterface() << " (" << std::boolalpha << this->M_useAdaptiveThicknessInterface << ")"
           //<< "\n     -- use regular phi (phi / |grad(phi)|) : " << std::boolalpha << this->M_useRegularPhi
           //<< "\n     -- Heaviside/Dirac projection method   : " << hdProjectionMethod
           //<< "\n     -- reinit initial value                : " << std::boolalpha << this->reinitInitialValue()
           //<< "\n     -- gradphi projection                  : " << gradPhiMethod
           //<< "\n     -- modgradphi projection               : " << modGradPhiMethod
           //<< "\n     -- curvature projection                : " << curvatureMethod
           //<< "\n     -- use gradient augmented              : " << std::boolalpha << this->M_useGradientAugmented
           //<< "\n     -- use stretch augmented               : " << std::boolalpha << this->M_useStretchAugmented

           //<< "\n   Reinitialization Parameters"
           //<< "\n     -- reinitialization method         : " << reinitMethod;
    //if( this->M_useGradientAugmented )
    //*_ostr << "\n     -- reinitialize gradient augmented : " << std::boolalpha << this->M_reinitGradientAugmented;
    //if( this->M_useGradientAugmented )
    //*_ostr << "\n     -- reinitialize stretch augmented  : " << std::boolalpha << this->M_reinitStretchAugmented;

    //if( this->projectorSMScalar() || this->projectorSMVectorial() )
    //*_ostr << "\n   Smoothers Parameters";
    //if( this->projectorSMScalar() )
    //*_ostr << "\n     -- scalar smoother    : " << scalarSmootherParameters;
    //if( this->projectorSMVectorial() )
    //*_ostr << "\n     -- vectorial smoother : " << vectorialSmootherParameters;

    //*_ostr << "\n   Space Discretization";
    //if( this->hasGeoFile() )
    //*_ostr << "\n     -- geo file name   : " << this->geoFile();
    //*_ostr << "\n     -- mesh file name  : " << this->meshFile()
           //<< "\n     -- nb elt in mesh  : " << this->mesh()->numGlobalElements()//numElements()
         ////<< "\n     -- nb elt in mesh  : " << this->mesh()->numElements()
         ////<< "\n     -- nb face in mesh : " << this->mesh()->numFaces()
           //<< "\n     -- hMin            : " << this->mesh()->hMin()
           //<< "\n     -- hMax            : " << this->mesh()->hMax()
           //<< "\n     -- hAverage        : " << this->mesh()->hAverage()
           //<< "\n     -- geometry order  : " << nOrderGeo
           //<< "\n     -- level set order : " << Order
           //<< "\n     -- nb dof          : " << this->functionSpace()->nDof() << " (" << this->functionSpace()->nLocalDof() << ")"
           //<< "\n     -- stabilization   : " << advectionStabilization;

    //*_ostr << "\n   Time Discretization"
           //<< "\n     -- initial time : " << this->timeStepBase()->timeInitial()
           //<< "\n     -- final time   : " << this->timeStepBase()->timeFinal()
           //<< "\n     -- time step    : " << this->timeStepBase()->timeStep()
           //<< "\n     -- order        : " << this->timeOrder()
           //<< "\n     -- restart mode : " << restartMode
           //<< "\n     -- save on disk : " << std::boolalpha << this->timeStepBase()->saveInFile();
    //if ( this->timeStepBase()->saveFreq() )
    //*_ostr << "\n     -- freq save : " << this->timeStepBase()->saveFreq()
           //<< "\n     -- file format save : " << this->timeStepBase()->fileFormat();

    //*_ostr << "\n   Exporter"
           //<< "\n     -- type            : " << exporterType
           //<< "\n     -- high order visu : " << hovisuMode
           //<< "\n     -- freq save       : " << exporterFreq
           //<< "\n     -- fields exported : " << exportedFields

           //<< "\n   Processors"
           //<< "\n     -- number of proc environment : " << Environment::worldComm().globalSize()
           //<< "\n     -- environment rank           : " << Environment::worldComm().rank()
           //<< "\n     -- global rank                : " << this->worldComm().globalRank()
           //<< "\n     -- local rank                 : " << this->worldComm().localRank()

           //<< "\n   Numerical Solver"
           //<< "\n     -- solver : " << M_advectionToolbox->solverName();

//#if 0
    //if ( this->algebraicFactory() )
    //*_ostr << this->algebraicFactory()->getInfo()->str();
//#endif
    ////if (enable_reinit)
    ////{
        ////if (reinitmethod == "hj")
        ////{
            ////infos << "\n      * hj maximum iteration per reinit : " << hj_max_iter
                  ////<< "\n      * hj pseudo time step dtau : " << hj_dtau
                  ////<< "\n      * hj stabilization : SUPG"
                  ////<< "\n      * hj coeff stab : " << option( prefixvm(M_prefix,"hj-coeff-stab")).template as<double>()
                  ////<< "\n      * hj tolerence on dist to dist error : "<<hj_tol;
        ////}
        ////else
        ////{
            ////infos << "\n      * fm smoothing coefficient for ILP : " << Environment::vm()[prefixvm(M_prefix,"fm-smooth-coeff")].template as< double >();
        ////}
    ////}
    ////infos << "\n\n  Level set spaces :"
          ////<< "\n     -- scalar LS space ndof : "<< this->functionSpace()->nDof()
          ////<< "\n     -- vectorial LS ndof : "<< this->functionSpaceVectorial()->nDof()
          ////<< "\n     -- scalar P0 space ndof : "<< this->functionSpaceMarkers()->nDof()
          ////<<"\n||==============================================||\n\n";

    //*_ostr << "\n||==============================================||"
           //<< "\n||==============================================||"
           //<< "\n\n";

    //return _ostr;
//}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
LEVELSET_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr = super_type::getInfo();

    std::string restartMode = (this->doRestart())? "ON": "OFF";
    std::string advectionStabilization = soption( _name="stabilization.method", _prefix=this->prefix() );

    *_ostr << "\n||==============================================||"
           << "\n||----------Info : LevelSet Advection-----------||"
           << "\n||==============================================||"
           << "\n   Advection"
           << "\n     -- stabilization   : " << advectionStabilization
           << "\n   Time Discretization"
           << "\n     -- initial time : " << this->timeStepBase()->timeInitial()
           << "\n     -- final time   : " << this->timeStepBase()->timeFinal()
           << "\n     -- time step    : " << this->timeStepBase()->timeStep()
           << "\n     -- order        : " << this->timeOrder()
           << "\n     -- restart mode : " << restartMode
           << "\n     -- save on disk : " << std::boolalpha << this->timeStepBase()->saveInFile();
    if ( this->timeStepBase()->saveFreq() )
    *_ostr << "\n     -- freq save : " << this->timeStepBase()->saveFreq()
           << "\n     -- file format save : " << this->timeStepBase()->fileFormat();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}

//----------------------------------------------------------------------------//
// Export results
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
std::set<std::string>
LEVELSET_CLASS_TEMPLATE_TYPE::postProcessSaveAllFieldsAvailable() const
{
    std::set<std::string> postProcessSaveAllFieldsAvailable = super_type::postProcessSaveAllFieldsAvailable();
    postProcessSaveAllFieldsAvailable.insert( { "advection-velocity", "stretch", "backwardcharacteristics", "cauchygreeninvariant1", "cauchygreeninvariant2" } );
    return postProcessSaveAllFieldsAvailable;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
std::set<std::string>
LEVELSET_CLASS_TEMPLATE_TYPE::postProcessExportsAllFieldsAvailable() const
{
    return super_type::postProcessExportsAllFieldsAvailable();
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    auto mfields = this->modelFields();
    //this->exportResults( time, this->symbolsExpr() );
    this->exportResults( time, this->symbolsExpr(mfields), mfields, this->allMeasuresQuantities() );
}

//----------------------------------------------------------------------------//
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::saveCurrent() const
{
    bool doSave = (this->timeStepBDF()->iteration() % this->timeStepBDF()->saveFreq() == 0);

    if (!doSave) return;

    if( this->worldComm().isMasterRank() )
    {
        fs::ofstream ofs( fs::path(this->rootRepository()) / fs::path( prefixvm(this->prefix(), "itersincereinit") ) );

        boost::archive::text_oarchive oa( ofs );
        oa << BOOST_SERIALIZATION_NVP( M_vecIterSinceRedistanciation );
    }
    this->worldComm().barrier();
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Update markers
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMarkerCrossedElements()
{
    // return a "marker" on the elements traversed by the interface between phio and phi
    // ie : mark as 1 is one of the dof of the element has  (phi * phio < 0)
    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();

    auto rangeElts = mesh->elementsWithProcessId( mesh->worldComm().localRank() );
    auto it_elt = std::get<0>( rangeElts );
    auto en_elt = std::get<1>( rangeElts );
    if (it_elt == en_elt) return;

    auto phi = this->phiPtr();
    auto phio = this->phiPreviousTimeStepPtr();

    auto prod = vf::project(_space=this->functionSpace(),
                            _range=this->rangeMeshElements(),
                            _expr=idv(phio) * idv(phi) );


    for (; it_elt!=en_elt; it_elt++)
    {
        auto const& elt = boost::unwrap_ref( *it_elt );
        bool mark_elt = false;
        for (int j=0; j<ndofv ; j++)
        {
            if (prod.localToGlobal(elt.id(), j, 0) <= 0.)
            {
                mark_elt = true;
                break;
            }
        }
        if( mark_elt )
            M_markerCrossedElements->assign(elt.id(), 0, 0, 1);
        else
            M_markerCrossedElements->assign(elt.id(), 0, 0, 0);
    }
}

} // FeelModels
} // Feel
