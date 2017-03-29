#include <feel/feelmodels/levelset/levelset.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>
#include <feel/feelmodels/levelset/reinitializer_hj.hpp>

namespace Feel {
namespace FeelModels {

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
std::map<std::string, typename LEVELSET_CLASS_TEMPLATE_TYPE::ShapeType>
LEVELSET_CLASS_TEMPLATE_TYPE::ShapeTypeMap = {
    {"sphere", ShapeType::SPHERE},
    {"ellipse", ShapeType::ELLIPSE}
};

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
LEVELSET_CLASS_TEMPLATE_TYPE::LevelSet( 
        std::string const& prefix,
        WorldComm const& worldComm,
        std::string const& subPrefix,
        std::string const& rootRepository ) 
:
    super_type( prefix, worldComm, subPrefix, rootRepository ),
    M_doUpdateDirac(true),
    M_doUpdateHeaviside(true),
    M_doUpdateNormal(true),
    M_doUpdateCurvature(true),
    M_doUpdateGradPhi(true),
    M_doUpdateModGradPhi(true),
    M_doUpdateMarkers(true),
    //M_periodicity(periodicityLS),
    M_reinitializerIsUpdatedForUse(false),
    M_hasReinitialized(false),
    M_hasReinitializedSmooth(false),
    M_iterSinceReinit(0)
{
    this->setFilenameSaveInfo( prefixvm(this->prefix(),"Levelset.info") );
    //-----------------------------------------------------------------------------//
    // Set advection model
    this->setModelName( "Advection" );
    //-----------------------------------------------------------------------------//
    // Load parameters
    this->loadParametersFromOptionsVm();
    // Load initial value
    this->loadConfigICFile();
    // Load boundary conditions
    this->loadConfigBCFile();
    // Load post-process
    this->loadConfigPostProcess();
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
        WorldComm const& worldComm,
        std::string const& subPrefix,
        std::string const& rootRepository )
{
    self_ptrtype new_ls( new self_type(prefix, worldComm, subPrefix, rootRepository) );
    return new_ls;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::build()
{
    this->log("LevelSet", "build", "start");

    super_type::build();
    this->createFunctionSpaces();
    this->createInterfaceQuantities();
    this->createReinitialization();
    this->createOthers();

    this->log("LevelSet", "build", "finish");
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::build( mesh_ptrtype const& mesh )
{
    this->log("LevelSet", "build (from mesh)", "start");

    super_type::build( mesh );
    this->createFunctionSpaces();
    this->createInterfaceQuantities();
    this->createReinitialization();
    this->createOthers();

    this->log("LevelSet", "build (from mesh)", "finish");
}

//----------------------------------------------------------------------------//
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::init()
{
    this->log("LevelSet", "init", "start");

    // Set levelset initial value
    this->initLevelsetValue();
    // Init levelset advection
    super_type::init( true, this->shared_from_this() );
    M_timeOrder = this->timeStepBDF()->timeOrder(); 

    // Init modGradPhi advection
    if( M_useGradientAugmented )
    {
        M_modGradPhiAdvection->init();
    }
    // Init stretch advection
    if( M_useStretchAugmented )
    {
        M_stretchAdvection->init();
    }

    // Init iterSinceReinit
    if( this->doRestart() )
    {
        // Reload saved iterSinceReinit data
        auto iterSinceReinitPath = fs::path(this->rootRepository()) / fs::path( prefixvm(this->prefix(), "itersincereinit") );
        if( fs::exists( iterSinceReinitPath ) )
        {
            fs::ifstream ifs( iterSinceReinitPath );
            boost::archive::text_iarchive ia( ifs );
            ia >> BOOST_SERIALIZATION_NVP( M_vecIterSinceReinit );
            M_iterSinceReinit = M_vecIterSinceReinit.back();
        }
        else
        {
            // If iterSinceReinit not found, we assume that last step reinitialized by default
            M_iterSinceReinit = 0;
        }
    }
    else
    {
            M_vecIterSinceReinit.push_back( M_iterSinceReinit );
    }
    // Adjust BDF order with iterSinceReinit
    if( M_iterSinceReinit < M_timeOrder )
    {
        this->timeStepBDF()->setTimeOrder( M_iterSinceReinit + 1 );
    }

    // Init post-process
    this->initPostProcess();

    this->log("LevelSet", "init", "finish");
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::initLevelsetValue()
{
    this->log("LevelSet", "initLevelsetValue", "start");

    bool hasInitialValue = false;

    auto phi_init = this->functionSpace()->element();
    phi_init.setConstant( std::numeric_limits<value_type>::max() );

    this->modelProperties().parameters().updateParameterValues();
    if( !this->M_icDirichlet.empty() )
    {
        M_icDirichlet.setParameterValues( this->modelProperties().parameters().toParameterValues() );

        for( auto const& iv : M_icDirichlet )
        {
            if( marker(iv).empty() )
            {
                phi_init = vf::project(
                        _space=this->functionSpace(),
                        _range=elements(this->mesh()),
                        _expr=expression(iv),
                        _geomap=this->geomap()
                        );
            }
            else
            {
                phi_init = vf::project(
                        _space=this->functionSpace(),
                        _range=markedelements(this->mesh(), marker(iv)),
                        _expr=expression(iv),
                        _geomap=this->geomap()
                        );
            }
        }

        hasInitialValue = true;
    }

    if( !this->M_icShapes.empty() )
    {
        // If phi_init already has a value, ensure that it is a proper distance function
        if( hasInitialValue )
        {
            // ILP on phi_init
            auto gradPhiInit = this->projectorL2Vectorial()->derivate( idv(phi_init) );
            auto modGradPhiInit = M_smootherFM->project( sqrt( trans(idv(gradPhiInit))*idv(gradPhiInit) ) );
            phi_init = vf::project( 
                _space=this->functionSpace(),
                _range=elements(this->mesh()),
                _expr=idv(phi_init) / idv(modGradPhiInit)
                );
            // Reinitialize phi_init
            phi_init = this->reinitializerFM()->run( phi_init );
        }
        // Add shapes
        for( auto const& shape: M_icShapes )
        {
            this->addShape( shape, phi_init );
        }

        hasInitialValue = true;
    }

    if( hasInitialValue )
    {
        this->setInitialValue( phi_init );
    }

    if( M_useGradientAugmented )
    {
        // Initialize modGradPhi
        if( M_reinitInitialValue )
        {
            M_modGradPhiAdvection->fieldSolutionPtr()->setConstant(1.);
        }
        else
        {
            auto gradPhi = this->gradPhi();
            *(M_modGradPhiAdvection->fieldSolutionPtr()) = vf::project( 
                    _space=this->functionSpace(),
                    _range=elements(this->mesh()),
                    _expr=sqrt( trans(idv(gradPhi))*idv(gradPhi) ) 
                    );
        }
    }
    if( M_useStretchAugmented )
    {
        // Initialize stretch modGradPhi
        M_stretchAdvection->fieldSolutionPtr()->setConstant(1.);
    }

    this->log("LevelSet", "initLevelsetValue", "finish");
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::addShape( 
        std::pair<ShapeType, parameter_map> const& shape, 
        element_levelset_type & phi 
        )
{
    ShapeType shapeType = shape.first;
    parameter_map const& shapeParams = shape.second;

    switch(shapeType)
    {
        case ShapeType::SPHERE:
        {
            auto X = Px() - shapeParams.dget("xc");
            auto Y = Py() - shapeParams.dget("yc");
            auto Z = Pz() - shapeParams.dget("zc"); 
            auto R = shapeParams.dget("radius");
            phi = vf::project(
                    _space=this->functionSpace(),
                    _range=elements(this->mesh()),
                    _expr=vf::min( idv(phi), sqrt(X*X+Y*Y+Z*Z)-R ),
                    _geomap=this->geomap()
                    );
        }
        break;

        case ShapeType::ELLIPSE:
        {
            // TODO
        }
        break;
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_levelset_type
LEVELSET_CLASS_TEMPLATE_TYPE::interfaceRectangularFunction( element_levelset_ptrtype const& p ) const
{
    //auto phi = idv(this->phi());
    auto phi = idv(p);
    double epsilon = M_thicknessInterfaceRectangularFunction;
    double epsilon_rect = 2.*epsilon;
    double epsilon_delta = (epsilon_rect - epsilon)/2.;
    double epsilon_zero = epsilon + epsilon_delta;

    auto R_expr =
        vf::chi( phi<-epsilon_rect )*vf::constant(0.0)
        +
        vf::chi( phi>=-epsilon_rect )*vf::chi( phi<=-epsilon )*
        1/2*(1 + (phi+epsilon_zero)/epsilon_delta + 1/M_PI*vf::sin( M_PI*(phi+epsilon_zero)/epsilon_delta ) )
        +
        vf::chi( phi>=-epsilon )*vf::chi( phi<=epsilon )*vf::constant(1.0)
        +
        vf::chi( phi>=epsilon )*vf::chi( phi<=epsilon_rect )*
        1/2*(1 - (phi-epsilon_zero)/epsilon_delta - 1/M_PI*vf::sin( M_PI*(phi-epsilon_zero)/epsilon_delta ) )
        +
        vf::chi(phi>epsilon_rect)*vf::constant(0.0)
        ;

    return vf::project( 
            this->functionSpace(), 
            elements(this->mesh()),
            R_expr
            );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    // Measures
    bool hasMeasureToExport = false;

    if( this->hasPostProcessMeasureExported( LevelSetMeasuresExported::Volume ) )
    {
        this->postProcessMeasuresIO().setMeasure( "volume", this->volume() );
        hasMeasureToExport = true;
    }
    if( this->hasPostProcessMeasureExported( LevelSetMeasuresExported::Perimeter ) )
    {
        this->postProcessMeasuresIO().setMeasure( "perimeter", this->perimeter() );
        hasMeasureToExport = true;
    }

    if ( hasMeasureToExport )
    {
        this->postProcessMeasuresIO().setParameter( "time", this->timeInitial() );
        // start or restart measure file
        if (!this->doRestart())
            this->postProcessMeasuresIO().start();
        else if ( !this->isStationary() )
            this->postProcessMeasuresIO().restart( "time", this->timeInitial() );
    }

    //// User-defined fields
    //if ( this->modelProperties().postProcess().find("Fields") != this->modelProperties().postProcess().end() )
    //{
        //for ( auto const& o : this->modelProperties().postProcess().find("Fields")->second )
        //{
            //if ( this->hasFieldUserScalar( o ) || this->hasFieldUserVectorial( o ) )
                //M_postProcessUserFieldExported.insert( o );
        //}
    //}
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::createFunctionSpaces()
{
    M_spaceLevelSetVec = space_levelset_vectorial_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
    M_spaceMarkers = space_markers_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::createInterfaceQuantities()
{
    if( Environment::vm().count( prefixvm(this->prefix(),"thickness-interface").c_str() ) )
        M_thicknessInterface = doption(prefixvm(this->prefix(),"thickness-interface"));
    else
        M_thicknessInterface = 1.5 * this->mesh()->hAverage();

    M_useAdaptiveThicknessInterface = boption(prefixvm(this->prefix(),"use-adaptive-thickness"));

    if( Environment::vm().count( prefixvm(this->prefix(),"thickness-interface-rectangular-function").c_str() ) )
        M_thicknessInterfaceRectangularFunction = doption(prefixvm(this->prefix(),"thickness-interface-rectangular-function")); 
    else
        M_thicknessInterfaceRectangularFunction = M_thicknessInterface;

    M_heaviside.reset( new element_levelset_type(this->functionSpace(), "Heaviside") );
    M_dirac.reset( new element_levelset_type(this->functionSpace(), "Dirac") );
    M_levelsetNormal.reset( new element_levelset_vectorial_type(this->functionSpaceVectorial(), "Normal") );
    M_levelsetCurvature.reset( new element_levelset_type(this->functionSpace(), "Curvature") );

    if( M_useGradientAugmented )
    {
        M_modGradPhiAdvection = modgradphi_advection_type::New(
                prefixvm(this->prefix(), "modgradphi-advection"),
                this->worldComm()
                );
        M_modGradPhiAdvection->setModelName( "Advection-Reaction" );
        M_modGradPhiAdvection->build( this->functionSpace() );
        M_modGradPhiAdvection->timeStepBDF()->setOrder( this->timeStepBDF()->bdfOrder() );

        M_modGradPhiAdvection->getExporter()->setDoExport( boption( _name="do_export_modgradphi-advection", _prefix=this->prefix() ) );
    }
    if( M_useStretchAugmented )
    {
        M_stretchAdvection = stretch_advection_type::New(
                prefixvm(this->prefix(), "stretch-advection"),
                this->worldComm()
                );
        M_stretchAdvection->setModelName( "Advection-Reaction" );
        M_stretchAdvection->build( this->functionSpace() );
        M_stretchAdvection->timeStepBDF()->setOrder( this->timeStepBDF()->bdfOrder() );

        M_stretchAdvection->getExporter()->setDoExport( boption( _name="do_export_stretch-advection", _prefix=this->prefix() ) );
    }
    for( std::string const& bcMarker: this->M_bcMarkersInflow )
    {
        if( M_useGradientAugmented )
            M_modGradPhiAdvection->addMarkerInflowBC( bcMarker );
        if( M_useStretchAugmented )
            M_stretchAdvection->addMarkerInflowBC( bcMarker );
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::createReinitialization()
{
    switch( M_reinitMethod )
    { 
        case LevelSetReinitMethod::FM :
        {
            if( !M_reinitializer )
                M_reinitializer = this->reinitializerFM();

            if( M_strategyBeforeFM == ILP )
            {
                M_backend_smooth = backend(
                        _name=prefixvm(this->prefix(), "smoother-fm"),
                        _worldcomm=this->worldComm()
                        );
                M_smootherFM = projector(
                        this->functionSpace(),/*domainSpace*/
                        this->functionSpace(),/*imageSpace*/
                        M_backend_smooth,
                        DIFF,
                        doption(_name="fm-smooth-coeff", _prefix=this->prefix())
                        );
            }
        }
        break;
        case LevelSetReinitMethod::HJ :
        {
            if( !M_reinitializer )
                M_reinitializer = this->reinitializerHJ();
            
            double thickness_heaviside;
            if( Environment::vm( _name="thickness-heaviside", _prefix=prefixvm(this->prefix(), "reinit-hj")).defaulted() )
            {
                thickness_heaviside =  M_thicknessInterface;
            }
            else
            {
                thickness_heaviside =  doption( _name="thickness-heaviside", _prefix=prefixvm(this->prefix(), "reinit-hj") );
            }
            boost::dynamic_pointer_cast<reinitializerHJ_type>(M_reinitializer)->setThicknessHeaviside( M_thicknessInterface );
        }
        break;
    }

    M_reinitializerIsUpdatedForUse = true;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::createOthers()
{
    M_projectorL2 = projector(this->functionSpace(), this->functionSpace(), backend(_name=prefixvm(this->prefix(), "projector-l2"), _worldcomm=this->worldComm()) );
    M_projectorL2Vec = projector(this->functionSpaceVectorial(), this->functionSpaceVectorial(), backend(_name=prefixvm(this->prefix(), "projector-l2-vec"), _worldcomm=this->worldComm()) );

    //if( M_doSmoothCurvature )
    //{
        //M_smoother = projector( 
                //this->functionSpace() , this->functionSpace(), 
                //backend(_name=prefixvm(this->prefix(),"smoother")), 
                //DIFF, 
                //this->mesh()->hAverage()*doption(_name="smooth-coeff", _prefix=this->prefix())/Order,
                //30);
    //}
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_levelset_vectorial_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::gradPhi() const
{
    if( !M_levelsetGradPhi )
        M_levelsetGradPhi.reset( new element_levelset_vectorial_type(this->functionSpaceVectorial(), "GradPhi") );

    if( M_doUpdateGradPhi )
       const_cast<self_type*>(this)->updateGradPhi(); 

    return M_levelsetGradPhi;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::modGradPhi() const
{
    if( M_useGradientAugmented )
    {
        return M_modGradPhiAdvection->fieldSolutionPtr();
    }
    else
    {
        if( !M_levelsetModGradPhi )
            M_levelsetModGradPhi.reset( new element_levelset_type(this->functionSpace(), "ModGradPhi") );

        if( M_doUpdateModGradPhi )
            const_cast<self_type*>(this)->updateModGradPhi();

        return M_levelsetModGradPhi;
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_stretch_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::stretch() const
{
    if( !M_levelsetStretch )
        M_levelsetStretch.reset( new element_stretch_type(this->M_stretchAdvection->functionSpace(), "Stretch") );

    if( M_useStretchAugmented )
    {
        *M_levelsetStretch = *(M_stretchAdvection->fieldSolutionPtr());
    }
    else
    {
        *M_levelsetStretch = *(this->modGradPhi());
    }

    M_levelsetStretch->add(-1.);

    return M_levelsetStretch;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::heaviside() const
{
    if( !M_heaviside )
        M_heaviside.reset( new element_levelset_type(this->functionSpace(), "Heaviside") );

    if( M_doUpdateHeaviside )
       const_cast<self_type*>(this)->updateHeaviside();

    return M_heaviside;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::dirac() const
{
    if( !M_dirac )
        M_dirac.reset( new element_levelset_type(this->functionSpace(), "Dirac") );

    if( M_doUpdateDirac )
       const_cast<self_type*>(this)->updateDirac();

    return M_dirac;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_levelset_vectorial_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::normal() const
{
    if( !M_levelsetNormal )
        M_levelsetNormal.reset( new element_levelset_vectorial_type(this->functionSpaceVectorial(), "Normal") );

    if( M_doUpdateNormal )
       const_cast<self_type*>(this)->updateNormal(); 

    return M_levelsetNormal;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::curvature() const
{
    if( !M_levelsetCurvature )
        M_levelsetCurvature.reset( new element_levelset_type(this->functionSpace(), "Curvature") );

    if( M_doUpdateCurvature )
       const_cast<self_type*>(this)->updateCurvature();

    return M_levelsetCurvature;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::setStrategyBeforeFm( int strat )
{
    if (M_reinitializerIsUpdatedForUse)
        LOG(INFO)<<" !!!  WARNING !!! : setStrategyBeforeFm set after the fast marching has been actually initialized ! \n";
    M_strategyBeforeFM = (strategy_before_FM_type) strat;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    super_type::loadParametersFromOptionsVm();

    M_useRegularPhi = boption(_name=prefixvm(this->prefix(),"use-regularized-phi"));
    M_useHeavisideDiracNodalProj = boption(_name=prefixvm(this->prefix(),"h-d-nodal-proj"));

    std::string reinitmethod = soption( _name="reinit-method", _prefix=this->prefix() );
    if( reinitmethod == "fm" )
        M_reinitMethod = LevelSetReinitMethod::FM;
    else if( reinitmethod == "hj" )
        M_reinitMethod = LevelSetReinitMethod::HJ;
    else
        CHECK( false ) << reinitmethod << " is not a valid reinitialization method\n";

    //M_useSmoothReinitialization = boption( _name="use-smooth-reinit", _prefix=this->prefix() );

    M_useMarkerDiracAsMarkerDoneFM = boption( _name="use-marker2-as-done", _prefix=prefixvm(this->prefix(), "reinit-fm") );

    M_strategyBeforeFM = (strategy_before_FM_type) ioption(prefixvm(this->prefix(),"fm-init-first-elts-strategy"));

    M_reinitInitialValue = boption( _name="reinit-initial-value", _prefix=this->prefix() );

    M_doSmoothGradient = boption( _name="smooth-gradient", _prefix=this->prefix() );

    if( Environment::vm( _name="smooth-curvature", _prefix=this->prefix()).defaulted() && Order < 2 )
        M_doSmoothCurvature = true;
    else
        M_doSmoothCurvature = boption( _name="smooth-curvature", _prefix=this->prefix() );

    M_useGradientAugmented = boption( _name="use-gradient-augmented", _prefix=this->prefix() );
    M_reinitGradientAugmented = boption( _name="reinit-gradient-augmented", _prefix=this->prefix() );

    M_useStretchAugmented = boption( _name="use-stretch-augmented", _prefix=this->prefix() );
    M_reinitStretchAugmented = boption( _name="reinit-stretch-augmented", _prefix=this->prefix() );

    //M_doExportAdvection = boption(_name="export-advection", _prefix=this->prefix());
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::loadConfigICFile()
{
    auto const& initialConditions = this->modelProperties().initialConditions();

    this->M_icDirichlet = initialConditions.getScalarFields( std::string(this->prefix()), "Dirichlet" );
    
    // Shapes
    for( std::string const& icShape: initialConditions.markers( this->prefix(), "shapes") )
    {
        parameter_map shapeParameterMap;

        auto shapeTypeRead = initialConditions.sparam( this->prefix(), "shapes", icShape, "shape" );
        auto shapeTypeIt = ShapeTypeMap.find(shapeTypeRead.second);
        if( shapeTypeIt != ShapeTypeMap.end() )
        {
            switch(shapeTypeIt->second)
            {
                case ShapeType::SPHERE:
                {
                    auto xcRead = initialConditions.dparam( this->prefix(), "shapes", icShape, "xc" );
                    CHECK(xcRead.first) << icShape << " xc not provided\n";
                    auto ycRead = initialConditions.dparam( this->prefix(), "shapes", icShape, "yc" );
                    CHECK(ycRead.first) << icShape << " yc not provided\n";
                    auto zcRead = initialConditions.dparam( this->prefix(), "shapes", icShape, "zc" );
                    CHECK(zcRead.first || nDim < 3) << icShape << " zc not provided\n";
                    auto radiusRead = initialConditions.dparam( this->prefix(), "shapes", icShape, "radius" );
                    CHECK(radiusRead.first) << icShape << " radius not provided\n";

                    shapeParameterMap["id"] = icShape;
                    shapeParameterMap["xc"] = xcRead.second;
                    shapeParameterMap["yc"] = ycRead.second;
                    shapeParameterMap["zc"] = zcRead.second;
                    shapeParameterMap["radius"] = radiusRead.second;
                }
                break;

                case ShapeType::ELLIPSE:
                {
                    auto xcRead = initialConditions.dparam( this->prefix(), "shapes", icShape, "xc" );
                    CHECK(xcRead.first) << icShape << " xc not provided\n";
                    auto ycRead = initialConditions.dparam( this->prefix(), "shapes", icShape, "yc" );
                    CHECK(ycRead.first) << icShape << " yc not provided\n";
                    auto zcRead = initialConditions.dparam( this->prefix(), "shapes", icShape, "zc" );
                    CHECK(zcRead.first || nDim < 3) << icShape << " zc not provided\n";
                    auto radiusRead = initialConditions.dparam( this->prefix(), "shapes", icShape, "radius" );
                    CHECK(radiusRead.first) << icShape << " radius not provided\n";

                    shapeParameterMap["id"] = icShape;
                    shapeParameterMap["xc"] = xcRead.second;
                    shapeParameterMap["yc"] = ycRead.second;
                    shapeParameterMap["zc"] = zcRead.second;
                    shapeParameterMap["radius"] = radiusRead.second;
                    // TODO
                }
                break;
            }

            M_icShapes.push_back( std::make_pair(shapeTypeIt->second, shapeParameterMap) );
        }
        else
        {
            CHECK(false) << "invalid shape type in " << icShape << std::endl;
        }
    } 
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    M_bcMarkersInflow.clear();

    for( std::string const& bcMarker: this->modelProperties().boundaryConditions().markers( this->prefix(), "inflow" ) )
    {
        if( std::find(M_bcMarkersInflow.begin(), M_bcMarkersInflow.end(), bcMarker) == M_bcMarkersInflow.end() )
            M_bcMarkersInflow.push_back( bcMarker );
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::loadConfigPostProcess()
{
    auto& modelPostProcess = this->modelProperties().postProcess();

    auto physicalQuantities = modelPostProcess.pTree().get_child_optional("PhysicalQuantities");
    if( physicalQuantities )
    {
        for( auto& i: *physicalQuantities )
        {
            auto i_str = i.second.template get_value<std::string>();
            modelPostProcess["PhysicalQuantities"].push_back( i_str  );
            LOG(INFO) << "add to postprocess physical quantity " << i_str;
        }
    }

    if ( modelPostProcess.find("PhysicalQuantities") != modelPostProcess.end() )
    {
        for ( auto const& o : modelPostProcess.find("PhysicalQuantities")->second )
        {
            if( o == "volume" || o == "all" )
                this->M_postProcessMeasuresExported.insert( LevelSetMeasuresExported::Volume );
            if( o == "perimeter" || o == "all" )
                this->M_postProcessMeasuresExported.insert( LevelSetMeasuresExported::Perimeter );
        }
    }

    // Load Fields from JSON
    if ( modelPostProcess.find("Fields") != modelPostProcess.end() )
    {
        for( auto const& o : modelPostProcess.find("Fields")->second )
        {
            if( o == "gradphi" || o == "all" )
                this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::GradPhi );
            if( o == "modgradphi" || o == "all" )
                this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::ModGradPhi );
        }
    }
    // Overwrite with options from CFG
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_gradphi").c_str()) )
        if ( boption(_name="do_export_gradphi",_prefix=this->prefix()) )
            this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::GradPhi );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_modgradphi").c_str()) )
        if ( boption(_name="do_export_modgradphi",_prefix=this->prefix()) )
            this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::ModGradPhi );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::reinitializer_ptrtype 
LEVELSET_CLASS_TEMPLATE_TYPE::buildReinitializer( 
        LevelSetReinitMethod method, 
        space_levelset_ptrtype const& space,
        std::string const& prefix )
{
    switch( method )
    { 
        case LevelSetReinitMethod::FM :
        {
            return reinitializer_ptrtype(
                    new ReinitializerFM<space_levelset_type>( space, prefixvm(prefix, "reinit-fm") ) 
                    );
        }
        break;
        case LevelSetReinitMethod::HJ :
        {
            return reinitializer_ptrtype(
                    new ReinitializerHJ<space_levelset_type>( space, prefixvm(prefix, "reinit-hj") )
                    );
        }
        break;
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Update levelset-dependent functions
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateGradPhi()
{
    auto phi = this->phi();
    if( M_doSmoothGradient )
    {
        //*M_levelsetGradPhi = this->smootherVectorial()->derivate( idv(phi) );
        *M_levelsetGradPhi = this->smootherVectorial()->project( trans(gradv(phi)) );
    }
    else
    {
        //*M_levelsetGradPhi = this->projectorL2Vectorial()->derivate( idv(phi) );
        *M_levelsetGradPhi = this->projectorL2Vectorial()->project( _expr=trans(gradv(phi)) );
        //*M_levelsetGradPhi = vf::project( 
                //_space=this->functionSpaceVectorial(),
                //_range=elements(this->mesh()),
                //_expr=trans(gradv(phi)) 
                //);
    }

    M_doUpdateGradPhi = false;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateModGradPhi()
{
    auto gradPhi = this->gradPhi();
    *M_levelsetModGradPhi = vf::project( 
            _space=this->functionSpace(),
            _range=elements(this->mesh()),
            _expr=sqrt( trans(idv(gradPhi))*idv(gradPhi) ) 
            );
    //*M_levelsetModGradPhi = this->projectorL2()->project( _expr=sqrt( trans(idv(gradPhi))*idv(gradPhi) ) );

    M_doUpdateModGradPhi = false;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateDirac()
{
    this->log("LevelSet", "updateDirac", "start");

    // derivative of Heaviside function
    auto eps0 = this->thicknessInterface();
    auto eps_elt = this->functionSpace()->element();

    if( M_useAdaptiveThicknessInterface )
    {
        auto gradPhi = this->gradPhi();
        auto gradPhiX = vf::project(
                _space=this->functionSpace(),
                _range=elements(this->mesh()),
                _expr=idv(gradPhi->comp(Component::X))
                );
        auto gradPhiY = vf::project(
                _space=this->functionSpace(),
                _range=elements(this->mesh()),
                _expr=idv(gradPhi->comp(Component::Y))
                );
#if FEELPP_DIM == 3
        auto gradPhiZ = vf::project(
                _space=this->functionSpace(),
                _range=elements(this->mesh()),
                _expr=idv(gradPhi->comp(Component::Z))
                );
#endif
        eps_elt = vf::project(
                _space=this->functionSpace(),
                _range=elements(this->mesh()),
                _expr=(vf::abs(idv(gradPhiX))+vf::abs(idv(gradPhiY))
#if FEELPP_DIM == 3
                    + vf::abs(idv(gradPhiZ))
#endif
                    )*cst(eps0)/idv(this->modGradPhi())
                );
    }
    else
    {
        eps_elt.setConstant(eps0); 
    }

    auto eps = idv(eps_elt);

    if (M_useRegularPhi)
    {
        //auto psi = idv(this->phi()) / sqrt( gradv(this->phi()) * trans(gradv(this->phi())) );
        auto psi = idv(this->phi()) / idv(this->modGradPhi());
        auto D_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/(2*eps) *( 1 + cos(M_PI*psi/eps) )
            +
            vf::chi(psi>eps)*vf::constant(0.0);

        if ( M_useHeavisideDiracNodalProj )
            *M_dirac = vf::project( this->functionSpace(), elements(this->mesh()), D_expr );
        else
            *M_dirac = M_projectorL2->project(D_expr);
    }
    else
    {
        auto psi = idv(this->phi()) ;
        auto D_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/(2*eps) *( 1 + cos(M_PI*psi/eps) )
            +
            vf::chi(psi>eps)*vf::constant(0.0);

        if (M_useHeavisideDiracNodalProj)
            *M_dirac = vf::project( this->functionSpace(), elements(this->mesh()), D_expr );
        else
            *M_dirac = M_projectorL2->project(D_expr);
    }

    M_doUpdateDirac = false;

    this->log("LevelSet", "updateDirac", "finish");
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateHeaviside()
{ 
    this->log("LevelSet", "updateHeaviside", "start");

    auto eps = this->thicknessInterface();

    if (M_useRegularPhi)
    {
        auto psi = idv(this->phi()) / vf::sqrt( gradv(this->phi()) * trans(gradv(this->phi())) );
        auto H_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/2*(1 + psi/eps + 1/M_PI*vf::sin( M_PI*psi/eps ) )
            +
            vf::chi(psi>eps)*vf::constant(1.0);

        if (M_useHeavisideDiracNodalProj)
            *M_heaviside = vf::project(this->functionSpace(), elements(this->mesh()), H_expr);
        else
            *M_heaviside = M_projectorL2->project(H_expr);
    }
    else
    {
        auto psi = idv(this->phi());
        auto H_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            1/2*(1 + psi/eps + 1/M_PI*vf::sin( M_PI*psi/eps ) )
            +
            vf::chi(psi>eps)*vf::constant(1.0);

        if (M_useHeavisideDiracNodalProj)
            *M_heaviside = vf::project(this->functionSpace(), elements(this->mesh()), H_expr);
        else
            *M_heaviside = M_projectorL2->project(H_expr);
    }

    M_doUpdateHeaviside = false;

    this->log("LevelSet", "updateHeaviside", "finish");
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateNormal()
{
    this->log("LevelSet", "updateNormal", "start");

    auto phi = this->phi();
    //*M_levelsetNormal = M_projectorL2Vec->project( _expr=trans(gradv(phi)) / sqrt(gradv(phi) * trans(gradv(phi))) );
    auto gradPhi = this->gradPhi();
    *M_levelsetNormal = vf::project( 
            _space=this->functionSpaceVectorial(),
            _range=elements(this->mesh()),
            //_expr=trans(gradv(phi)) / sqrt(gradv(phi) * trans(gradv(phi))) 
            _expr=idv(gradPhi) / sqrt(trans(idv(gradPhi)) * idv(gradPhi)) 
            );

    M_doUpdateNormal = false;

    this->log("LevelSet", "updateNormal", "finish");
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateCurvature()
{
    this->log("LevelSet", "updateCurvature", "start");

    if( M_doSmoothCurvature )
    {
        *M_levelsetCurvature = this->smoother()->project( _expr=divv(this->normal()) );
    }
    else
    {
        *M_levelsetCurvature = this->projectorL2()->project( _expr=divv(this->normal()) );
        //*M_levelsetCurvature = vf::project( 
                //_space=this->functionSpace(),
                //_range=elements(this->mesh()),
                //_expr=divv(this->normal())
                //);
    }

    M_doUpdateCurvature = false;

    this->log("LevelSet", "updateCurvature", "finish");
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::projector_levelset_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::smoother() const
{
    if( !M_smoother )
    {
        M_smoother = projector( 
                this->functionSpace() , this->functionSpace(), 
                backend(_name=prefixvm(this->prefix(),"smoother"), _worldcomm=this->worldComm()), 
                DIFF, 
                this->mesh()->hAverage()*doption(_name="smooth-coeff", _prefix=prefixvm(this->prefix(),"smoother"))/Order,
                30);
    }
    return M_smoother; 
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::projector_levelset_vectorial_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::smootherVectorial() const
{
    if( !M_smootherVectorial )
        M_smootherVectorial = projector( 
                this->functionSpaceVectorial() , this->functionSpaceVectorial(), 
                backend(_name=prefixvm(this->prefix(),"smoother-vec"), _worldcomm=this->worldComm()), 
                DIFF, 
                this->mesh()->hAverage()*doption(_name="smooth-coeff", _prefix=prefixvm(this->prefix(),"smoother-vec"))/Order,
                30);
    return M_smootherVectorial; 
}
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Markers accessors
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerInterface()
{
    if( !M_markerInterface )
        M_markerInterface.reset( new element_markers_type(M_spaceMarkers, "MarkerInterface") );

    if( M_doUpdateMarkers )
       this->updateMarkerInterface(); 

    return M_markerInterface;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerDirac()
{
    if( !M_markerDirac )
        M_markerDirac.reset( new element_markers_type(M_spaceMarkers, "MarkerDirac") );

    if( M_doUpdateMarkers )
       this->updateMarkerDirac(); 

    return M_markerDirac;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerHeaviside(bool invert, bool cut_at_half)
{
    if( !M_markerHeaviside )
        M_markerHeaviside.reset( new element_markers_type(M_spaceMarkers, "MarkerHeaviside") );

    if( M_doUpdateMarkers )
       this->updateMarkerHeaviside( invert, cut_at_half );

    return M_markerHeaviside;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerCrossedElements()
{
    if( !M_markerCrossedElements )
        M_markerCrossedElements.reset( new element_markers_type(M_spaceMarkers, "MarkerCrossedElements") );

    if( M_doUpdateMarkers )
       this->updateMarkerCrossedElements(); 

    return M_markerCrossedElements;
}

//----------------------------------------------------------------------------//
// Utility distances
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype
LEVELSET_CLASS_TEMPLATE_TYPE::distToBoundary()
{
    element_levelset_ptrtype distToBoundary( new element_levelset_type(this->functionSpace(), "DistToBoundary") );

    // Retrieve the elements touching the boundary
    auto boundaryelts = boundaryelements( this->mesh() );

    // Mark the elements in myrange
    auto mymark = vf::project(
            _space=this->functionSpaceMarkers(),
            _range=boundaryelts,
            _expr=cst(1)
            );

    // Update mesh marker2
    this->mesh()->updateMarker2( mymark );

    // Init phi0 with h on marked2 elements
    auto phi0 = vf::project(
            _space=this->functionSpace(),
            _range=boundaryelts,
            _expr=h()
            );
    phi0.on( _range=boundaryfaces(this->mesh()), _expr=cst(0.) );

    // Run FM using marker2 as marker DONE
    this->reinitializerFM()->setUseMarker2AsMarkerDone( true );
    *distToBoundary = this->reinitializerFM()->run( phi0 );

    return distToBoundary;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype
LEVELSET_CLASS_TEMPLATE_TYPE::distToMarkedFaces( boost::any const& marker )
{
    element_levelset_ptrtype distToMarkedFaces( new element_levelset_type(this->functionSpace(), "DistToMarkedFaces") );

    typedef boost::reference_wrapper<typename MeshTraits<mesh_type>::element_type const> element_ref_type;
    typedef std::vector<element_ref_type> cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    // Retrieve the elements touching the marked faces
    auto mfaces = markedfaces( this->mesh(), marker );
    for( auto const& face: mfaces )
    {
        if( face.isConnectedTo0() )
            myelts->push_back(boost::cref(face.element0()));
        if( face.isConnectedTo1() )
            myelts->push_back(boost::cref(face.element1()));
    }

    auto myrange = boost::make_tuple(
            mpl::size_t<MESH_ELEMENTS>(), myelts->begin(), myelts->end(), myelts
            );

    // Mark the elements in myrange
    auto mymark = this->functionSpaceMarkers()->element();
    mymark.on( _range=myrange, _expr=cst(1) );

    // Update mesh marker2
    this->mesh()->updateMarker2( mymark );

    // Init phi0 with h on marked2 elements
    auto phi0 = vf::project(
            _space=this->functionSpace(),
            _range=myrange,
            _expr=h()
            );
    phi0.on( _range=mfaces, _expr=cst(0.) );

    // Run FM using marker2 as marker DONE
    this->reinitializerFM()->setUseMarker2AsMarkerDone( true );
    *distToMarkedFaces = this->reinitializerFM()->run( phi0 );

    return distToMarkedFaces;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype
LEVELSET_CLASS_TEMPLATE_TYPE::distToMarkedFaces( std::initializer_list<boost::any> marker )
{
    element_levelset_ptrtype distToMarkedFaces( new element_levelset_type(this->functionSpace(), "DistToMarkedFaces") );

    typedef boost::reference_wrapper<typename MeshTraits<mesh_type>::element_type const> element_ref_type;
    typedef std::vector<element_ref_type> cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    // Retrieve the elements touching the marked faces
    auto mfaces_list = markedfaces( this->mesh(), marker );
    for( auto const& mfaces: mfaces_list )
    {
        for( auto const& face: mfaces )
        {
            if( face.isConnectedTo0() )
                myelts->push_back(boost::cref(face.element0()));
            if( face.isConnectedTo1() )
                myelts->push_back(boost::cref(face.element1()));
        }
    }

    auto myrange = boost::make_tuple(
            mpl::size_t<MESH_ELEMENTS>(), myelts->begin(), myelts->end(), myelts
            );

    // Mark the elements in myrange
    auto mymark = this->functionSpaceMarkers()->element();
    mymark.on( _range=myrange, _expr=cst(1) );

    // Update mesh marker2
    this->mesh()->updateMarker2( mymark );

    // Init phi0 with h on marked2 elements
    auto phi0 = vf::project(
            _space=this->functionSpace(),
            _range=myrange,
            _expr=h()
            );
    phi0.on( _range=mfaces_list, _expr=cst(0.) );

    // Run FM using marker2 as marker DONE
    this->reinitializerFM()->setUseMarker2AsMarkerDone( true );
    *distToMarkedFaces = this->reinitializerFM()->run( phi0 );

    return distToMarkedFaces;
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
    if( this->M_bcMarkersInflow.empty() ) return;

    this->log("LevelSet","updateBCStrongDirichletLinearPDE","start" );

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto const& phi = this->fieldSolution();
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
                    _expr=(idv(this->timeStepBDF()->polyDeriv())-trans(idv(gradPhi))*idv(this->fieldAdvectionVelocity()))/(this->timeStepBDF()->polyDerivCoefficient(0))
              );
    }

    this->log("LevelSet","updateBCStrongDirichletLinearPDE","finish" );

}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("LevelSet", "solve", "start");
    this->timerTool("Solve").start();

    // Solve phi
    super_type::solve();

    if( M_useGradientAugmented )
    {
        // Solve modGradPhi
        //auto modGradPhi = M_modGradPhiAdvection->fieldSolutionPtr();
        auto u = this->fieldAdvectionVelocityPtr();
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
        auto u = this->fieldAdvectionVelocityPtr();
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
    // Update interface-related quantities
    this->updateInterfaceQuantities();

    // Reset hasReinitialized
    M_hasReinitialized = false;
    M_hasReinitializedSmooth = false;

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("LevelSet","solve","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    double current_time = this->timeStepBDF()->time();

    if( M_hasReinitialized )
        M_iterSinceReinit = 0;
    else
        ++M_iterSinceReinit;

    M_vecIterSinceReinit.push_back( M_iterSinceReinit );

    this->saveCurrent();

    super_type::updateTimeStep();
    if( M_useGradientAugmented )
        M_modGradPhiAdvection->updateTimeStep();
    if( M_useStretchAugmented )
        M_stretchAdvection->updateTimeStep();

    if( M_iterSinceReinit < M_timeOrder )
    {
        this->timeStepBDF()->setTimeOrder( M_iterSinceReinit + 1 );
        if( M_useGradientAugmented )
            M_modGradPhiAdvection->timeStepBDF()->setTimeOrder( M_iterSinceReinit + 1 );
        if( M_useStretchAugmented )
            M_stretchAdvection->timeStepBDF()->setTimeOrder( M_iterSinceReinit + 1 );
        //this->timeStepBDF()->setTimeInitial( current_time ); 
        //this->timeStepBDF()->restart();
        //this->timeStepBDF()->setTimeInitial( this->timeInitial() );
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateInterfaceQuantities()
{
    M_doUpdateDirac = true;
    M_doUpdateHeaviside = true;
    M_doUpdateNormal = true;
    M_doUpdateCurvature = true;
    M_doUpdateMarkers = true;
    M_doUpdateGradPhi = true;
    M_doUpdateModGradPhi = true;
}

//----------------------------------------------------------------------------//
// Reinitialization
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::reinitialize( bool useSmoothReinit )
{ 
    this->log("LevelSet", "reinitialize", "start");
    this->timerTool("Reinit").start();

    if( !M_reinitializerIsUpdatedForUse )
        this->createReinitialization();

    auto phi = this->phi();
    auto phiReinit = this->functionSpace()->elementPtr();

    if ( M_reinitMethod == LevelSetReinitMethod::FM )
    {
        if ( M_useMarkerDiracAsMarkerDoneFM )
        {
            this->mesh()->updateMarker2( *this->markerDirac() );
        }

        switch (M_strategyBeforeFM)
        {
            case ILP :
            {
                // save the smoothed gradient magnitude of phi
                //auto modgradphi = M_smootherFM->project( vf::min(vf::max(vf::sqrt(inner(gradv(phi), gradv(phi))), 0.92), 2.) );
                //auto gradPhi = idv(this->gradPhi());
                auto gradPhi = trans(gradv(phi));
                auto modgradphi = M_smootherFM->project( sqrt( trans(gradPhi)*gradPhi ) );
                
                //*phi = vf::project(this->functionSpace(), elements(this->mesh()), idv(phi)/idv(modgradphi) );
                *phiReinit = vf::project(this->functionSpace(), elements(this->mesh()), idv(phi)/idv(modgradphi) );
            }
            break;

            case HJ_EQ :
            {
                CHECK(false) << "TODO\n";
                //*phi = *explicitHJ(max_iter, dtau, tol);
            }
            break;
            case NONE :
            {
                *phiReinit = *phi;
            }
            break;

            default:
            {
                CHECK(false)<<"no strategy chosen to initialize first elements before fast marching\n"
                            <<"please, consider setting the option fm-init-first-elts-strategy\n";
            }
            break;
        } // switch M_strategyBeforeFM

        // Fast Marching Method
        boost::dynamic_pointer_cast<reinitializerFM_type>( M_reinitializer )->setUseMarker2AsMarkerDone( 
                M_useMarkerDiracAsMarkerDoneFM 
                );

        LOG(INFO)<< "reinit with FMM done"<<std::endl;
    } // Fast Marching

    else if ( M_reinitMethod == LevelSetReinitMethod::HJ )
    {
        *phiReinit = *phi;
        // TODO
    } // Hamilton-Jacobi

    //*phi = M_reinitializer->run( *phi );
    if( useSmoothReinit )
    {
        *phiReinit = this->smoother()->project(
                _expr=idv(*phiReinit)
                );
    }

    *phi = M_reinitializer->run( *phiReinit );

    //*phiReinit = M_reinitializer->run( *phiReinit );
    //if( useSmoothReinit )
    //{
        ////auto R = this->interfaceRectangularFunction(phiReinit);
        //auto R = this->interfaceRectangularFunction();
        //*phi = vf::project(
                //_space=this->functionSpace(),
                //_range=elements(this->mesh()),
                //_expr=idv(phi)*idv(R) + idv(phiReinit)*(1.-idv(R))
                //);
    //}
    //else
    //{
        //*phi = *phiReinit;
    //}

    if( M_useGradientAugmented && M_reinitGradientAugmented )
    {
        auto sol = M_modGradPhiAdvection->fieldSolutionPtr();
        sol->setConstant(1.);
    }
    if( M_useStretchAugmented && M_reinitStretchAugmented )
    {
        auto R = this->interfaceRectangularFunction();
        auto sol = M_stretchAdvection->fieldSolutionPtr();
        *sol = vf::project(
                _space=M_stretchAdvection->functionSpace(),
                _range=elements(M_stretchAdvection->mesh()),
                _expr = 1. + (idv(sol)-1.)*idv(R)
                );
    }

    if( useSmoothReinit )
    {
        M_hasReinitializedSmooth = true;
        M_hasReinitialized = true;
    }
    else
        M_hasReinitialized = true;

    double timeElapsed = this->timerTool("Reinit").stop();
    this->log("LevelSet","reinitialize","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::reinitializerFM_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::reinitializerFM( bool buildOnTheFly )
{
    if( !M_reinitializerFM && buildOnTheFly )
    {
        M_reinitializerFM.reset( 
                new ReinitializerFM<space_levelset_type>( this->functionSpace(), prefixvm(this->prefix(), "reinit-fm") ) 
                );
    }

    return M_reinitializerFM;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::reinitializerHJ_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::reinitializerHJ( bool buildOnTheFly )
{
    if( !M_reinitializerHJ && buildOnTheFly )
    {
        M_reinitializerHJ.reset( 
                new ReinitializerHJ<space_levelset_type>( this->functionSpace(), prefixvm(this->prefix(), "reinit-hj") ) 
                );
    }

    return M_reinitializerHJ;
}

//----------------------------------------------------------------------------//
// Initial value
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::setInitialValue(element_levelset_type const& phiv, bool doReinitialize)
{
    this->log("LevelSet", "setInitialValue", "start");

    *this->phi() = phiv;

    if (doReinitialize)
    {
        this->reinitialize();
        // The initial reinitialization is not a "real" one
        M_hasReinitialized = false;
    }

    this->updateInterfaceQuantities();

    this->log("LevelSet", "setInitialValue", "finish");
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
boost::shared_ptr<std::ostringstream>
LEVELSET_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::string advectionStabilization = soption( _name="advec-stab-method", _prefix=this->prefix() );

    std::string hdProjectionMethod = (this->M_useHeavisideDiracNodalProj)? "nodal": "L2";

    std::string reinitMethod;
    std::string reinitmethod = soption( _name="reinit-method", _prefix=this->prefix() );
    if( reinitmethod == "fm" )
        reinitMethod = "Fast-Marching";
    else if( reinitmethod == "hj" )
        reinitMethod = "Hamilton-Jacobi";

    std::string scalarSmootherParameters;
    double scalarSmootherCoeff = this->smoother()->epsilon() * Order / this->mesh()->hAverage();
    scalarSmootherParameters = "coeff (h*c/order) = " 
        + std::to_string(this->smoother()->epsilon())
        + " (" + std::to_string(this->mesh()->hAverage()) + " * " + std::to_string(scalarSmootherCoeff) + " / " + std::to_string(Order) + ")"
        ;
    std::string vectorialSmootherParameters;
    double vectorialSmootherCoeff = this->smootherVectorial()->epsilon() * Order / this->mesh()->hAverage();
    vectorialSmootherParameters = "coeff (h*c/order) = " 
        + std::to_string(this->smootherVectorial()->epsilon())
        + " (" + std::to_string(this->mesh()->hAverage()) + " * " + std::to_string(vectorialSmootherCoeff) + " / " + std::to_string(Order) + ")"
        ;

    std::string restartMode = (this->doRestart())? "ON": "OFF";

    std::string exporterType = this->M_exporter->type();
    std::string hovisuMode = "OFF";
    int exporterFreq = this->M_exporter->freq();
    std::string exportedFields;
    if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::GradPhi ) )
        exportedFields = (exportedFields.empty())? "GradPhi": exportedFields+" - GradPhi";
    if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::ModGradPhi ) )
        exportedFields = (exportedFields.empty())? "ModGradPhi": exportedFields+" - ModGradPhi";
    if ( this->M_useStretchAugmented )
        exportedFields = (exportedFields.empty())? "Stretch": exportedFields+" - Stretch";

    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||---------------Info : LevelSet----------------||"
           << "\n||==============================================||"
           << "\n   Prefix          : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository()
           << "\n   Dim             : " << nDim
           << "\n   Order           : " << Order
           << "\n   Periodicity     : " << M_periodicity.isPeriodic()

           << "\n   Level Set Parameters"
           << "\n     -- thickness interface (use adaptive)  : " << this->thicknessInterface() << " (" << std::boolalpha << this->M_useAdaptiveThicknessInterface << ")"
           << "\n     -- use regular phi (phi / |grad(phi)|) : " << std::boolalpha << this->M_useRegularPhi
           << "\n     -- Heaviside/Dirac projection method   : " << hdProjectionMethod
           << "\n     -- reinit initial value                : " << std::boolalpha << this->M_reinitInitialValue
           << "\n     -- smooth gradient                     : " << std::boolalpha << this->M_doSmoothGradient
           << "\n     -- smooth curvature                    : " << std::boolalpha << this->M_doSmoothCurvature
           << "\n     -- use gradient augmented              : " << std::boolalpha << this->M_useGradientAugmented
           << "\n     -- use stretch augmented               : " << std::boolalpha << this->M_useStretchAugmented

           << "\n   Reinitialization Parameters"
           << "\n     -- reinitialization method         : " << reinitMethod;
    if( this->M_useGradientAugmented )
    *_ostr << "\n     -- reinitialize gradient augmented : " << std::boolalpha << this->M_reinitGradientAugmented;
    if( this->M_useGradientAugmented )
    *_ostr << "\n     -- reinitialize stretch augmented  : " << std::boolalpha << this->M_reinitStretchAugmented;

    *_ostr << "\n   Smoothers Parameters"
           << "\n     -- scalar smoother    : " << scalarSmootherParameters
           << "\n     -- vectorial smoother : " << vectorialSmootherParameters;

    *_ostr << "\n   Space Discretization";
    if( this->hasGeofileStr() )
    *_ostr << "\n     -- geo file name   : " << this->geofileStr();
    *_ostr << "\n     -- mesh file name  : " << this->mshfileStr()
           << "\n     -- nb elt in mesh  : " << this->mesh()->numGlobalElements()//numElements()
         //<< "\n     -- nb elt in mesh  : " << this->mesh()->numElements()
         //<< "\n     -- nb face in mesh : " << this->mesh()->numFaces()
           << "\n     -- geometry order  : " << nOrderGeo
           << "\n     -- level set order : " << Order
           << "\n     -- nb dof          : " << this->functionSpace()->nDof() << " (" << this->functionSpace()->nLocalDof() << ")"
           << "\n     -- stabilization   : " << advectionStabilization;

    *_ostr << "\n   Time Discretization"
           << "\n     -- initial time : " << this->timeStepBase()->timeInitial()
           << "\n     -- final time   : " << this->timeStepBase()->timeFinal()
           << "\n     -- time step    : " << this->timeStepBase()->timeStep()
           << "\n     -- order        : " << this->timeStepBDF()->timeOrder()
           << "\n     -- restart mode : " << restartMode
           << "\n     -- save on disk : " << std::boolalpha << this->timeStepBase()->saveInFile();
    if ( this->timeStepBase()->saveFreq() )
    *_ostr << "\n     -- freq save : " << this->timeStepBase()->saveFreq()
           << "\n     -- file format save : " << this->timeStepBase()->fileFormat();

    *_ostr << "\n   Exporter"
           << "\n     -- type            : " << exporterType
           << "\n     -- high order visu : " << hovisuMode
           << "\n     -- freq save       : " << exporterFreq
           << "\n     -- fields exported : " << exportedFields

           << "\n   Processors"
           << "\n     -- number of proc environment : " << Environment::worldComm().globalSize()
           << "\n     -- environment rank           : " << Environment::worldComm().rank()
           << "\n     -- global rank                : " << this->worldComm().globalRank()
           << "\n     -- local rank                 : " << this->worldComm().localRank()

           << "\n   Numerical Solver"
           << "\n     -- solver : " << this->solverName();

    if ( this->algebraicFactory() )
    *_ostr << this->algebraicFactory()->getInfo()->str();
    //if (enable_reinit)
    //{
        //if (reinitmethod == "hj")
        //{
            //infos << "\n      * hj maximum iteration per reinit : " << hj_max_iter
                  //<< "\n      * hj pseudo time step dtau : " << hj_dtau
                  //<< "\n      * hj stabilization : SUPG"
                  //<< "\n      * hj coeff stab : " << option( prefixvm(M_prefix,"hj-coeff-stab")).template as<double>()
                  //<< "\n      * hj tolerence on dist to dist error : "<<hj_tol;
        //}
        //else
        //{
            //infos << "\n      * fm smoothing coefficient for ILP : " << Environment::vm()[prefixvm(M_prefix,"fm-smooth-coeff")].template as< double >();
        //}
    //}
    //infos << "\n\n  Level set spaces :"
          //<< "\n     -- scalar LS space ndof : "<< this->functionSpace()->nDof()
          //<< "\n     -- vectorial LS ndof : "<< this->functionSpaceVectorial()->nDof()
          //<< "\n     -- scalar P0 space ndof : "<< this->functionSpaceMarkers()->nDof()
          //<<"\n||==============================================||\n\n";

    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n\n";

    return _ostr;
}

//----------------------------------------------------------------------------//
// Export results
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::exportResultsImpl( double time )
{
    this->log("LevelSet","exportResults", "start");
    this->timerTool("PostProcessing").start();
 
    if ( !this->M_exporter->doExport() ) return;

    this->M_exporter->step( time )->add( prefixvm(this->prefix(),"Dirac"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Dirac")),
                                   *this->dirac() );

    this->M_exporter->step( time )->add( prefixvm(this->prefix(),"Heaviside"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Heaviside")),
                                   *this->heaviside() );

    this->M_exporter->step( time )->add( prefixvm(this->prefix(),"Normal"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Normal")),
                                   *this->normal() );

    this->M_exporter->step( time )->add( prefixvm(this->prefix(),"Curvature"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Curvature")),
                                   *this->curvature() );

    if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::GradPhi ) )
    {
        this->M_exporter->step( time )->add( prefixvm(this->prefix(),"GradPhi"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"GradPhi")),
                                       *this->gradPhi() );
    }
    if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::ModGradPhi ) )
    {
        this->M_exporter->step( time )->add( prefixvm(this->prefix(),"ModGradPhi"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"ModGradPhi")),
                                       *this->modGradPhi() );
    }

    if( M_useGradientAugmented )
    {
        M_modGradPhiAdvection->exportResults( time );
    }
    if( M_useStretchAugmented )
    {
        M_stretchAdvection->exportResults( time );
        this->M_exporter->step( time )->add( prefixvm(this->prefix(),"InterfaceRectangularF"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"InterfaceRectangularF")),
                                       this->interfaceRectangularFunction() );
        this->M_exporter->step( time )->add( prefixvm(this->prefix(),"Stretch"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Stretch")),
                                       *this->stretch() );
    }

    super_type::exportResultsImpl( time );

    this->timerTool("PostProcessing").stop("exportResults");
    this->log("LevelSet","exportResults", "finish");
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::exportMeasuresImpl( double time )
{
    bool hasMeasureToExport = false;

    if( this->hasPostProcessMeasureExported( LevelSetMeasuresExported::Volume ) )
    {
        this->postProcessMeasuresIO().setMeasure( "volume", this->volume() );
        hasMeasureToExport = true;
    }
    if( this->hasPostProcessMeasureExported( LevelSetMeasuresExported::Perimeter ) )
    {
        this->postProcessMeasuresIO().setMeasure( "perimeter", this->perimeter() );
        hasMeasureToExport = true;
    }

    if( hasMeasureToExport )
    {
        this->postProcessMeasuresIO().setParameter( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
bool
LEVELSET_CLASS_TEMPLATE_TYPE::hasPostProcessMeasureExported( 
        LevelSetMeasuresExported const& measure) const
{
    return M_postProcessMeasuresExported.find(measure) != M_postProcessMeasuresExported.end();
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
bool
LEVELSET_CLASS_TEMPLATE_TYPE::hasPostProcessFieldExported( 
        LevelSetFieldsExported const& field) const
{
    return M_postProcessFieldsExported.find(field) != M_postProcessFieldsExported.end();
}

//----------------------------------------------------------------------------//
// Physical quantities
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
double
LEVELSET_CLASS_TEMPLATE_TYPE::volume() const
{
    double volume = integrate(
            _range=elements(this->mesh()),
            _expr=(1-idv(this->heaviside())) 
            ).evaluate()(0,0);
            //_expr=vf::chi( idv(this->phi())<0.0) ).evaluate()(0,0); // gives very noisy results

    return volume;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
double
LEVELSET_CLASS_TEMPLATE_TYPE::perimeter() const
{
    double perimeter = integrate(
            _range=elements(this->mesh()),
            _expr=idv(this->dirac())
            ).evaluate()(0,0);

    return perimeter;
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
        oa << BOOST_SERIALIZATION_NVP( M_vecIterSinceReinit );
    }
    this->worldComm().barrier();
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Update markers
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMarkerInterface()
{
    /* returns a marker (P0_type) on the elements crossed by the levelset
       ie :
       - if the element has not phi on all the dof of the same sign -> the element is mark as 1
       - if phi on the dof of the element are of the same sign -> mark as 0
      */

    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();
    auto phi = this->phi();

    auto it_elt = mesh->beginElementWithProcessId(mesh->worldComm().localRank());
    auto en_elt = mesh->endElementWithProcessId(mesh->worldComm().localRank());
    if (it_elt == en_elt) return;

    for (; it_elt!=en_elt; it_elt++)
    {
        int nbplus = 0;
        int nbminus = 0;

        for (int j=0; j<ndofv ; j++)
        {
            if (phi->localToGlobal(it_elt->id(), j, 0) >= 0.)
                nbplus++;
            else
                nbminus++;
        }

        //if elt crossed by interface
        if ( (nbminus != ndofv) && (nbplus!=ndofv) )
            M_markerInterface->assign(it_elt->id(), 0, 0, 1);
        else
            M_markerInterface->assign(it_elt->id(), 0, 0, 0);
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMarkerDirac()
{
    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();

    auto it_elt = mesh->beginElementWithProcessId(mesh->worldComm().localRank());
    auto en_elt = mesh->endElementWithProcessId(mesh->worldComm().localRank());
    if (it_elt == en_elt) return;

    double dirac_cut = this->dirac()->max() / 10.;

    for (; it_elt!=en_elt; it_elt++)
    {
        bool mark_elt = false;
        for (int j=0; j<ndofv; j++)
        {
            if ( std::abs( this->dirac()->localToGlobal(it_elt->id(), j, 0) ) > dirac_cut )
            {
                mark_elt = true;
                break; //don't need to do the others dof
            }
        }
        if( mark_elt )
            M_markerDirac->assign(it_elt->id(), 0, 0, 1);
        else
            M_markerDirac->assign(it_elt->id(), 0, 0, 0);
    }
}//markerDirac


LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMarkerHeaviside(bool invert, bool cut_at_half)
{
    /* returns P0 element having :
    if invert == true : 1 on elements inside Heaviside function (where H is smaller than epsilon on at least 1 dof)
    if invert == false : 1 on elements where Heaviside function greater than epsilon
    if cut_in_out == true : marker = 1 where H>0.5 and 0 where H<0.5
    */

    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();

    auto it_elt = mesh->beginElementWithProcessId(mesh->worldComm().localRank());
    auto en_elt = mesh->endElementWithProcessId(mesh->worldComm().localRank());
    if (it_elt == en_elt) return;

    double cut;
    if (cut_at_half) cut = 0.5;
    else             cut = invert ? 1e-3 : 0.999;

    for (; it_elt!=en_elt; it_elt++)
    {
        bool mark_elt = false;
        for (int j=0; j<ndofv; j++)
        {
            if ( std::abs( this->heaviside()->localToGlobal(it_elt->id(), j, 0) ) > cut )
            {
                mark_elt = true;
                break;
            }
        }
        if( mark_elt )
            M_markerHeaviside->assign(it_elt->id(), 0, 0, (invert)?0:1);
        else
            M_markerHeaviside->assign(it_elt->id(), 0, 0, (invert)?1:0);
    }
} 

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMarkerCrossedElements()
{
    // return a "marker" on the elements traversed by the interface between phio and phi
    // ie : mark as 1 is one of the dof of the element has  (phi * phio < 0)
    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();

    auto it_elt = mesh->beginElementWithProcessId(mesh->worldComm().localRank());
    auto en_elt = mesh->endElementWithProcessId(mesh->worldComm().localRank());
    if (it_elt == en_elt) return;

    auto phi = this->phi();
    auto phio = this->phio();

    auto prod = vf::project(this->functionSpace(), elements(mesh),
                            idv(phio) * idv(phi) );


    for (; it_elt!=en_elt; it_elt++)
    {
        bool mark_elt = false;
        for (int j=0; j<ndofv ; j++)
        {
            if (prod.localToGlobal(it_elt->id(), j, 0) <= 0.)
            {
                mark_elt = true;
                break;
            }
        }
        if( mark_elt )
            M_markerCrossedElements->assign(it_elt->id(), 0, 0, 1);
        else
            M_markerCrossedElements->assign(it_elt->id(), 0, 0, 0);
    }
}

} // FeelModels
} // Feel
