#include <feel/feelmodels/levelset/levelset.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>
#include <feel/feelmodels/levelset/reinitializer_hj.hpp>

#include <feel/feelmodels/levelset/cauchygreeninvariantsexpr.hpp>
#include <feel/feelmodels/levelset/levelsetdeltaexpr.hpp>

#include <boost/assign/list_of.hpp>
#include <boost/assign/list_inserter.hpp>

namespace Feel {
namespace FeelModels {

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
std::map<std::string, typename LEVELSET_CLASS_TEMPLATE_TYPE::ShapeType>
LEVELSET_CLASS_TEMPLATE_TYPE::ShapeTypeMap = {
    {"sphere", ShapeType::SPHERE},
    {"ellipse", ShapeType::ELLIPSE}
};

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
const typename LEVELSET_CLASS_TEMPLATE_TYPE::fastmarchinginitializationmethodidmap_type
LEVELSET_CLASS_TEMPLATE_TYPE::FastMarchingInitializationMethodIdMap = boost::assign::list_of< typename LEVELSET_CLASS_TEMPLATE_TYPE::fastmarchinginitializationmethodidmap_type::relation >
    ( "none", FastMarchingInitializationMethod::NONE )
    ( "ilp", FastMarchingInitializationMethod::ILP )
    ( "smoothed_ilp", FastMarchingInitializationMethod::SMOOTHED_ILP )
    ( "hj", FastMarchingInitializationMethod::HJ_EQ )
    ( "il-hj", FastMarchingInitializationMethod::IL_HJ_EQ )
;

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
    M_doUpdateInterfaceElements(true),
    M_doUpdateSmootherInterface(true),
    M_doUpdateSmootherInterfaceVectorial(true),
    M_doUpdateNormal(true),
    M_doUpdateCurvature(true),
    M_doUpdateGradPhi(true),
    M_doUpdateModGradPhi(true),
    M_doUpdateSubmeshDirac(true),
    M_doUpdateSubmeshOuter(true),
    M_doUpdateSubmeshInner(true),
    M_advectionToolbox( new advection_toolbox_type( prefix, worldComm, subPrefix, rootRepository ) ),
    M_doUpdateMarkers(true),
    M_doUpdateCauchyGreenTensor(true),
    M_doUpdateCauchyGreenInvariant1(true),
    M_doUpdateCauchyGreenInvariant2(true),
    //M_periodicity(periodicityLS),
    M_reinitializerIsUpdatedForUse(false),
    M_hasReinitialized(false),
    M_hasReinitializedSmooth(false),
    M_iterSinceReinit(0)
{
    this->setFilenameSaveInfo( prefixvm(this->prefix(),"Levelset.info") );
    //-----------------------------------------------------------------------------//
    // Set advection model
    M_advectionToolbox->setModelName( "Advection" );
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

//----------------------------------------------------------------------------//
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::init()
{
    this->log("LevelSet", "init", "start");

    // Init levelset advection
    if( !this->doRestart() )
    {
        // Set levelset initial value
        this->initLevelsetValue();
    }
    M_advectionToolbox->init();
    M_timeOrder = this->timeStepBDF()->timeOrder();
    this->updateTime( M_advectionToolbox->currentTime() );
    if (this->doRestart())
        this->setTimeInitial( M_advectionToolbox->timeInitial() );

    M_initialVolume = this->volume();

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
    // Init backward characteristics advection
    if( M_useCauchyAugmented )
    {
        M_backwardCharacteristicsAdvection->init();
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

    if ( M_doSmoothGradient )
        auto smootherbuild = this->smootherVectorial();
    if ( M_doSmoothCurvature )
        auto smootherbuild = this->smoother();

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
            //auto modGradPhiInit = this->smoother()->project( sqrt( trans(idv(gradPhiInit))*idv(gradPhiInit) ) );
                
            phi_init = vf::project( 
                _space=this->functionSpace(),
                _range=elements(this->mesh()),
                //_expr=idv(phi_init) / idv(modGradPhiInit)
                _expr=idv(phi_init) / sqrt( trans(idv(gradPhiInit))*idv(gradPhiInit) )
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
    if( M_useCauchyAugmented )
    {
        // Initialize backward characteristics
        if( M_hasInitialBackwardCharacteristics )
        {
            *(M_backwardCharacteristicsAdvection->fieldSolutionPtr()) = vf::project(
                    _space=M_backwardCharacteristicsAdvection->functionSpace(),
                    _range=elements(M_backwardCharacteristicsAdvection->mesh()),
                    _expr=M_initialBackwardCharacteristics
                    );
        }
        else
        {
            *(M_backwardCharacteristicsAdvection->fieldSolutionPtr()) = vf::project(
                    _space=M_backwardCharacteristicsAdvection->functionSpace(),
                    _range=elements(M_backwardCharacteristicsAdvection->mesh()),
                    _expr=vf::P()
                    );
        }
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
        0.5*(1 + (phi+epsilon_zero)/epsilon_delta + 1/M_PI*vf::sin( M_PI*(phi+epsilon_zero)/epsilon_delta ) )
        +
        vf::chi( phi>=-epsilon )*vf::chi( phi<=epsilon )*vf::constant(1.0)
        +
        vf::chi( phi>=epsilon )*vf::chi( phi<=epsilon_rect )*
        0.5*(1 - (phi-epsilon_zero)/epsilon_delta - 1/M_PI*vf::sin( M_PI*(phi-epsilon_zero)/epsilon_delta ) )
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
    if (this->doRestart() && this->restartPath().empty() )
    {
        if ( M_exporter->doExport() ) M_exporter->restart(this->timeInitial());
    }

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
    if( this->hasPostProcessMeasureExported( LevelSetMeasuresExported::Position_COM ) )
    {
        auto com = this->positionCOM();
        std::vector<double> vecCOM = { com(0,0) };
        if( nDim > 1 ) vecCOM.push_back( com(1,0) );
        if( nDim > 2 ) vecCOM.push_back( com(2,0) );
        this->postProcessMeasuresIO().setMeasureComp( "position_com", vecCOM );
        hasMeasureToExport = true;
    }
    if( this->hasPostProcessMeasureExported( LevelSetMeasuresExported::Velocity_COM ) )
    {
        auto ucom = this->velocityCOM();
        std::vector<double> vecUCOM = { ucom(0,0) };
        if( nDim > 1 ) vecUCOM.push_back( ucom(1,0) );
        if( nDim > 2 ) vecUCOM.push_back( ucom(2,0) );
        this->postProcessMeasuresIO().setMeasureComp( "velocity_com", vecUCOM );
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
LEVELSET_CLASS_TEMPLATE_TYPE::createFunctionSpaces( bool buildSpaceMarkersExtendedDofTable )
{
    M_spaceLevelSetVec = space_levelset_vectorial_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
    M_spaceMarkers = space_markers_type::New( 
            _mesh=this->mesh(), _worldscomm=this->worldsComm(),
            _extended_doftable=std::vector<bool>(1, buildSpaceMarkersExtendedDofTable)
            );
    if( M_useCauchyAugmented )
        M_spaceTensor2Symm = space_tensor2symm_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
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
    if( M_useCauchyAugmented )
    {
        M_backwardCharacteristicsAdvection = backwardcharacteristics_advection_type::New(
                prefixvm(this->prefix(), "backward-characteristics-advection"),
                this->worldComm()
                );
        M_backwardCharacteristicsAdvection->setModelName( "Advection" );
        M_backwardCharacteristicsAdvection->build( this->functionSpaceVectorial() );
        M_backwardCharacteristicsAdvection->timeStepBDF()->setOrder( this->timeStepBDF()->bdfOrder() );

        M_backwardCharacteristicsAdvection->getExporter()->setDoExport( boption( _name="do_export_backward-characteristics-advection", _prefix=this->prefix() ) );

        M_leftCauchyGreenTensor_K.reset( new element_tensor2symm_type(this->functionSpaceTensor2Symm(), "LeftCauchyGreenTensor_K") );
        M_leftCauchyGreenTensor_KN.reset( new element_levelset_vectorial_type(this->functionSpaceVectorial(), "LeftCauchyGreenTensor_KN") );
        M_leftCauchyGreenTensor.reset( new element_tensor2symm_type(this->functionSpaceTensor2Symm(), "LeftCauchyGreenTensor") );
        M_cauchyGreenInvariant1.reset( new element_cauchygreen_invariant_type(this->functionSpace(), "CauchyGreenI1(TrC)") );
        M_cauchyGreenInvariant2.reset( new element_cauchygreen_invariant_type(this->functionSpace(), "CauchyGreenI2(TrCofC)") );
    }
    for( std::string const& bcMarker: this->M_bcMarkersInflow )
    {
        if( M_useGradientAugmented )
            M_modGradPhiAdvection->addMarkerInflowBC( bcMarker );
        if( M_useStretchAugmented )
            M_stretchAdvection->addMarkerInflowBC( bcMarker );
        if( M_useCauchyAugmented )
            M_backwardCharacteristicsAdvection->addMarkerInflowBC( bcMarker );
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
            boost::dynamic_pointer_cast<reinitializerHJ_type>(M_reinitializer)->setThicknessHeaviside( thickness_heaviside );
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
    if( M_useCauchyAugmented )
    {
        M_projectorL2Tensor2Symm = projector( 
                this->functionSpaceTensor2Symm() , this->functionSpaceTensor2Symm(), 
                backend(_name=prefixvm(this->prefix(),"projector-l2-tensor2symm"), _worldcomm=this->worldComm())
                );
    }

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
void
LEVELSET_CLASS_TEMPLATE_TYPE::createExporters()
{
    std::string geoExportType = "static";//this->geoExportType();//change_coords_only, change, static
    M_exporter = exporter( 
            _mesh=this->mesh(),
            _name="ExportLS",
            _geo=geoExportType,
            _path=this->exporterPath() 
            );
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
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_backwardcharacteristics_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::backwardCharacteristics() const
{
    if( M_useCauchyAugmented )
    {
        return M_backwardCharacteristicsAdvection->fieldSolutionPtr();
    }
    else
    {
        throw std::logic_error( this->prefix()+".use-cauchy-augmented option must be true to use backward characteristics" );
    }
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
LEVELSET_CLASS_TEMPLATE_TYPE::setFastMarchingInitializationMethod( FastMarchingInitializationMethod m )
{
    if (M_reinitializerIsUpdatedForUse)
        LOG(INFO)<<" !!!  WARNING !!! : fastMarchingInitializationMethod set after the fast marching has been actually initialized ! \n";
    M_fastMarchingInitializationMethod = m;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
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

    const std::string fm_init_method = soption( _name="fm-initialization-method", _prefix=this->prefix() );
    CHECK(FastMarchingInitializationMethodIdMap.left.count(fm_init_method)) << fm_init_method <<" is not in the list of possible fast-marching initialization methods\n";
    M_fastMarchingInitializationMethod = FastMarchingInitializationMethodIdMap.left.at(fm_init_method);

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
            if( o == "position_com" || o == "all" )
                this->M_postProcessMeasuresExported.insert( LevelSetMeasuresExported::Position_COM );
            if( o == "velocity_com" || o == "all" )
                this->M_postProcessMeasuresExported.insert( LevelSetMeasuresExported::Velocity_COM );
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
            if( o == "backwardcharacteristics" )
                this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::BackwardCharacteristics );
            if( o == "cauchygreeninvariant1" )
                this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::CauchyGreenInvariant1 );
            if( o == "cauchygreeninvariant2" )
                this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::CauchyGreenInvariant2 );
        }
    }
    // Overwrite with options from CFG
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_gradphi").c_str()) )
        if ( boption(_name="do_export_gradphi",_prefix=this->prefix()) )
            this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::GradPhi );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_modgradphi").c_str()) )
        if ( boption(_name="do_export_modgradphi",_prefix=this->prefix()) )
            this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::ModGradPhi );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_backwardcharacteristics").c_str()) )
        if ( boption(_name="do_export_backwardcharacteristics",_prefix=this->prefix()) )
            this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::BackwardCharacteristics );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_cauchygreeninvariant1").c_str()) )
        if ( boption(_name="do_export_cauchygreeninvariant1",_prefix=this->prefix()) )
            this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::CauchyGreenInvariant1 );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_cauchygreeninvariant2").c_str()) )
        if ( boption(_name="do_export_cauchygreeninvariant2",_prefix=this->prefix()) )
            this->M_postProcessFieldsExported.insert( LevelSetFieldsExported::CauchyGreenInvariant2 );
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
    this->log("LevelSet", "updateGradPhi", "start");
    this->timerTool("UpdateInterfaceData").start();

    auto phi = this->phi();
    if( M_doSmoothGradient )
    {
        this->log("LevelSet", "updateGradPhi", "perform smooth projection");
        //*M_levelsetGradPhi = this->smootherVectorial()->derivate( idv(phi) );
        *M_levelsetGradPhi = this->smootherVectorial()->project( trans(gradv(phi)) );
    }
    else
    {
        //*M_levelsetGradPhi = this->projectorL2Vectorial()->derivate( idv(phi) );
        this->log("LevelSet", "updateGradPhi", "perform L2 projection");
        *M_levelsetGradPhi = this->projectorL2Vectorial()->project( _expr=trans(gradv(phi)) );
        //*M_levelsetGradPhi = vf::project( 
                //_space=this->functionSpaceVectorial(),
                //_range=elements(this->mesh()),
                //_expr=trans(gradv(phi)) 
                //);
    }

    M_doUpdateGradPhi = false;
    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSet", "updateGradPhi", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateModGradPhi()
{
    this->log("LevelSet", "updateModGradPhi", "start");
    this->timerTool("UpdateInterfaceData").start();

    auto gradPhi = this->gradPhi();
    *M_levelsetModGradPhi = vf::project( 
            _space=this->functionSpace(),
            _range=elements(this->mesh()),
            _expr=sqrt( trans(idv(gradPhi))*idv(gradPhi) ) 
            );
    //*M_levelsetModGradPhi = this->projectorL2()->project( _expr=sqrt( trans(idv(gradPhi))*idv(gradPhi) ) );

    M_doUpdateModGradPhi = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSet", "updateModGradPhi", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateDirac()
{
    this->log("LevelSet", "updateDirac", "start");
    this->timerTool("UpdateInterfaceData").start();

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
        //auto D_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            //+
            //vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            //1/(2*eps) *( 1 + cos(M_PI*psi/eps) )
            //+
            //vf::chi(psi>eps)*vf::constant(0.0);

        if ( M_useHeavisideDiracNodalProj )
            *M_dirac = vf::project( this->functionSpace(), elements(this->mesh()),
                   Feel::FeelModels::levelsetDelta(psi, eps0) );
        else
            *M_dirac = M_projectorL2->project( Feel::FeelModels::levelsetDelta(psi, eps0) );
    }
    else
    {
        auto psi = idv(this->phi()) ;
        //auto D_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            //+
            //vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            //1/(2*eps) *( 1 + cos(M_PI*psi/eps) )
            //+
            //vf::chi(psi>eps)*vf::constant(0.0);

        if ( M_useHeavisideDiracNodalProj )
            *M_dirac = vf::project( this->functionSpace(), elements(this->mesh()),
                   Feel::FeelModels::levelsetDelta(psi, eps0) );
        else
            *M_dirac = M_projectorL2->project( Feel::FeelModels::levelsetDelta(psi, eps0) );
    }

    M_doUpdateDirac = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSet", "updateDirac", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateHeaviside()
{ 
    this->log("LevelSet", "updateHeaviside", "start");
    this->timerTool("UpdateInterfaceData").start();

    auto eps = this->thicknessInterface();

    if (M_useRegularPhi)
    {
        auto psi = idv(this->phi()) / idv(this->modGradPhi());
        auto H_expr = vf::chi( psi<-eps )*vf::constant(0.0)
            +
            vf::chi( psi>=-eps )*vf::chi( psi<=eps )*
            0.5*(1 + psi/eps + 1/M_PI*vf::sin( M_PI*psi/eps ) )
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
            0.5*(1 + psi/eps + 1/M_PI*vf::sin( M_PI*psi/eps ) )
            +
            vf::chi(psi>eps)*vf::constant(1.0);

        if (M_useHeavisideDiracNodalProj)
            *M_heaviside = vf::project(this->functionSpace(), elements(this->mesh()), H_expr);
        else
            *M_heaviside = M_projectorL2->project(H_expr);
    }

    M_doUpdateHeaviside = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSet", "updateHeaviside", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateNormal()
{
    this->log("LevelSet", "updateNormal", "start");
    this->timerTool("UpdateInterfaceData").start();

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

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSet", "updateNormal", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateCurvature()
{
    this->log("LevelSet", "updateCurvature", "start");
    this->timerTool("UpdateInterfaceData").start();

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

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSet", "updateCurvature", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
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

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::projector_levelset_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::smootherInterface() const
{
    if( !M_smootherInterface || M_doUpdateSmootherInterface )
    {
        auto const spaceInterface = self_type::space_levelset_type::New( 
                _mesh=this->mesh(),
                _range=this->interfaceElements(),
                _worldscomm=this->worldsComm()
                );
        M_smootherInterface = Feel::projector( 
                spaceInterface, spaceInterface,
                backend(_name=prefixvm(this->prefix(),"smoother"), _worldcomm=this->worldComm(), _rebuild=true), 
                DIFF,
                this->mesh()->hAverage()*doption(_name="smooth-coeff", _prefix=prefixvm(this->prefix(),"smoother"))/Order,
                30);
        M_doUpdateSmootherInterface = false;
    }
    return M_smootherInterface;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::projector_levelset_vectorial_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::smootherInterfaceVectorial() const
{
    if( !M_smootherInterfaceVectorial || M_doUpdateSmootherInterfaceVectorial )
    {
        auto const spaceInterfaceVectorial = self_type::space_levelset_vectorial_type::New( 
                _mesh=this->mesh(),
                _range=this->interfaceElements(),
                _worldscomm=this->worldsComm()
                );
        M_smootherInterfaceVectorial = Feel::projector(
                spaceInterfaceVectorial, spaceInterfaceVectorial,
                backend(_name=prefixvm(this->prefix(),"smoother-vec"), _worldcomm=this->worldComm(), _rebuild=true),
                DIFF, 
                this->mesh()->hAverage()*doption(_name="smooth-coeff", _prefix=prefixvm(this->prefix(),"smoother-vec"))/Order,
                30);
        M_doUpdateSmootherInterfaceVectorial = false;
    }
    return M_smootherInterfaceVectorial;
}
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Markers accessors
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerInterface() const
{
    if( !M_markerInterface )
        M_markerInterface.reset( new element_markers_type(M_spaceMarkers, "MarkerInterface") );

    if( M_doUpdateMarkers )
       const_cast<self_type*>(this)->updateMarkerInterface(); 

    return M_markerInterface;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerDirac() const
{
    if( !M_markerDirac )
        M_markerDirac.reset( new element_markers_type(M_spaceMarkers, "MarkerDirac") );

    if( M_doUpdateMarkers )
       const_cast<self_type*>(this)->updateMarkerDirac(); 

    return M_markerDirac;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerOuter( double cut ) const
{
    if( !M_markerOuter )
        M_markerOuter.reset( new element_markers_type(M_spaceMarkers, "MarkerOuter") );

    if( M_doUpdateMarkers || cut != M_markerOuterCut )
    {
       const_cast<self_type*>(this)->markerHeavisideImpl( M_markerOuter, false, cut );
       M_markerOuterCut = cut;
    }

    return M_markerOuter;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerInner( double cut ) const
{
    if( !M_markerInner )
        M_markerInner.reset( new element_markers_type(M_spaceMarkers, "MarkerInner") );

    if( M_doUpdateMarkers || cut != M_markerInnerCut )
    {
       const_cast<self_type*>(this)->markerHeavisideImpl( M_markerInner, true, cut );
       M_markerInnerCut = cut;
    }

    return M_markerInner;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::markerCrossedElements() const
{
    if( !M_markerCrossedElements )
        M_markerCrossedElements.reset( new element_markers_type(M_spaceMarkers, "MarkerCrossedElements") );

    if( M_doUpdateMarkers )
       const_cast<self_type*>(this)->updateMarkerCrossedElements(); 

    return M_markerCrossedElements;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::range_elements_type
LEVELSET_CLASS_TEMPLATE_TYPE::interfaceElements() const
{
    if( this->M_doUpdateInterfaceElements )
    {
        mesh_ptrtype const& mesh = this->mesh();

        auto it_elt = mesh->beginOrderedElement();
        auto en_elt = mesh->endOrderedElement();

        const rank_type pid = mesh->worldCommElements().localRank();
        const int ndofv = space_levelset_type::fe_type::nDof;

        double thickness = 2*this->thicknessInterface();
        elements_reference_wrapper_ptrtype interfaceElts( new elements_reference_wrapper_type );

        for (; it_elt!=en_elt; it_elt++)
        {
            auto const& elt = boost::unwrap_ref( *it_elt );
            if ( elt.processId() != pid )
                continue;
            bool mark_elt = false;
            for (int j=0; j<ndofv; j++)
            {
                if ( std::abs( this->phi()->localToGlobal(elt.id(), j, 0) ) <= thickness )
                {
                    mark_elt = true;
                    break; //don't need to do the others dof
                }
            }
            if( mark_elt )
                interfaceElts->push_back( boost::cref(elt) );
        }

        M_interfaceElements = boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                interfaceElts->begin(),
                interfaceElts->end(),
                interfaceElts
                );

        M_doUpdateInterfaceElements = false;
    }

    return M_interfaceElements;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::range_elements_type
LEVELSET_CLASS_TEMPLATE_TYPE::outerElementsRange( double cut )
{
    mesh_ptrtype const& mesh = this->mesh();

    element_levelset_type phi = this->functionSpace()->element();
    if( this->M_useRegularPhi )
        phi.on( _range=elements(this->mesh()), _expr=idv(this->phi()) / idv(this->modGradPhi()) );
    else
        phi = *(this->phi());

    auto it_elt = mesh->beginOrderedElement();
    auto en_elt = mesh->endOrderedElement();

    const rank_type pid = mesh->worldCommElements().localRank();
    const int ndofv = space_levelset_type::fe_type::nDof;

    elements_reference_wrapper_ptrtype outerElts( new elements_reference_wrapper_type );

    for (; it_elt!=en_elt; it_elt++)
    {
        auto const& elt = boost::unwrap_ref( *it_elt );
        if ( elt.processId() != pid )
            continue;
        bool mark_elt = true;
        for (int j=0; j<ndofv; j++)
        {
            if ( phi.localToGlobal(elt.id(), j, 0) < cut )
            {
                mark_elt = false;
                break; //don't need to do the others dof
            }
        }
        if( mark_elt )
            outerElts->push_back( boost::cref(elt) );
    }

    range_elements_type outerElements = boost::make_tuple( 
            mpl::size_t<MESH_ELEMENTS>(),
            outerElts->begin(),
            outerElts->end(),
            outerElts
            );

    return outerElements;
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

    typedef boost::reference_wrapper<typename MeshTraits<mymesh_type>::element_type const> element_ref_type;
    typedef std::vector<element_ref_type> cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    // Retrieve the elements touching the marked faces
    auto mfaces = markedfaces( this->mesh(), marker );
    for( auto const& faceWrap: mfaces )
    {
        auto const& face = unwrap_ref( faceWrap );
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

    typedef boost::reference_wrapper<typename MeshTraits<mymesh_type>::element_type const> element_ref_type;
    typedef std::vector<element_ref_type> cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    // Retrieve the elements touching the marked faces
    //auto mfaces_list = markedfaces( this->mesh(), marker );
    auto mfaces = markedfaces( this->mesh(), marker );
    for( auto const& faceWrap: mfaces )
    {
        auto const& face = unwrap_ref( faceWrap );
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

}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("LevelSet", "solve", "start");
    this->timerTool("Solve").start();

    // Solve phi
    M_advectionToolbox->solve();

    if( M_useGradientAugmented )
    {
        // Solve modGradPhi
        //auto modGradPhi = M_modGradPhiAdvection->fieldSolutionPtr();
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

    // Update interface-related quantities
    this->updateInterfaceQuantities();

    // Correct volume if requested
    if( this->M_fixVolume )
    {
        auto const& phi = this->phi();
        double lambda = ( this->volume() - this->M_initialVolume ) / this->perimeter();
        phi->add( lambda );
        // Request update interface-related quantities again since phi has changed
        // Note that updateInterfaceQuantities has lazy evaluation
        this->updateInterfaceQuantities();
    }

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

    M_advectionToolbox->updateTimeStep();
    if( M_useGradientAugmented )
        M_modGradPhiAdvection->updateTimeStep();
    if( M_useStretchAugmented )
        M_stretchAdvection->updateTimeStep();
    if( M_useCauchyAugmented )
        M_backwardCharacteristicsAdvection->updateTimeStep();

    this->updateTime( M_advectionToolbox->currentTime() );

    if( M_iterSinceReinit < M_timeOrder )
    {
        this->timeStepBDF()->setTimeOrder( M_iterSinceReinit + 1 );
        if( M_useGradientAugmented )
            M_modGradPhiAdvection->timeStepBDF()->setTimeOrder( M_iterSinceReinit + 1 );
        if( M_useStretchAugmented )
            M_stretchAdvection->timeStepBDF()->setTimeOrder( M_iterSinceReinit + 1 );
        if( M_useCauchyAugmented )
            M_backwardCharacteristicsAdvection->timeStepBDF()->setTimeOrder( M_iterSinceReinit + 1 );
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
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateLeftCauchyGreenTensor()
{
    DCHECK( this->M_useCauchyAugmented ) << this->prefix()+".use-cauchy-augmented option must be true to use Cauchy-Green tensor";

    this->log("LevelSet", "updateLeftCauchyGreenTensor", "start");
    this->timerTool("UpdateInterfaceData").start();

#if 0
    // Create interface projector L2
    auto const interfaceElts = this->interfaceElements();
    auto const spaceTensor2SymmInterface = self_type::space_tensor2symm_type::New( 
            _mesh=this->mesh(),
            _range=interfaceElts,
            _worldscomm=this->worldsComm()
            );
    auto const projectorL2Tensor2SymmInterface = Feel::projector(
            spaceTensor2SymmInterface, spaceTensor2SymmInterface,
            backend(_name=prefixvm(this->prefix(),"projector-l2-tensor2symm"), _worldcomm=this->worldComm(), _rebuild=true )
            );

    auto Y = M_backwardCharacteristicsAdvection->fieldSolutionPtr();
    auto gradY = projectorL2Tensor2SymmInterface->project(
            _expr=gradv(Y)
            );
    auto invGradY = this->functionSpaceTensor2Symm()->element();
    invGradY.on(
            _range=interfaceElts,
            _expr=inv(idv(gradY))
            );
    // K = (gradY)^-1 (gradY)^-T
    auto const& N = this->N();
    M_leftCauchyGreenTensor_K->zero();
    M_leftCauchyGreenTensor_K->on(
            _range=interfaceElts,
            _expr=idv(invGradY)*trans(idv(invGradY))
            );
    auto const& K = *M_leftCauchyGreenTensor_K;
    M_leftCauchyGreenTensor_KN->zero();
    M_leftCauchyGreenTensor_KN->on(
            _range=interfaceElts,
            _expr=idv(K)*idv(this->N())
            );
    auto const& KN = *M_leftCauchyGreenTensor_KN;
    M_leftCauchyGreenTensor->zero();
    M_leftCauchyGreenTensor->on(
            _range=interfaceElts,
            _expr=idv(K) - idv(KN)*trans(idv(KN))/(trans(idv(N))*idv(KN))
            );
#else
    auto Y = M_backwardCharacteristicsAdvection->fieldSolutionPtr();
    auto gradY = this->projectorL2Tensor2Symm()->project(
            _expr=gradv(Y)
            );
    auto invGradY = vf::project(
            _space=this->functionSpaceTensor2Symm(),
            _range=elements(this->mesh()),
            _expr=inv(idv(gradY))
            );
    // K = (gradY)^-1 (gradY)^-T
    auto const& N = this->N();
    M_leftCauchyGreenTensor_K->on(
            _range=elements(this->mesh()),
            _expr=idv(invGradY)*trans(idv(invGradY))
            );
    auto const& K = *M_leftCauchyGreenTensor_K;
    M_leftCauchyGreenTensor_KN->on(
            _range=elements(this->mesh()),
            _expr=idv(K)*idv(this->N())
            );
    auto const& KN = *M_leftCauchyGreenTensor_KN;
    M_leftCauchyGreenTensor->on(
            _range=elements(this->mesh()),
            _expr=idv(K) - idv(KN)*trans(idv(KN))/(trans(idv(N))*idv(KN))
            );
#endif

    M_doUpdateCauchyGreenTensor = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSet", "updateLeftCauchyGreenTensor", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
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
                _range=this->interfaceElements(),
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
                _range=this->interfaceElements(),
                _expr=Feel::vf::FeelModels::cauchyGreenInvariant1Expr( A )
                );
#endif
        M_doUpdateCauchyGreenInvariant1 = false;

        double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
        this->log("LevelSet", "cauchyGreenInvariant1", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
    }

    return M_cauchyGreenInvariant1;
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
                _range=this->interfaceElements(),
                _expr=det(idv(K))/(trans(idv(N))*idv(KN))
                );
//#elif 1 // New implementation TrA / (2 sqrt(cofA))
#elif 1 // Old implementation TrA/2
        auto A = idv(this->leftCauchyGreenTensor());
        auto trA = trace(A);
        M_cauchyGreenInvariant2->zero();
        M_cauchyGreenInvariant2->on(
                _range=this->interfaceElements(),
                _expr=0.5 * trA
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
    M_doUpdateDirac = true;
    M_doUpdateHeaviside = true;
    M_doUpdateInterfaceElements = true;
    M_doUpdateSmootherInterface = true;
    M_doUpdateSmootherInterfaceVectorial = true;
    M_doUpdateNormal = true;
    M_doUpdateCurvature = true;
    M_doUpdateMarkers = true;
    M_doUpdateGradPhi = true;
    M_doUpdateModGradPhi = true;
    M_doUpdateSubmeshDirac = true;
    M_doUpdateSubmeshOuter = true;
    M_doUpdateSubmeshInner = true;
    M_doUpdateCauchyGreenTensor = true;
    M_doUpdateCauchyGreenInvariant1 = true;
    M_doUpdateCauchyGreenInvariant2 = true;
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

        switch (M_fastMarchingInitializationMethod)
        {
            case FastMarchingInitializationMethod::ILP :
            {
                auto const gradPhi = idv(this->gradPhi());
                
                *phiReinit = vf::project(
                        this->functionSpace(), 
                        elements(this->mesh()), 
                        idv(phi)/ sqrt( trans(gradPhi) * gradPhi )
                        );
            }
            break;

            case FastMarchingInitializationMethod::SMOOTHED_ILP :
            {
                // save the smoothed gradient magnitude of phi
                //auto modgradphi = M_smootherFM->project( vf::min(vf::max(vf::sqrt(inner(gradv(phi), gradv(phi))), 0.92), 2.) );
                //auto gradPhi = idv(this->gradPhi());
                auto gradPhi = trans(gradv(phi));
                //auto modgradphi = M_smootherFM->project( sqrt( trans(gradPhi)*gradPhi ) );
                auto modgradphi = this->smoother()->project( sqrt( trans(gradPhi)*gradPhi ) );
                
                *phiReinit = vf::project(
                        this->functionSpace(), 
                        elements(this->mesh()), 
                        idv(phi)/idv(modgradphi) 
                        );
            }
            break;

            case FastMarchingInitializationMethod::HJ_EQ :
            {
                CHECK(false) << "TODO\n";
                //*phi = *explicitHJ(max_iter, dtau, tol);
            }
            break;
            case FastMarchingInitializationMethod::IL_HJ_EQ :
            {
                CHECK(false) << "TODO\n";
                //*phi = *explicitHJ(max_iter, dtau, tol);
            }
            break;
            case FastMarchingInitializationMethod::NONE :
            {
                *phiReinit = *phi;
            }
            break;
        } // switch M_fastMarchingInitializationMethod

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
    std::string advectionStabilization = soption( _name="stabilization.method", _prefix=this->prefix() );

    std::string hdProjectionMethod = (this->M_useHeavisideDiracNodalProj)? "nodal": "L2";

    std::string reinitMethod;
    std::string reinitmethod = soption( _name="reinit-method", _prefix=this->prefix() );
    if( reinitmethod == "fm" )
    {
        reinitMethod = "Fast-Marching";
        std::string fmInitMethod = FastMarchingInitializationMethodIdMap.right.at( this->M_fastMarchingInitializationMethod );
        reinitMethod += " (" + fmInitMethod + ")";
    }
    else if( reinitmethod == "hj" )
        reinitMethod = "Hamilton-Jacobi";

    std::string scalarSmootherParameters;
    if ( M_smoother )
    {
        double scalarSmootherCoeff = this->smoother()->epsilon() * Order / this->mesh()->hAverage();
        scalarSmootherParameters = "coeff (h*c/order) = " 
            + std::to_string(this->smoother()->epsilon())
            + " (" + std::to_string(this->mesh()->hAverage()) + " * " + std::to_string(scalarSmootherCoeff) + " / " + std::to_string(Order) + ")"
            ;
    }
    std::string vectorialSmootherParameters;

    if ( M_smootherVectorial )
    {
        double vectorialSmootherCoeff = this->smootherVectorial()->epsilon() * Order / this->mesh()->hAverage();
        vectorialSmootherParameters = "coeff (h*c/order) = " 
            + std::to_string(this->smootherVectorial()->epsilon())
            + " (" + std::to_string(this->mesh()->hAverage()) + " * " + std::to_string(vectorialSmootherCoeff) + " / " + std::to_string(Order) + ")"
            ;
    }

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

    if( M_smoother || M_smootherVectorial )
    *_ostr << "\n   Smoothers Parameters";
    if( M_smoother )
    *_ostr << "\n     -- scalar smoother    : " << scalarSmootherParameters;
    if( M_smootherVectorial )
    *_ostr << "\n     -- vectorial smoother : " << vectorialSmootherParameters;

    *_ostr << "\n   Space Discretization";
    if( this->hasGeofileStr() )
    *_ostr << "\n     -- geo file name   : " << this->geofileStr();
    *_ostr << "\n     -- mesh file name  : " << this->mshfileStr()
           << "\n     -- nb elt in mesh  : " << this->mesh()->numGlobalElements()//numElements()
         //<< "\n     -- nb elt in mesh  : " << this->mesh()->numElements()
         //<< "\n     -- nb face in mesh : " << this->mesh()->numFaces()
           << "\n     -- hMin            : " << this->mesh()->hMin()
           << "\n     -- hMax            : " << this->mesh()->hMax()
           << "\n     -- hAverage        : " << this->mesh()->hAverage()
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
           << "\n     -- solver : " << M_advectionToolbox->solverName();

#if 0
    if ( this->algebraicFactory() )
    *_ostr << this->algebraicFactory()->getInfo()->str();
#endif
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
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::mesh_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::submeshDirac() const
{
    if( M_doUpdateSubmeshDirac )
    {
        this->mesh()->updateMarker2( *this->markerDirac() );
        M_submeshDirac = createSubmesh( this->mesh(), marked2elements( this->mesh(), 1 ) );
        M_doUpdateSubmeshDirac = false;
    }
    return M_submeshDirac;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::mesh_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::submeshOuter( double cut ) const
{
    if( M_doUpdateSubmeshOuter || cut != M_markerOuterCut )
    {
        this->mesh()->updateMarker2( *this->markerOuter( cut ) );
        M_submeshOuter = createSubmesh( this->mesh(), marked2elements( this->mesh(), 1 ) );
        M_doUpdateSubmeshOuter = false;
    }
    return M_submeshOuter;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSET_CLASS_TEMPLATE_TYPE::mesh_ptrtype const&
LEVELSET_CLASS_TEMPLATE_TYPE::submeshInner( double cut ) const
{
    if( M_doUpdateSubmeshInner || cut != M_markerInnerCut )
    {
        this->mesh()->updateMarker2( *this->markerInner( cut ) );
        M_submeshInner = createSubmesh( this->mesh(), marked2elements( this->mesh(), 1 ) );
        M_doUpdateSubmeshInner = false;
    }
    return M_submeshInner;
}

//----------------------------------------------------------------------------//
// Export results
LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->exportResultsImpl( time );
    this->exportMeasuresImpl( time );
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::exportResultsImpl( double time )
{
    this->log("LevelSet","exportResults", "start");
    this->timerTool("PostProcessing").start();

    if( this->M_doExportAdvection )
        this->M_advectionToolbox->exportResults( time );

    if ( !this->M_exporter->doExport() ) return;

    this->M_exporter->step( time )->add( prefixvm(this->prefix(),"Phi"),
                                         prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Phi")),
                                         *this->phi() );
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
    if( M_useCauchyAugmented )
    {
        M_backwardCharacteristicsAdvection->exportResults( time );
    }
    if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::BackwardCharacteristics ) )
    {
        this->M_exporter->step( time )->add( prefixvm(this->prefix(),"BackwardCharacteristics"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"BackwardCharacteristics")),
                                       M_backwardCharacteristicsAdvection->fieldSolution() );
    }
    if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::CauchyGreenInvariant1 ) )
    {
        this->M_exporter->step( time )->add( prefixvm(this->prefix(),"CauchyGreenInvariant1(SqrtTrCofC)"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"CauchyGreenInvariant1(SqrtTrCofC)")),
                                       *this->cauchyGreenInvariant1() );
    }
    if ( this->hasPostProcessFieldExported( LevelSetFieldsExported::CauchyGreenInvariant2 ) )
    {
        this->M_exporter->step( time )->add( prefixvm(this->prefix(),"CauchyGreenInvariant2(TrC_2SqrtTrCofC)"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"CauchyGreenInvariant2(TrC_2SqrtTrCofC)")),
                                       *this->cauchyGreenInvariant2() );
    }

    this->M_exporter->save();

    double tElapsed = this->timerTool("PostProcessing").stop("exportResults");
    this->log("LevelSet","exportResults", (boost::format("finish in %1% s")%tElapsed).str() );
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
    if( this->hasPostProcessMeasureExported( LevelSetMeasuresExported::Position_COM ) )
    {
        auto com = this->positionCOM();
        std::vector<double> vecCOM = { com(0,0) };
        if( nDim > 1 ) vecCOM.push_back( com(1,0) );
        if( nDim > 2 ) vecCOM.push_back( com(2,0) );
        this->postProcessMeasuresIO().setMeasureComp( "position_com", vecCOM );
        hasMeasureToExport = true;
    }
    if( this->hasPostProcessMeasureExported( LevelSetMeasuresExported::Velocity_COM ) )
    {
        auto ucom = this->velocityCOM();
        std::vector<double> vecUCOM = { ucom(0,0) };
        if( nDim > 1 ) vecUCOM.push_back( ucom(1,0) );
        if( nDim > 2 ) vecUCOM.push_back( ucom(2,0) );
        this->postProcessMeasuresIO().setMeasureComp( "velocity_com", vecUCOM );
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

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
auto
LEVELSET_CLASS_TEMPLATE_TYPE::positionCOM() const
{
    auto com = integrate( 
            _range=elements(this->mesh()), 
            _expr=vf::P() * (1.-idv(this->H()))
            ).evaluate();
    com = com / this->volume();

    return com;
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
auto
LEVELSET_CLASS_TEMPLATE_TYPE::velocityCOM() const
{
    auto ucom = integrate( 
            _range=elements(this->mesh()), 
            _expr=idv(M_advectionToolbox->fieldAdvectionVelocity()) * (1.-idv(this->H()))
            ).evaluate();
    ucom /= this->volume();

    return ucom;
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

    auto rangeElts = mesh->elementsWithProcessId( mesh->worldComm().localRank() );
    auto it_elt = std::get<0>( rangeElts );
    auto en_elt = std::get<1>( rangeElts );
    if (it_elt == en_elt) return;

    for (; it_elt!=en_elt; it_elt++)
    {
        auto const& elt = boost::unwrap_ref( *it_elt );
        int nbplus = 0;
        int nbminus = 0;

        for (int j=0; j<ndofv ; j++)
        {
            if (phi->localToGlobal(elt.id(), j, 0) >= 0.)
                nbplus++;
            else
                nbminus++;
        }

        //if elt crossed by interface
        if ( (nbminus != ndofv) && (nbplus!=ndofv) )
            M_markerInterface->assign(elt.id(), 0, 0, 1);
        else
            M_markerInterface->assign(elt.id(), 0, 0, 0);
    }
}

LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::updateMarkerDirac()
{
    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();

    auto rangeElts = mesh->elementsWithProcessId( mesh->worldComm().localRank() );
    auto it_elt = std::get<0>( rangeElts );
    auto en_elt = std::get<1>( rangeElts );
    if (it_elt == en_elt) return;

    double dirac_cut = this->dirac()->max() / 10.;

    for (; it_elt!=en_elt; it_elt++)
    {
        auto const& elt = boost::unwrap_ref( *it_elt );
        bool mark_elt = false;
        for (int j=0; j<ndofv; j++)
        {
            if ( std::abs( this->dirac()->localToGlobal(elt.id(), j, 0) ) > dirac_cut )
            {
                mark_elt = true;
                break; //don't need to do the others dof
            }
        }
        if( mark_elt )
            M_markerDirac->assign(elt.id(), 0, 0, 1);
        else
            M_markerDirac->assign(elt.id(), 0, 0, 0);
    }
}//markerDirac


LEVELSET_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSET_CLASS_TEMPLATE_TYPE::markerHeavisideImpl( element_markers_ptrtype const& marker, bool invert, double cut )
{
    /* returns P0 element having :
    if invert == true : 1 on elements inside Heaviside function (where H is smaller than epsilon on at least 1 dof)
    if invert == false : 1 on elements where Heaviside function greater than epsilon
    if cut_in_out == true : marker = 1 where H>0.5 and 0 where H<0.5
    */

    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();

    auto rangeElts = mesh->elementsWithProcessId( mesh->worldComm().localRank() );
    auto it_elt = std::get<0>( rangeElts );
    auto en_elt = std::get<1>( rangeElts );
    if (it_elt == en_elt) return;

    for (; it_elt!=en_elt; it_elt++)
    {
        auto const& elt = boost::unwrap_ref( *it_elt );
        bool mark_elt = false;
        for (int j=0; j<ndofv; j++)
        {
            if ( std::abs( this->heaviside()->localToGlobal(elt.id(), j, 0) ) > cut )
            {
                mark_elt = true;
                break;
            }
        }
        if( mark_elt )
            marker->assign(elt.id(), 0, 0, (invert)?0:1);
        else
            marker->assign(elt.id(), 0, 0, (invert)?1:0);
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

    auto rangeElts = mesh->elementsWithProcessId( mesh->worldComm().localRank() );
    auto it_elt = std::get<0>( rangeElts );
    auto en_elt = std::get<1>( rangeElts );
    if (it_elt == en_elt) return;

    auto phi = this->phi();
    auto phio = this->phio();

    auto prod = vf::project(this->functionSpace(), elements(mesh),
                            idv(phio) * idv(phi) );


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
