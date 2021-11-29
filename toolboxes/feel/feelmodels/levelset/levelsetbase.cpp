#include <feel/feelmodels/levelset/levelsetbase.hpp>

#include <feel/feelmodels/levelset/levelsetdeltaexpr.hpp>
#include <feel/feelmodels/levelset/levelsetheavisideexpr.hpp>

#include <boost/assign/list_of.hpp>
#include <boost/assign/list_inserter.hpp>

namespace Feel {

// Prevent LevelSetRedistanciation* instantiations
extern template class LevelSetRedistanciationFM< 
    typename FeelModels::LEVELSETBASE_CLASS_TYPE::space_levelset_type 
    >;
extern template class LevelSetRedistanciationHJ< 
    typename FeelModels::LEVELSETBASE_CLASS_TYPE::space_levelset_type 
    >;

namespace FeelModels {

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
LEVELSETBASE_CLASS_TEMPLATE_TYPE::LevelSetBase( 
        std::string const& prefix,
        std::string const& keyword,
        worldcomm_ptr_t const& worldComm,
        std::string const& subPrefix,
        ModelBaseRepository const& modelRep ) 
:
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep ),
    M_doUpdateDirac(true),
    M_doUpdateHeaviside(true),
    M_doUpdateInterfaceElements(true),
    M_doUpdateRangeDiracElements(true),
    M_doUpdateInterfaceFaces(true),
    M_doUpdateSmootherInterface(true),
    M_doUpdateSmootherInterfaceVectorial(true),
    M_doUpdateNormal(true),
    M_doUpdateCurvature(true),
    M_doUpdateGradPhi(true),
    M_doUpdateModGradPhi(true),
    M_doUpdatePhiPN(true),
    M_doUpdateDistance(true),
    M_doUpdateDistanceNormal(true),
    M_doUpdateDistanceCurvature(true),
    M_doUpdateSubmeshDirac(true),
    M_doUpdateSubmeshOuter(true),
    M_doUpdateSubmeshInner(true),
    M_doUpdateMarkers(true),
    M_buildExtendedDofSpace(false),
    M_useCurvatureDiffusion(false),
    //M_periodicity(periodicityLS),
    M_redistanciationIsUpdatedForUse(false),
    M_hasRedistanciated(false)
{
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

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::self_ptrtype
LEVELSETBASE_CLASS_TEMPLATE_TYPE::New(
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
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::createMesh()
{
    this->log("LevelSetBase","createMesh","start");
    this->timerTool("Constructor").start();

    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( this->mesh() ) << "mesh generation fail";

    double tElapsed = this->timerTool("Constructor").stop("create");
    this->log("LevelSetBase","createMesh", (boost::format("finish in %1% s") %tElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::init()
{
    this->log("LevelSetBase", "init", "start");

    // Mesh and space manager
    if( !M_spaceManager )
    {
        if( !this->mesh() )
        {
            // Create mesh
            this->createMesh();
        }
        M_spaceManager = std::make_shared<levelset_space_manager_type>( this->mesh() );
    }
    else
    {
        this->setMesh( M_spaceManager->mesh() );
    }
    // Tool manager
    if( !M_toolManager )
        M_toolManager = std::make_shared<levelset_tool_manager_type>( M_spaceManager, this->prefix() );
    // Function spaces
    this->functionSpaceManager()->setBuildExtendedDofTable( this->buildExtendedDofSpace() );
    this->createFunctionSpaces();
    // Tools
    this->createInterfaceQuantities();
    this->createTools();
    this->createRedistanciation();

    // Initial value
    //if( !this->doRestart() )
    //{
    // Set levelset initial value
    // Note: we initialise the levelset even in case of restart to
    // compute the initial geometrical quantities (volume, perimeter, ...)
    this->initLevelsetValue();
    //}

    // Init user-defined functions
    this->initUserFunctions();
    // Init post-process
    this->initPostProcess();

    this->log("LevelSetBase", "init", "finish");
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::initLevelsetValue()
{
    this->log("LevelSetBase", "initLevelsetValue", "start");

    if( !M_initialPhi ) // look for JSON initial values
    {
        bool hasInitialValue = !this->modelProperties().initialConditions().empty();

        if( hasInitialValue )
        {
            auto phiInit = this->functionSpace()->elementPtr();
            phiInit->setConstant( std::numeric_limits<value_type>::max() );
            std::vector<element_levelset_ptrtype> icLevelSetFields = { phiInit };

            this->modelProperties().parameters().updateParameterValues();
            auto paramValues = this->modelProperties().parameters().toParameterValues();
            this->modelProperties().initialConditions().setParameterValues( paramValues );

            this->updateInitialConditions( this->prefix(), this->rangeMeshElements(), this->symbolsExpr(), icLevelSetFields );

            ModelInitialConditionTimeSet const& icts = this->modelProperties().initialConditions().get( this->prefix() );
            for( auto const& [time,icByType]: icts )
            {
                auto itFindIcShapes = icByType.find( "Shapes" );
                if( itFindIcShapes != icByType.end() )
                {
                    for( auto const& icShape: itFindIcShapes->second )
                    {
                        this->addShape( icShape.pTree(), *phiInit );
                    }
                }
            }

            this->setInitialValue( phiInit );
        }
    }

    // Synchronize with current phi
    if( M_initialPhi ) // user-provided initial value
    {
        *M_phi = *M_initialPhi;
    }
    else // no initial value
    {
        M_phi->zero();
    }

    this->updateInterfaceQuantities();

    M_initialVolume = this->volume();
    M_initialPerimeter = this->perimeter();

    this->log("LevelSetBase", "initLevelsetValue", "finish");
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::addShape( 
        pt::ptree const& shapeParameters,
        element_levelset_type & phi 
        )
{
    boost::optional<std::string> shapeStr = shapeParameters.get_optional<std::string>( "shape" );
    CHECK( shapeStr ) << "missing shape type\n";
    auto const& shapeTypeIt = ShapeTypeMap.find( *shapeStr );
    if( shapeTypeIt != ShapeTypeMap.end() )
    {
        ShapeType shapeType = shapeTypeIt->second;

        switch(shapeType)
        {
            case ShapeType::SPHERE:
            {
                auto X = Px() - shapeParameters.get("xc", 0.);
                auto Y = Py() - shapeParameters.get("yc", 0.);
                auto Z = Pz() - shapeParameters.get("zc", 0.); 
                auto R = shapeParameters.get("radius", 0.);
                phi = vf::project(
                        _space=this->functionSpace(),
                        _range=elements(this->functionSpace()->mesh()),
                        _expr=vf::min( idv(phi), sqrt(X*X+Y*Y+Z*Z)-R ),
                        _geomap=this->geomap()
                        );
            }
            break;

            case ShapeType::ELLIPSE:
            {
                auto X = Px() - shapeParameters.get("xc", 0.);
                auto Y = Py() - shapeParameters.get("yc", 0.);
                auto Z = Pz() - shapeParameters.get("zc", 0.);
                double A = shapeParameters.get("a", 0.);
                double B = shapeParameters.get("b", 0.);
                double C = shapeParameters.get("c", 0.);
                double psi = shapeParameters.get("psi", 0.);
                double theta = shapeParameters.get("theta", 0.);
                // Apply inverse ZYX TaitBryan rotation
                double cosPsi = std::cos(psi); double sinPsi = std::sin(psi);
                double cosTheta = std::cos(theta); double sinTheta = std::sin(theta);
                auto Xp = cosTheta*(cosPsi*X+sinPsi*Y) + sinTheta*Z;
                auto Yp = -sinPsi*X + cosPsi*Y;
                auto Zp = -sinTheta*(cosPsi*X+sinPsi*Y) + cosTheta*Z;
                // Project
                phi = vf::project(
                        _space=this->functionSpace(),
                        _range=elements(this->functionSpace()->mesh()),
                        _expr=vf::min( idv(phi), sqrt(Xp*Xp+Yp*Yp*(A*A)/(B*B)+Zp*Zp*(A*A)/(C*C))-A ),
                        _geomap=this->geomap()
                        );
            }
            break;
        }
    }
    else
    {
        CHECK( false ) << "unknown shape type " << *shapeStr << std::endl;
    }
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::interfaceRectangularFunction( element_levelset_type const& p ) const
{
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
            _space=this->functionSpace(), 
            _range=this->rangeMeshElements(),
            _expr=R_expr
            );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::initUserFunctions()
{
    for ( auto const& modelfunc : this->modelProperties().functions() )
    {
        auto const& funcData = modelfunc.second;
        std::string funcName = funcData.name();

        if ( funcData.isScalar() )
        {
            if ( this->hasFieldUserScalar( funcName ) )
                continue;
            M_fieldsUserScalar[funcName] = this->functionSpaceScalar()->elementPtr();
        }
        else if ( funcData.isVectorial2() )
        {
            if ( nDim != 2 ) continue;
            if ( this->hasFieldUserVectorial( funcName ) )
                continue;
            M_fieldsUserVectorial[funcName] = this->functionSpaceVectorial()->elementPtr();
        }
        else if ( funcData.isVectorial3() )
        {
            if ( nDim != 3 ) continue;
            if ( this->hasFieldUserVectorial( funcName ) )
                continue;
            M_fieldsUserVectorial[funcName] = this->functionSpaceVectorial()->elementPtr();
        }
    }

    // update custom field given by registerCustomField
    for ( auto & [name,uptr] : M_fieldsUserScalar )
        if ( !uptr )
            uptr = this->functionSpaceScalar()->elementPtr();
    for ( auto & [name,uptr] : M_fieldsUserVectorial )
        if ( !uptr )
            uptr = this->functionSpaceVectorial()->elementPtr();

    this->updateUserFunctions();
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateUserFunctions( bool onlyExprWithTimeSymbol )
{
    if ( this->modelProperties().functions().empty() )
        return;

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().functions().setParameterValues( paramValues );
    for ( auto const& modelfunc : this->modelProperties().functions() )
    {
        auto const& funcData = modelfunc.second;
        if ( onlyExprWithTimeSymbol && !funcData.hasSymbol("t") )
            continue;

        std::string funcName = funcData.name();
        if ( funcData.isScalar() )
        {
            CHECK( this->hasFieldUserScalar( funcName ) ) << "user function " << funcName << "not registered";
            M_fieldsUserScalar[funcName]->on(_range=this->rangeMeshElements(),_expr=funcData.expressionScalar() );
        }
        else if ( funcData.isVectorial2() )
        {
            if constexpr ( nDim == 2 )
            {
                CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
                M_fieldsUserVectorial[funcName]->on(_range=this->rangeMeshElements(),_expr=funcData.expressionVectorial2() );
            }
        }
        else if ( funcData.isVectorial3() )
        {
            if constexpr ( nDim == 3 )
            {
                CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
                M_fieldsUserVectorial[funcName]->on(_range=this->rangeMeshElements(),_expr=funcData.expressionVectorial3() );
            }
        }
    }
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
std::set<std::string>
LEVELSETBASE_CLASS_TEMPLATE_TYPE::postProcessSaveAllFieldsAvailable() const
{
    std::set<std::string> postProcessSaveAllFieldsAvailable( { "phi", "dirac", "heaviside", "normal", "curvature", "gradphi", "modgradphi", "distance", "distance-normal", "distance-curvature" } );
    // Add possible user defined fields
    for( auto const& f : this->modelProperties().postProcess().exports( this->keyword() ).fields() )
    {
        if ( this->hasFieldUserScalar( f ) || this->hasFieldUserVectorial( f ) )
        {
            postProcessSaveAllFieldsAvailable.insert( f );
        }
    }
    return postProcessSaveAllFieldsAvailable;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
std::set<std::string>
LEVELSETBASE_CLASS_TEMPLATE_TYPE::postProcessExportsAllFieldsAvailable() const
{
    return this->postProcessSaveAllFieldsAvailable();
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    auto mfields = this->modelFields();
    this->exportResults( time, this->symbolsExpr(mfields), mfields, this->allMeasuresQuantities() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::initPostProcessExportsAndMeasures()
{
    // Update post-process expressions
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    this->setPostProcessExportsAllFieldsAvailable( this->postProcessExportsAllFieldsAvailable() );
    this->setPostProcessSaveAllFieldsAvailable( this->postProcessSaveAllFieldsAvailable() );
    this->setPostProcessExportsPidName( "pid" );
    super_type::initPostProcess();

    // Point measures
    auto geospace = std::make_shared<GeometricSpace<mesh_type>>( this->mesh() );
    auto geospaceWithNames = std::make_pair( std::set<std::string>({this->keyword()}), geospace );
    auto meshes = hana::make_tuple( geospaceWithNames );
    M_measurePointsEvaluation = std::make_shared<measure_points_evaluation_type>( meshes );
    for ( auto /*const*/& evalPoints : this->modelProperties().postProcess().measuresPoint( this->keyword() ) )
    {
        evalPoints.updateForUse( this->symbolsExpr() );
        M_measurePointsEvaluation->init( evalPoints );
    }
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->initPostProcessExportsAndMeasures();
    this->createPostProcessExporters();
    this->createPostProcessMeasures();
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::createFunctionSpaces()
{
    if( M_useSpaceIsoPN )
    {
        this->functionSpaceManager()->createFunctionSpaceIsoPN();
        M_spaceLevelset = this->functionSpaceManager()->functionSpaceScalarIsoPN();
        M_spaceVectorial = this->functionSpaceManager()->functionSpaceVectorialIsoPN();
        M_spaceMarkers = this->functionSpaceManager()->functionSpaceMarkersIsoPN();
    }
    else
    {
        this->functionSpaceManager()->createFunctionSpaceDefault();
        M_spaceLevelset = this->functionSpaceManager()->functionSpaceScalar();
        M_spaceVectorial = this->functionSpaceManager()->functionSpaceVectorial();
        M_spaceMarkers = this->functionSpaceManager()->functionSpaceMarkers();
    }

    M_rangeMeshElements = elements( M_spaceLevelset->mesh() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::createInterfaceQuantities()
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

    M_phi.reset( new element_levelset_type(this->functionSpace(), "Phi") );
    M_heaviside.reset( new element_levelset_type(this->functionSpace(), "Heaviside") );
    M_dirac.reset( new element_levelset_type(this->functionSpace(), "Dirac") );
    M_levelsetNormal.reset( new element_vectorial_type(this->functionSpaceVectorial(), "Normal") );
    M_levelsetCurvature.reset( new element_levelset_type(this->functionSpace(), "Curvature") );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::createRedistanciation()
{
    switch( M_redistanciationMethod )
    { 
        case LevelSetDistanceMethod::NONE :
        // Nothing to do - Remove warning
        break;
        case LevelSetDistanceMethod::FASTMARCHING :
        {
            if( !M_redistanciationFM )
                this->createRedistanciationFM();
        }
        break;
        case LevelSetDistanceMethod::HAMILTONJACOBI :
        {
            if( !M_redistanciationHJ )
                this->createRedistanciationHJ();
            
            double thickness_heaviside;
            if( Environment::vm( _name="thickness-heaviside", _prefix=prefixvm(this->prefix(), "redist-hj")).defaulted() )
            {
                thickness_heaviside =  M_thicknessInterface;
            }
            else
            {
                thickness_heaviside =  doption( _name="thickness-heaviside", _prefix=prefixvm(this->prefix(), "redist-hj") );
            }
            M_redistanciationHJ->setThicknessHeaviside( thickness_heaviside );
        }
        break;
        case LevelSetDistanceMethod::RENORMALISATION :
        // Nothing to do - Remove warning
        break;
    }

    M_redistanciationIsUpdatedForUse = true;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::createRedistanciationFM()
{
    M_redistanciationFM.reset( 
            new LevelSetRedistanciationFM<space_levelset_type>( 
                this->functionSpace(), prefixvm(this->prefix(), "redist-fm") 
                ) 
            );
    M_redistanciationFM->setProjectorL2( this->projectorL2Scalar() );
    M_redistanciationFM->setProjectorSM( this->projectorSMScalar() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::createRedistanciationHJ()
{
    M_redistanciationHJ.reset( 
            new LevelSetRedistanciationHJ<space_levelset_type>( this->functionSpace(), prefixvm(this->prefix(), "redist-hj") ) 
            );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::createTools()
{
    if( M_useSpaceIsoPN )
    {
        this->toolManager()->createProjectorL2IsoPN();
        M_projectorL2Scalar = this->toolManager()->projectorL2ScalarIsoPN();
        M_projectorL2Vectorial = this->toolManager()->projectorL2VectorialIsoPN();

        this->toolManager()->createProjectorSMIsoPN();
        M_projectorSMScalar = this->toolManager()->projectorSMScalarIsoPN();
        M_projectorSMVectorial = this->toolManager()->projectorSMVectorialIsoPN();
    }
    else
    {
        this->toolManager()->createProjectorL2Default();
        M_projectorL2Scalar = this->toolManager()->projectorL2Scalar();
        M_projectorL2Vectorial = this->toolManager()->projectorL2Vectorial();

        this->toolManager()->createProjectorSMDefault();
        M_projectorSMScalar = this->toolManager()->projectorSMScalar();
        M_projectorSMVectorial = this->toolManager()->projectorSMVectorial();
    }

    if( this->useCurvatureDiffusion() )
    {
        this->toolManager()->createCurvatureDiffusion();
    }
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::createPostProcessExporters()
{
    if ( !this->postProcessExportsFields().empty() )
    {
        std::string geoExportType = "static";//this->geoExportType();//change_coords_only, change, static
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

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::createPostProcessMeasures()
{
    // Start measures export
    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasures().restart( this->timeInitial() );
    }
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_PN_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::phiPN() const
{
    CHECK( M_useSpaceIsoPN ) << "use-space-iso-pn must be enabled to use phiPN \n";

    if( !M_levelsetPhiPN )
        M_levelsetPhiPN.reset( new element_levelset_PN_type(this->functionSpaceManager()->functionSpaceScalarPN(), "PhiPN") );

    if( M_doUpdatePhiPN )
       const_cast<self_type*>(this)->updatePhiPN(); 

    return M_levelsetPhiPN;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_vectorial_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::gradPhi( bool update ) const
{
    if( !M_levelsetGradPhi )
        M_levelsetGradPhi.reset( new element_vectorial_type(this->functionSpaceVectorial(), "GradPhi") );

    if( update )
        this->updateGradPhi();

    return M_levelsetGradPhi;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::modGradPhi( bool update ) const
{
    if( !M_levelsetModGradPhi )
        M_levelsetModGradPhi.reset( new element_levelset_type(this->functionSpace(), "ModGradPhi") );

    if( update )
        this->updateModGradPhi();

    return M_levelsetModGradPhi;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::distance( bool update ) const
{
    if( !M_distance )
        M_distance.reset( new element_levelset_type(this->functionSpace(), "Distance") );

    if( update )
        this->updateDistance();

    return M_distance;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::heaviside( bool update ) const
{
    if( !M_heaviside )
        M_heaviside.reset( new element_levelset_type(this->functionSpace(), "Heaviside") );

    if( update )
        this->updateHeaviside();

    return M_heaviside;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::levelset_delta_expr_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::diracExpr() const
{
    std::string deltaImpl = "generic";
    if ( this->M_useRegularPhi )
        deltaImpl = "renorm";
    return levelsetDelta(
            _element=this->phiElt(),
            _thickness=this->thicknessInterface(),
            _use_adaptive_thickness=this->M_useAdaptiveThicknessInterface,
            _impl=deltaImpl
            );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::dirac( bool update ) const
{
    if( !M_dirac )
        M_dirac.reset( new element_levelset_type(this->functionSpace(), "Dirac") );

    if ( update )
        this->updateDirac();

    return M_dirac;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_vectorial_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::normal( bool update ) const
{
    if( !M_levelsetNormal )
        M_levelsetNormal.reset( new element_vectorial_type(this->functionSpaceVectorial(), "Normal") );

    if( update )
        this->updateNormal();

    return M_levelsetNormal;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::curvature( bool update ) const
{
    if( !M_levelsetCurvature )
        M_levelsetCurvature.reset( new element_levelset_type(this->functionSpace(), "Curvature") );

    if( update )
        this->updateCurvature();

    return M_levelsetCurvature;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_vectorial_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::distanceNormal( bool update ) const
{
    if( !M_distanceNormal )
        M_distanceNormal.reset( new element_vectorial_type(this->functionSpaceVectorial(), "DistanceNormal") );

    if( update )
        this->updateDistanceNormal();

    return M_distanceNormal;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::distanceCurvature( bool update ) const
{
    if( !M_distanceCurvature )
        M_distanceCurvature.reset( new element_levelset_type(this->functionSpace(), "DistanceCurvature") );

    if( update )
        this->updateDistanceCurvature();

    return M_distanceCurvature;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    M_useRegularPhi = boption(_name=prefixvm(this->prefix(),"use-regularized-phi"));
    M_useHeavisideDiracNodalProj = boption(_name=prefixvm(this->prefix(),"h-d-nodal-proj"));

    std::string redistmethod = soption( _name="redist-method", _prefix=this->prefix() );
    CHECK( LevelSetDistanceMethodIdMap.count( redistmethod ) ) << redistmethod << " is not in the list of possible redistanciation methods\n";
    M_redistanciationMethod = LevelSetDistanceMethodIdMap.at( redistmethod );

    std::string distancemethod = soption( _name="distance-method", _prefix=this->prefix() );
    CHECK( LevelSetDistanceMethodIdMap.count( distancemethod ) ) << distancemethod << " is not in the list of possible redistanciation methods\n";
    M_distanceMethod = LevelSetDistanceMethodIdMap.at( distancemethod );

    M_redistInitialValue = boption( _name="redist-initial-value", _prefix=this->prefix() );

    const std::string gradPhiMethod = soption( _name="gradphi-method", _prefix=this->prefix() );
    CHECK(LevelSetDerivationMethodMap.left.count(gradPhiMethod)) << gradPhiMethod <<" is not in the list of possible gradphi derivation methods\n";
    M_gradPhiMethod = LevelSetDerivationMethodMap.left.at(gradPhiMethod);

    if( Environment::vm( _name="modgradphi-method", _prefix=this->prefix() ).defaulted() &&
        !Environment::vm( _name="gradphi-method", _prefix=this->prefix() ).defaulted() )
    {
        M_modGradPhiMethod = M_gradPhiMethod;
    }
    else
    {
        const std::string modGradPhiMethod = soption( _name="modgradphi-method", _prefix=this->prefix() );
        CHECK(LevelSetDerivationMethodMap.left.count(modGradPhiMethod)) << modGradPhiMethod <<" is not in the list of possible modgradphi derivation methods\n";
        M_modGradPhiMethod = LevelSetDerivationMethodMap.left.at(modGradPhiMethod);
    }

    const std::string curvatureMethod = soption( _name="curvature-method", _prefix=this->prefix() );
    CHECK(LevelSetCurvatureMethodMap.left.count(curvatureMethod)) << curvatureMethod <<" is not in the list of possible curvature methods\n";
    M_curvatureMethod = LevelSetCurvatureMethodMap.left.at(curvatureMethod);

    if( M_curvatureMethod == LevelSetCurvatureMethod::DIFFUSION_ORDER1 
            || M_curvatureMethod == LevelSetCurvatureMethod::DIFFUSION_ORDER2 )
        this->setUseCurvatureDiffusion( true );

    M_useSpaceIsoPN = boption( _name="use-space-iso-pn", _prefix=this->prefix() );
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Update levelset-dependent functions
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateGradPhi() const
{
    this->log("LevelSetBase", "updateGradPhi", "start");
    this->timerTool("UpdateInterfaceData").start();

    *M_levelsetGradPhi = this->grad( this->phiElt(), M_gradPhiMethod );
    M_doUpdateGradPhi = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSetBase", "updateGradPhi", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateModGradPhi() const
{
    this->log("LevelSetBase", "updateModGradPhi", "start");
    this->timerTool("UpdateInterfaceData").start();

    *M_levelsetModGradPhi = this->modGrad( this->phiElt(), M_modGradPhiMethod );
    M_doUpdateModGradPhi = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSetBase", "updateModGradPhi", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateDirac() const
{
    this->log("LevelSetBase", "updateDirac", "start");
    this->timerTool("UpdateInterfaceData").start();

    auto eps0 = this->thicknessInterface();

    //if( M_useAdaptiveThicknessInterface )
    //{
        //auto gradPhi = this->gradPhi();
        //auto gradPhiX = vf::project(
                //_space=this->functionSpace(),
                //_range=this->rangeMeshElements(),
                //_expr=idv(gradPhi->comp(Component::X))
                //);
        //auto gradPhiY = vf::project(
                //_space=this->functionSpace(),
                //_range=this->rangeMeshElements(),
                //_expr=idv(gradPhi->comp(Component::Y))
                //);
//#if FEELPP_DIM == 3
        //auto gradPhiZ = vf::project(
                //_space=this->functionSpace(),
                //_range=this->rangeMeshElements(),
                //_expr=idv(gradPhi->comp(Component::Z))
                //);
//#endif
        //auto eps_elt = this->functionSpace()->element();
        //eps_elt = vf::project(
                //_space=this->functionSpace(),
                //_range=this->rangeMeshElements(),
                //_expr=(vf::abs(idv(gradPhiX))+vf::abs(idv(gradPhiY))
//#if FEELPP_DIM == 3
                    //+ vf::abs(idv(gradPhiZ))
//#endif
                    //)*cst(eps0)/idv(this->modGradPhi())
                //);

        //auto eps = idv(eps_elt);
        //auto psi = idv(this->phi());

        //if ( M_useHeavisideDiracNodalProj )
            //*M_dirac = vf::project( this->functionSpace(), this->rangeMeshElements(),
                    //Feel::FeelModels::levelsetDelta(psi, eps) );
        //else
            //*M_dirac = M_projectorL2Scalar->project( Feel::FeelModels::levelsetDelta(psi, eps) );
    //}
    auto const& psi = this->phiElt();

    if ( M_useHeavisideDiracNodalProj )
        M_dirac->on(_range=this->rangeMeshElements(),_expr=this->diracExpr() );
    else
        *M_dirac = M_projectorL2Scalar->project( _expr=this->diracExpr() );

    M_doUpdateDirac = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSetBase", "updateDirac", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateHeaviside() const
{
    this->log("LevelSetBase", "updateHeaviside", "start");
    this->timerTool("UpdateInterfaceData").start();

    auto eps = this->thicknessInterface();

    if ( M_useRegularPhi )
    {
        auto psi = idv(this->phiElt()) / idv(this->modGradPhi());

        if ( M_useHeavisideDiracNodalProj )
            M_heaviside->on(_range=this->rangeMeshElements(),_expr=Feel::FeelModels::levelsetHeaviside(psi, cst(eps)));
        else
            *M_heaviside = M_projectorL2Scalar->project( _expr=Feel::FeelModels::levelsetHeaviside(psi, cst(eps)) );
    }
    else
    {
        auto psi = idv(this->phiElt());

        if ( M_useHeavisideDiracNodalProj )
            M_heaviside->on(_range=this->rangeMeshElements(),_expr=Feel::FeelModels::levelsetHeaviside(psi, cst(eps)) );
        else
            *M_heaviside = M_projectorL2Scalar->project( _expr=Feel::FeelModels::levelsetHeaviside(psi, cst(eps)) );
    }

    M_doUpdateHeaviside = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSetBase", "updateHeaviside", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updatePhiPN()
{
    this->log("LevelSetBase", "updatePhiPN", "start");
    this->timerTool("UpdateInterfaceData").start();

    auto const& phi = this->phiElt();
    this->functionSpaceManager()->opInterpolationScalarToPN()->apply( phi, *M_levelsetPhiPN );
    //*M_levelsetPhiPN = M_projectorL2P1PN->project(
            //_expr=idv(phi)
            //);

    M_doUpdatePhiPN = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSetBase", "updatePhiPN", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateNormal() const
{
    this->log("LevelSetBase", "updateNormal", "start");
    this->timerTool("UpdateInterfaceData").start();

    //auto const& phi = this->phiElt();
    //*M_levelsetNormal = M_projectorL2Vectorial->project( _expr=trans(gradv(phi)) / sqrt(gradv(phi) * trans(gradv(phi))) );
    auto gradPhi = this->gradPhi();
    M_levelsetNormal->on(_range=this->rangeMeshElements(), _expr=idv(gradPhi) / sqrt(trans(idv(gradPhi)) * idv(gradPhi)) );

    M_doUpdateNormal = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSetBase", "updateNormal", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateCurvature() const
{
    this->log("LevelSetBase", "updateCurvature", "start");
    this->timerTool("UpdateInterfaceData").start();

    switch( M_curvatureMethod )
    {
        case LevelSetCurvatureMethod::NODAL_PROJECTION:
        {
            this->log("LevelSetBase", "updateCurvature", "perform nodal projection");
            M_levelsetCurvature->on( _range=this->rangeMeshElements(), _expr=divv(this->normal()) );
        }
        break;
        case LevelSetCurvatureMethod::L2_PROJECTION:
        {
            this->log("LevelSetBase", "updateCurvature", "perform L2 projection");
            //*M_levelsetCurvature = this->projectorL2()->project( _expr=divv(this->normal()) );
            auto const& phi = this->phiElt();
            *M_levelsetCurvature = this->projectorL2()->derivate( gradv(phi) / sqrt(gradv(phi) * trans(gradv(phi))) );
        }
        break;
        case LevelSetCurvatureMethod::SMOOTH_PROJECTION:
        {
            this->log("LevelSetBase", "updateCurvature", "perform smooth projection");
            //*M_levelsetCurvature = this->smoother()->project( _expr=divv(this->normal()) );
            auto const& phi = this->phiElt();
            *M_levelsetCurvature = this->smoother()->derivate( gradv(phi) / sqrt(gradv(phi) * trans(gradv(phi))) );
        }
        break;
        //case LevelSetCurvatureMethod::PN_NODAL_PROJECTION:
        //{
            //this->log("LevelSetBase", "updateCurvature", "perform PN-nodal projection");
            //auto phiPN = this->phiPN();
            //auto normalPN = vf::project(
                    //_space=this->functionSpaceManager()->functionSpaceVectorialPN(),
                    //_range=this->functionSpaceManager()->rangeMeshPNElements(),
                    //_expr=trans(gradv(phiPN)) / sqrt(gradv(phiPN)*trans(gradv(phiPN)))
                    //);
            //auto curvaturePN = vf::project(
                    //_space=this->functionSpaceManager()->functionSpaceScalarPN(),
                    //_range=this->functionSpaceManager()->rangeMeshPNElements(),
                    //_expr=divv(normalPN)
                    //);

            //this->functionSpaceManager()->opInterpolationScalarFromPN()->apply( curvaturePN, *levelsetCurvature );
        //}
        break;
        case LevelSetCurvatureMethod::DIFFUSION_ORDER1:
        {
            this->log("LevelSetBase", "updateCurvature", "perform diffusion order1");
            *M_levelsetCurvature = this->toolManager()->curvatureDiffusion()->curvatureOrder1( this->distance() );
        }
        break;
        case LevelSetCurvatureMethod::DIFFUSION_ORDER2:
        {
            this->log("LevelSetBase", "updateCurvature", "perform diffusion order2");
            *M_levelsetCurvature = this->toolManager()->curvatureDiffusion()->curvatureOrder2( this->distance() );
        }
        break;
    }

    M_doUpdateCurvature = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSetBase", "updateCurvature", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateDistance() const
{
    this->log("LevelSetBase", "updateDistance", "start");
    this->timerTool("UpdateInterfaceData").start();

    *M_distance = this->redistanciate( this->phiElt(), M_distanceMethod );

    M_doUpdateDistance = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSetBase", "updateDistance", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateDistanceNormal() const
{
    this->log("LevelSetBase", "updateDistanceNormal", "start");
    this->timerTool("UpdateInterfaceData").start();

    auto const& phi = this->distance();
    auto const N_expr = trans(gradv(phi)) / sqrt( gradv(phi)*trans(gradv(phi)) );

    switch( M_gradPhiMethod )
    {
        case LevelSetDerivationMethod::NODAL_PROJECTION:
            this->log("LevelSetBase", "updateDistanceNormal", "perform nodal projection");
            M_distanceNormal->on( _range=this->rangeMeshElements(), _expr=N_expr );
            break;
        case LevelSetDerivationMethod::L2_PROJECTION:
            this->log("LevelSetBase", "updateDistanceNormal", "perform L2 projection");
            *M_distanceNormal = this->projectorL2Vectorial()->project( _expr=N_expr );
            break;
        case LevelSetDerivationMethod::SMOOTH_PROJECTION:
            this->log("LevelSetBase", "updateDistanceNormal", "perform smooth projection");
            *M_distanceNormal = this->smootherVectorial()->project( _expr=N_expr );
            break;
        //case LevelSetDerivationMethod::PN_NODAL_PROJECTION:
            //this->log("LevelSetBase", "updateDistanceNormal", "perform PN-nodal projection");
            //CHECK( false ) << "TODO: updateDistanceNormal with PN_NODAL_PROJECTION method\n";
            //auto phiPN = this->phiPN();
            //auto gradPhiPN = vf::project(
                    //_space=this->functionSpaceManager()->functionSpaceVectorialPN(),
                    //_range=this->functionSpaceManager()->rangeMeshPNElements(),
                    //_expr=trans(gradv(phiPN))
                    //);
            //this->functionSpaceManager()->opInterpolationVectorialFromPN()->apply( gradPhiPN, *distanceNormal );
            //break;
    }

    M_doUpdateDistanceNormal = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSetBase", "updateDistanceNormal", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateDistanceCurvature() const
{
    this->log("LevelSetBase", "updateDistanceCurvature", "start");
    this->timerTool("UpdateInterfaceData").start();

    switch( M_curvatureMethod )
    {
        case LevelSetCurvatureMethod::NODAL_PROJECTION:
        {
            this->log("LevelSetBase", "updateDistanceCurvature", "perform nodal projection");
            M_distanceCurvature->on( _range=this->rangeMeshElements(), _expr=divv(this->distanceNormal()) );
        }
        break;
        case LevelSetCurvatureMethod::L2_PROJECTION:
        {
            this->log("LevelSetBase", "updateDistanceCurvature", "perform L2 projection");
            //*M_distanceCurvature = this->projectorL2()->project( _expr=divv(this->distanceNormal()) );
            *M_distanceCurvature = this->projectorL2()->derivate( trans(idv(this->distanceNormal())) );
        }
        break;
        case LevelSetCurvatureMethod::SMOOTH_PROJECTION:
        {
            this->log("LevelSetBase", "updateDistanceCurvature", "perform smooth projection");
            *M_distanceCurvature = this->smoother()->project( _expr=divv(this->distanceNormal()) );
        }
        break;
        //case LevelSetCurvatureMethod::PN_NODAL_PROJECTION:
        //{
            //this->log("LevelSetBase", "updateDistanceCurvature", "perform PN-nodal projection");
            //CHECK( false ) << "TODO: updateDistanceCurvature with PN_NODAL_PROJECTION method\n";
            //auto phiPN = this->phiPN();
            //auto normalPN = vf::project(
                    //_space=this->functionSpaceManager()->functionSpaceVectorialPN(),
                    //_range=this->functionSpaceManager()->rangeMeshPNElements(),
                    //_expr=trans(gradv(phiPN)) / sqrt(gradv(phiPN)*trans(gradv(phiPN)))
                    //);
            //auto curvaturePN = vf::project(
                    //_space=this->functionSpaceManager()->functionSpaceScalarPN(),
                    //_range=this->functionSpaceManager()->rangeMeshPNElements(),
                    //_expr=divv(normalPN)
                    //);

            //this->functionSpaceManager()->opInterpolationScalarFromPN()->apply( curvaturePN, *distanceCurvature );
        //}
        //break;
        case LevelSetCurvatureMethod::DIFFUSION_ORDER1:
        {
            this->log("LevelSetBase", "updateDistanceCurvature", "perform diffusion order1");
            *M_distanceCurvature = this->toolManager()->curvatureDiffusion()->curvatureOrder1( this->distance() );
        }
        break;
        case LevelSetCurvatureMethod::DIFFUSION_ORDER2:
        {
            this->log("LevelSetBase", "updateDistanceCurvature", "perform diffusion order2");
            *M_distanceCurvature = this->toolManager()->curvatureDiffusion()->curvatureOrder2( this->distance() );
        }
        break;
    }

    M_doUpdateDistanceCurvature = false;

    double timeElapsed = this->timerTool("UpdateInterfaceData").stop();
    this->log("LevelSetBase", "updateDistanceCurvature", "finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::projector_levelset_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::smootherInterface() const
{
    if( !M_smootherInterface || M_doUpdateSmootherInterface )
    {
        auto const spaceInterface = self_type::space_levelset_type::New( 
                _mesh=this->mesh(),
                _range=this->rangeDiracElements(),
                _worldscomm=this->worldsComm()
                );
        M_smootherInterface = Feel::projector( 
                spaceInterface, spaceInterface,
                Feel::backend(_name=prefixvm(this->prefix(),"smoother"), _worldcomm=this->worldCommPtr(), _rebuild=true), 
                DIFF,
                this->mesh()->hAverage()*doption(_name="smooth-coeff", _prefix=prefixvm(this->prefix(),"smoother"))/Order,
                30);
        M_doUpdateSmootherInterface = false;
    }
    return M_smootherInterface;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::projector_levelset_vectorial_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::smootherInterfaceVectorial() const
{
    if( !M_smootherInterfaceVectorial || M_doUpdateSmootherInterfaceVectorial )
    {
        auto const spaceInterfaceVectorial = self_type::space_vectorial_type::New( 
                _mesh=this->mesh(),
                _range=this->rangeDiracElements(),
                _worldscomm=this->worldsComm()
                );
        M_smootherInterfaceVectorial = Feel::projector(
                spaceInterfaceVectorial, spaceInterfaceVectorial,
                Feel::backend(_name=prefixvm(this->prefix(),"smoother-vec"), _worldcomm=this->worldCommPtr(), _rebuild=true),
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
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::markerInterface() const
{
    if( !M_markerInterface )
        M_markerInterface.reset( new element_markers_type(M_spaceMarkers, "MarkerInterface") );

    if( M_doUpdateMarkers )
       const_cast<self_type*>(this)->updateMarkerInterface(); 

    return M_markerInterface;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::markerDirac() const
{
    if( !M_markerDirac )
        M_markerDirac.reset( new element_markers_type(M_spaceMarkers, "MarkerDirac") );

    if( M_doUpdateMarkers )
       const_cast<self_type*>(this)->updateMarkerDirac(); 

    return M_markerDirac;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::markerOuter( double cut ) const
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

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::markerInner( double cut ) const
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

//LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
//typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_markers_ptrtype const&
//LEVELSETBASE_CLASS_TEMPLATE_TYPE::markerCrossedElements() const
//{
    //if( !M_markerCrossedElements )
        //M_markerCrossedElements.reset( new element_markers_type(M_spaceMarkers, "MarkerCrossedElements") );

    //if( M_doUpdateMarkers )
       //const_cast<self_type*>(this)->updateMarkerCrossedElements(); 

    //return M_markerCrossedElements;
//}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::range_elements_type const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::rangeInterfaceElements() const
{
    if( this->M_doUpdateInterfaceElements )
    {
        M_interfaceElements = this->rangeInterfaceElementsImpl( this->phiElt() );
        M_doUpdateInterfaceElements = false;
    }

    return M_interfaceElements;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::range_elements_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::rangeInterfaceElementsImpl( element_levelset_type const& phi ) const
{
    mesh_ptrtype const& mesh = this->mesh();

    auto it_elt = mesh->beginOrderedElement();
    auto en_elt = mesh->endOrderedElement();

    const rank_type pid = mesh->worldCommElements().localRank();
    const int ndofv = space_levelset_type::fe_type::nDof;

    elements_reference_wrapper_ptrtype interfaceElts( new elements_reference_wrapper_type );

    for (; it_elt!=en_elt; it_elt++)
    {
        auto const& elt = boost::unwrap_ref( *it_elt );
        if ( elt.processId() != pid )
            continue;
        bool mark_elt = false;
        bool hasPositivePhi = false;
        bool hasNegativePhi = false;
        for (int j=0; j<ndofv; j++)
        {
            if ( phi.localToGlobal(elt.id(), j, 0) < 0. )
            {
                // phi < 0
                if( hasPositivePhi )
                {
                    mark_elt = true;
                    break; //don't need to do the others dof
                }
                hasNegativePhi = true;
            }
            if ( phi.localToGlobal(elt.id(), j, 0) > 0. )
            {
                // phi > 0
                if( hasNegativePhi )
                {
                    mark_elt = true;
                    break; //don't need to do the others dof
                }
                hasPositivePhi = true;
            }
        }
        if( mark_elt )
            interfaceElts->push_back( boost::cref(elt) );
    }

    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
            interfaceElts->begin(),
            interfaceElts->end(),
            interfaceElts
            );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::range_elements_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::rangeThickInterfaceElementsImpl( element_levelset_type const& phi, double thickness ) const
{
    mesh_ptrtype const& mesh = this->mesh();

    auto it_elt = mesh->beginOrderedElement();
    auto en_elt = mesh->endOrderedElement();

    const rank_type pid = mesh->worldCommElements().localRank();
    const int ndofv = space_levelset_type::fe_type::nDof;

    elements_reference_wrapper_ptrtype interfaceElts( new elements_reference_wrapper_type );

    for (; it_elt!=en_elt; it_elt++)
    {
        auto const& elt = boost::unwrap_ref( *it_elt );
        if ( elt.processId() != pid )
            continue;
        bool mark_elt = false;
        for (int j=0; j<ndofv; j++)
        {
            if ( std::abs(phi.localToGlobal(elt.id(), j, 0)) <= thickness )
            {
                mark_elt = true;
                break; //don't need to do the others dof
            }
        }
        if( mark_elt )
            interfaceElts->push_back( boost::cref(elt) );
    }

    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
            interfaceElts->begin(),
            interfaceElts->end(),
            interfaceElts
            );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::range_elements_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::rangeOuterElements( double cut ) const
{
    mesh_ptrtype const& mesh = this->mesh();

    element_levelset_type phi = this->functionSpace()->element();
    if( this->M_useRegularPhi )
        phi.on( _range=this->rangeMeshElements(), _expr=idv(this->phiElt()) / idv(this->modGradPhi()) );
    else
        phi = this->phiElt();

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

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::range_elements_type const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::rangeDiracElements() const
{
    if( M_doUpdateRangeDiracElements )
    {
        mesh_ptrtype const& mesh = this->mesh();

        auto it_elt = mesh->beginOrderedElement();
        auto en_elt = mesh->endOrderedElement();

        const rank_type pid = mesh->worldCommElements().localRank();
        const int ndofv = space_levelset_type::fe_type::nDof;

        double thickness = 2*this->thicknessInterface();
        elements_reference_wrapper_ptrtype diracElts( new elements_reference_wrapper_type );

        for (; it_elt!=en_elt; it_elt++)
        {
            auto const& elt = boost::unwrap_ref( *it_elt );
            if ( elt.processId() != pid )
                continue;
            bool mark_elt = false;
            for (int j=0; j<ndofv; j++)
            {
                if ( this->dirac()->localToGlobal(elt.id(), j, 0) > 0. )
                {
                    mark_elt = true;
                    break; //don't need to do the others dof
                }
            }
            if( mark_elt )
                diracElts->push_back( boost::cref(elt) );
        }

        M_rangeDiracElements = boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                diracElts->begin(),
                diracElts->end(),
                diracElts
                );

        M_doUpdateRangeDiracElements = false;
    }

    return M_rangeDiracElements;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::range_faces_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::rangeInterfaceFaces() const
{
    CHECK( Environment::isSequential() ) << "There is a bug to be fixed in parallel. Only run in sequential at the moment...\n";

    if( this->M_doUpdateInterfaceFaces )
    {
        mesh_ptrtype const& mesh = this->mesh();
        auto& markerIn = this->markerInner();
        markerIn->close();
        CHECK( markerIn->functionSpace()->extendedDofTable() ) << "interfaceFaces needs extended doftable in markers function space\n";

        auto it_face = mesh->beginFace();
        auto en_face = mesh->endFace();

        const rank_type pid = mesh->worldCommFaces().localRank();

        faces_reference_wrapper_ptrtype interfaceFaces( new faces_reference_wrapper_type );

        for (; it_face!=en_face; it_face++)
        {
            auto const& curFace = it_face->second; 
            if ( curFace.processId() != pid )
                continue;
            if ( !(curFace.isConnectedTo0() && curFace.isConnectedTo1()) )
                continue;
            bool isInnerElt0 = (markerIn->localToGlobal( curFace.element0().id(), 0, 0 ) > 1e-3);
            bool isInnerElt1 = (markerIn->localToGlobal( curFace.element1().id(), 0, 0 ) > 1e-3);

            if( (isInnerElt0 && !isInnerElt1) || (!isInnerElt0 && isInnerElt1) )
            {
                interfaceFaces->push_back( boost::cref(curFace) );
            }
        }

        M_interfaceFaces = boost::make_tuple( mpl::size_t<MESH_FACES>(),
                interfaceFaces->begin(),
                interfaceFaces->end(),
                interfaceFaces
                );

        M_doUpdateInterfaceFaces = false;
    }

    return M_interfaceFaces;
}

//----------------------------------------------------------------------------//
// Utility distances
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype
LEVELSETBASE_CLASS_TEMPLATE_TYPE::distToBoundary()
{
    element_levelset_ptrtype distToBoundary( new element_levelset_type(this->functionSpace(), "DistToBoundary") );

    // Retrieve the elements touching the boundary
    auto boundaryelts = boundaryelements( this->mesh() );

    // Init phi0 with h on marked2 elements
    auto phi0 = vf::project(
            _space=this->functionSpace(),
            _range=boundaryelts,
            _expr=h()
            );
    phi0.on( _range=boundaryfaces(this->mesh()), _expr=cst(0.) );

    // Run FM
    *distToBoundary = this->redistanciationFM()->run( phi0, boundaryelts );

    return distToBoundary;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype
LEVELSETBASE_CLASS_TEMPLATE_TYPE::distToMarkedFaces( boost::any const& marker )
{
    element_levelset_ptrtype distToMarkedFaces( new element_levelset_type(this->functionSpace(), "DistToMarkedFaces") );

    elements_reference_wrapper_ptrtype myelts( new elements_reference_wrapper_type );

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

    elements_reference_wrapper_t<mesh_type> myrange = boost::make_tuple(
            mpl::size_t<MESH_ELEMENTS>(), myelts->begin(), myelts->end(), myelts
            );

    // Init phi0 with h on marked2 elements
    auto phi0 = vf::project(
            _space=this->functionSpace(),
            _range=myrange,
            _expr=h()
            );
    phi0.on( _range=mfaces, _expr=cst(0.) );

    // Run FM using marker2 as marker DONE
    *distToMarkedFaces = this->redistanciationFM()->run( phi0, myrange );

    return distToMarkedFaces;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_ptrtype
LEVELSETBASE_CLASS_TEMPLATE_TYPE::distToMarkedFaces( std::initializer_list<boost::any> marker )
{
    element_levelset_ptrtype distToMarkedFaces( new element_levelset_type(this->functionSpace(), "DistToMarkedFaces") );

    elements_reference_wrapper_ptrtype myelts( new elements_reference_wrapper_type );

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

    elements_reference_wrapper_t<mesh_type> myrange = boost::make_tuple(
            mpl::size_t<MESH_ELEMENTS>(), myelts->begin(), myelts->end(), myelts
            );

    // Init phi0 with h on marked2 elements
    auto phi0 = vf::project(
            _space=this->functionSpace(),
            _range=myrange,
            _expr=h()
            );
    phi0.on( _range=mfaces, _expr=cst(0.) );

    // Run FM using marker2 as marker DONE
    *distToMarkedFaces = this->redistanciationFM()->run( phi0, myrange );

    return distToMarkedFaces;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateInterfaceQuantities()
{
    M_doUpdateDirac = true;
    M_doUpdateHeaviside = true;
    M_doUpdateInterfaceElements = true;
    M_doUpdateRangeDiracElements = true;
    M_doUpdateInterfaceFaces = true;
    M_doUpdateSmootherInterface = true;
    M_doUpdateSmootherInterfaceVectorial = true;
    M_doUpdateNormal = true;
    M_doUpdateCurvature = true;
    M_doUpdateMarkers = true;
    M_doUpdateGradPhi = true;
    M_doUpdateModGradPhi = true;
    M_doUpdatePhiPN = true;
    M_doUpdateDistance = true;
    M_doUpdateDistanceNormal = true;
    M_doUpdateDistanceCurvature = true;
    M_doUpdateSubmeshDirac = true;
    M_doUpdateSubmeshOuter = true;
    M_doUpdateSubmeshInner = true;
}

//----------------------------------------------------------------------------//
// Interface quantities helpers
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_vectorial_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::grad( element_levelset_type const& phi, LevelSetDerivationMethod method ) const
{
    this->log("LevelSetBase", "grad", "perform " + LevelSetDerivationMethodMap.right.at( method ) );
    return this->toolManager()->grad( phi, method );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::modGrad( element_levelset_type const& phi, LevelSetDerivationMethod method ) const
{
    this->log("LevelSetBase", "modGrad", "perform " + LevelSetDerivationMethodMap.right.at( method ) );
    return this->toolManager()->modGrad( phi, method );
}
#if 0
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_vectorial_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::grad( element_levelset_type const& phi, LevelSetDerivationMethod method ) const
{
    switch( method )
    {
        case LevelSetDerivationMethod::NODAL_PROJECTION:
            this->log("LevelSetBase", "grad", "perform nodal projection");
            return vf::project( 
                    _space=this->functionSpaceVectorial(),
                    _range=this->rangeMeshElements(),
                    _expr=trans(gradv(phi))
                    );
        case LevelSetDerivationMethod::L2_PROJECTION:
            this->log("LevelSetBase", "grad", "perform L2 projection");
            //return this->projectorL2Vectorial()->project( _expr=trans(gradv(phi)) );
            return this->projectorL2Vectorial()->derivate( idv(phi) );
        case LevelSetDerivationMethod::SMOOTH_PROJECTION:
            this->log("LevelSetBase", "grad", "perform smooth projection");
            return this->smootherVectorial()->project( _expr=trans(gradv(phi)) );
        case LevelSetDerivationMethod::PN_NODAL_PROJECTION:
            this->log("LevelSetBase", "grad", "perform PN-nodal projection");
            CHECK( M_useSpaceIsoPN ) << "use-space-iso-pn must be enabled to use PN_NODAL_PROJECTION \n";
            auto phiPN = this->functionSpaceManager()->opInterpolationScalarToPN()->operator()( phi );
            auto gradPhiPN = vf::project(
                    _space=this->functionSpaceManager()->functionSpaceVectorialPN(),
                    _range=this->functionSpaceManager()->rangeMeshPNElements(),
                    _expr=trans(gradv(phiPN))
                    );
            return this->functionSpaceManager()->opInterpolationVectorialFromPN()->operator()( gradPhiPN );
    }
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::modGrad( element_levelset_type const& phi, LevelSetDerivationMethod method ) const
{
    switch( method )
    {
        case LevelSetDerivationMethod::NODAL_PROJECTION:
            this->log("LevelSetBase", "modGrad", "perform nodal projection");
            return vf::project( 
                    _space=this->functionSpace(),
                    _range=this->rangeMeshElements(),
                    _expr=sqrt( gradv(phi)*trans(gradv(phi)) )
                    );
        case LevelSetDerivationMethod::L2_PROJECTION:
            this->log("LevelSetBase", "modGrad", "perform L2 projection");
            return this->projectorL2()->project( _expr=sqrt( gradv(phi)*trans(gradv(phi)) ) );
        case LevelSetDerivationMethod::SMOOTH_PROJECTION:
            this->log("LevelSetBase", "modGrad", "perform smooth projection");
            return this->smoother()->project( _expr=sqrt( gradv(phi)*trans(gradv(phi)) ) );
        case LevelSetDerivationMethod::PN_NODAL_PROJECTION:
            this->log("LevelSetBase", "modGrad", "perform PN-nodal projection");
            CHECK( M_useSpaceIsoPN ) << "use-space-iso-pn must be enabled to use PN_NODAL_PROJECTION \n";
            auto phiPN = this->functionSpaceManager()->opInterpolationScalarToPN()->operator()( phi );
            auto modGradPhiPN = vf::project(
                    _space=this->functionSpaceManager()->functionSpaceScalarPN(),
                    _range=this->functionSpaceManager()->rangeMeshPNElements(),
                    _expr=sqrt( gradv(phiPN)*trans(gradv(phiPN)) )
                    );
            return this->functionSpaceManager()->opInterpolationScalarFromPN()->operator()( modGradPhiPN );
    }
}
#endif

//----------------------------------------------------------------------------//
// Redistanciation
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::redistanciate()
{ 
    this->log("LevelSetBase", "redistanciate", "start");
    this->timerTool("Redist").start();

    auto phi = this->phiPtr();

    *phi = this->redistanciate( *phi, M_redistanciationMethod );

    M_hasRedistanciated = true;

    double timeElapsed = this->timerTool("Redist").stop();
    this->log("LevelSetBase","redistanciate","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::element_levelset_type
LEVELSETBASE_CLASS_TEMPLATE_TYPE::redistanciate( element_levelset_type const& phi, LevelSetDistanceMethod method ) const
{ 
    auto phiRedist = this->functionSpace()->elementPtr();

    switch( method )
    {
        case LevelSetDistanceMethod::NONE:
        {
            this->log("LevelSetBase", "redistanciate", "no redistanciation");
            *phiRedist = phi;
        }
        break;
        case LevelSetDistanceMethod::FASTMARCHING:
        {
            this->log("LevelSetBase", "redistanciate", "run fast-marching");
            *phiRedist = this->redistanciationFM()->run( phi );
        } // Fast Marching
        break;

        case LevelSetDistanceMethod::HAMILTONJACOBI:
        {
            // TODO
            this->log("LevelSetBase", "redistanciate", "run hamilton-jacobi");
            *phiRedist = this->redistanciationHJ()->run( phi );
        } // Hamilton-Jacobi
        break;

        case LevelSetDistanceMethod::RENORMALISATION:
        {
            //auto R = this->interfaceRectangularFunction( phi );
            this->log("LevelSetBase", "redistanciate", "renormalise");
            auto const modGradPhi = this->modGrad( phi );
            phiRedist->on(
                    _range=this->rangeMeshElements(),
                    //_expr = idv(phi) * ( 1. + ( 1./idv(modGradPhi)-1. )*idv(R) )
                    _expr = idv(phi) / idv(modGradPhi)
                    );
        }
        break;
    }
    
    return *phiRedist;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::redistanciationFM_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::redistanciationFM( bool buildOnTheFly ) const
{
    if( !M_redistanciationFM && buildOnTheFly )
    {
        const_cast<self_type*>(this)->createRedistanciationFM();
    }

    return M_redistanciationFM;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::redistanciationHJ_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::redistanciationHJ( bool buildOnTheFly ) const
{
    if( !M_redistanciationHJ && buildOnTheFly )
    {
        const_cast<self_type*>(this)->createRedistanciationHJ();
    }

    return M_redistanciationHJ;
}

//----------------------------------------------------------------------------//
// Initial value
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::setInitialValue(element_levelset_ptrtype const& phiv, bool doRedistanciate)
{
    this->log("LevelSetBase", "setInitialValue", "start");

    if( !M_initialPhi )
        M_initialPhi.reset( new element_levelset_type(this->functionSpace(), "InitialPhi") );

    if ( doRedistanciate )
    {
        this->log("LevelSetBase", "setInitialValue", "redistanciate");
        *M_initialPhi = this->redistanciate( *phiv, M_redistanciationMethod );
    }
    else
    {
        *M_initialPhi = *phiv;
    }

    this->log("LevelSetBase", "setInitialValue", "finish");
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
LEVELSETBASE_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::string advectionStabilization = soption( _name="stabilization.method", _prefix=this->prefix() );

    std::string hdProjectionMethod = (this->M_useHeavisideDiracNodalProj)? "nodal": "L2";

    const std::string gradPhiMethod = LevelSetDerivationMethodMap.right.at(this->M_gradPhiMethod);
    const std::string modGradPhiMethod = LevelSetDerivationMethodMap.right.at(this->M_modGradPhiMethod);
    const std::string curvatureMethod = LevelSetCurvatureMethodMap.right.at(this->M_curvatureMethod);

    std::string redistMethod;
    std::string redistmethod = soption( _name="redist-method", _prefix=this->prefix() );
    if( redistmethod == "fm" )
    {
        redistMethod = "Fast-Marching";
        std::string fmInitMethod = redistanciationFM_type::FastMarchingInitialisationMethodMap.right.at( this->redistanciationFM( false )->fastMarchingInitialisationMethod() );
        redistMethod += " (" + fmInitMethod + ")";
    }
    else if( redistmethod == "hj" )
        redistMethod = "Hamilton-Jacobi";

    std::string scalarSmootherParameters;
    if ( M_projectorSMScalar )
    {
        double scalarSmootherCoeff = this->smoother()->epsilon() * Order / this->mesh()->hAverage();
        scalarSmootherParameters = "coeff (h*c/order) = " 
            + std::to_string(this->smoother()->epsilon())
            + " (" + std::to_string(this->mesh()->hAverage()) + " * " + std::to_string(scalarSmootherCoeff) + " / " + std::to_string(Order) + ")"
            ;
    }
    std::string vectorialSmootherParameters;

    if ( M_projectorSMVectorial )
    {
        double vectorialSmootherCoeff = this->smootherVectorial()->epsilon() * Order / this->mesh()->hAverage();
        vectorialSmootherParameters = "coeff (h*c/order) = " 
            + std::to_string(this->smootherVectorial()->epsilon())
            + " (" + std::to_string(this->mesh()->hAverage()) + " * " + std::to_string(vectorialSmootherCoeff) + " / " + std::to_string(Order) + ")"
            ;
    }

    std::string restartMode = (this->doRestart())? "ON": "OFF";

    std::string exporterType, hovisuMode;
    int exporterFreq;
    hovisuMode = "OFF";
    if( M_exporter )
    {
        exporterType = this->M_exporter->type();
        exporterFreq = this->M_exporter->freq();
    }

    std::string exportedFields;
    for( std::string const& field : this->postProcessExportsFields() )
        exportedFields = (exportedFields.empty())? field : exportedFields + " - " + field;

    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||------------Info : LevelSet Base--------------||"
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
           << "\n     -- redistanciate initial value         : " << std::boolalpha << this->M_redistInitialValue
           << "\n     -- gradphi projection                  : " << gradPhiMethod
           << "\n     -- modgradphi projection               : " << modGradPhiMethod
           << "\n     -- curvature projection                : " << curvatureMethod

           << "\n   Redistanciation Parameters"
           << "\n     -- redistanciation method          : " << redistMethod;

    if( M_projectorSMScalar || M_projectorSMVectorial )
    *_ostr << "\n   Smoothers Parameters";
    if( M_projectorSMScalar )
    *_ostr << "\n     -- scalar smoother    : " << scalarSmootherParameters;
    if( M_projectorSMVectorial )
    *_ostr << "\n     -- vectorial smoother : " << vectorialSmootherParameters;

    *_ostr << "\n   Space Discretization";
#if 0
    if( this->hasGeoFile() )
    *_ostr << "\n     -- geo file name   : " << this->geoFile();
    *_ostr << "\n     -- mesh file name  : " << this->meshFile();
#endif
    *_ostr << "\n     -- nb elt in mesh  : " << this->mesh()->numGlobalElements()//numElements()
         //<< "\n     -- nb elt in mesh  : " << this->mesh()->numElements()
         //<< "\n     -- nb face in mesh : " << this->mesh()->numFaces()
           << "\n     -- hMin            : " << this->mesh()->hMin()
           << "\n     -- hMax            : " << this->mesh()->hMax()
           << "\n     -- hAverage        : " << this->mesh()->hAverage()
           << "\n     -- geometry order  : " << nOrderGeo
           << "\n     -- level set order : " << Order
           << "\n     -- nb dof          : " << this->functionSpace()->nDof() << " (" << this->functionSpace()->nLocalDof() << ")"
           ;

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
           ;

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
           << "\n";

    return _ostr;
}

//----------------------------------------------------------------------------//
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::mesh_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::submeshDirac() const
{
    if( M_doUpdateSubmeshDirac )
    {
        this->mesh()->updateMarker2( *this->markerDirac() );
        M_submeshDirac = createSubmesh( _mesh=this->mesh(), _range=marked2elements( this->mesh(), 1 ) );
        M_doUpdateSubmeshDirac = false;
    }
    return M_submeshDirac;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::mesh_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::submeshOuter( double cut ) const
{
    if( M_doUpdateSubmeshOuter || cut != M_markerOuterCut )
    {
        this->mesh()->updateMarker2( *this->markerOuter( cut ) );
        M_submeshOuter = createSubmesh( _mesh=this->mesh(), _range=marked2elements( this->mesh(), 1 ) );
        M_doUpdateSubmeshOuter = false;
    }
    return M_submeshOuter;
}

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETBASE_CLASS_TEMPLATE_TYPE::mesh_ptrtype const&
LEVELSETBASE_CLASS_TEMPLATE_TYPE::submeshInner( double cut ) const
{
    if( M_doUpdateSubmeshInner || cut != M_markerInnerCut )
    {
        this->mesh()->updateMarker2( *this->markerInner( cut ) );
        M_submeshInner = createSubmesh( _mesh=this->mesh(), _range=marked2elements( this->mesh(), 1 ) );
        M_doUpdateSubmeshInner = false;
    }
    return M_submeshInner;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Update markers
LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateMarkerInterface()
{
    /* returns a marker (P0_type) on the elements crossed by the levelset
       ie :
       - if the element has not phi on all the dof of the same sign -> the element is mark as 1
       - if phi on the dof of the element are of the same sign -> mark as 0
      */

    const int ndofv = space_levelset_type::fe_type::nDof;

    mesh_ptrtype const& mesh = this->mesh();
    auto phi = this->phiPtr();

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

LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::updateMarkerDirac()
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


LEVELSETBASE_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETBASE_CLASS_TEMPLATE_TYPE::markerHeavisideImpl( element_markers_ptrtype const& marker, bool invert, double cut )
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

} // FeelModels
} // Feel
