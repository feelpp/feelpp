/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelfilters/savegmshmesh.hpp>
//#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
//#include <feel/feelvf/inv.hpp>
#include <feel/feelpde/operatorpcd.hpp>
#include <feel/feells/reinit_fms.hpp>

#include <feel/feelmodels/modelmesh/markedmeshtool.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>


namespace Feel {
namespace FeelModels {

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::FluidMechanics( std::string const& prefix, std::string const& keyword,
                                                    worldcomm_ptr_t const& worldComm,
                                                    std::string const& subPrefix,
                                                    ModelBaseRepository const& modelRep )
    :
    super_type( prefix,keyword,worldComm,subPrefix, modelRep ),
    ModelBase( prefix,keyword,worldComm,subPrefix, modelRep ),
    ModelPhysics<nDim>( "fluid" ),
    M_applyMovingMeshBeforeSolve( true )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".FluidMechanicsConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".FluidMechanicsSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".FluidMechanicsPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".FluidMechanicsTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::self_ptrtype
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::New( std::string const& prefix, std::string const& keyword,
                                         worldcomm_ptr_t const& worldComm, std::string const& subPrefix,
                                         ModelBaseRepository const& modelRep )
{
    return std::make_shared<self_type>( prefix, keyword, worldComm, subPrefix, modelRep );

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::string
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::expandStringFromSpec( std::string const& expr )
{
    std::string res = expr;
    boost::replace_all( res, "$fluid_u_order", (boost::format("%1%")%nOrderVelocity).str() );
    boost::replace_all( res, "$fluid_p_order", (boost::format("%1%")%nOrderPressure).str() );
    boost::replace_all( res, "$fluid_geo_order", (boost::format("%1%")%nOrderGeo).str() );
    std::string fluidTag = (boost::format("P%1%P%2%G%3%")%nOrderVelocity %nOrderPressure %nOrderGeo ).str();
    boost::replace_all( res, "$fluid_tag", fluidTag );
    return res;
}

// add members instatantiations need by static function expandStringFromSpec
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
const uint16_type FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::nOrderVelocity;
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
const uint16_type FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::nOrderPressure;
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
const uint16_type FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::nOrderGeo;



//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::loadMesh( mesh_ptrtype __mesh )
{
    this->log("FluidMechanics","loadMesh", "start");
    //-----------------------------------------------------------------------------//
    // create or reload mesh
    if (this->doRestart() && !__mesh)
        this->initMesh();
    else
        this->setMesh(  __mesh );
    //-----------------------------------------------------------------------------//
    this->log("FluidMechanics","loadMesh", "finish");
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    this->log("FluidMechanics","loadParameterFromOptionsVm", "start");

    //--------------------------------------------------------------//
    // exporters options
    M_isHOVisu = nOrderGeo > 2;
    if ( Environment::vm().count(prefixvm(this->prefix(),"hovisu").c_str()) )
        M_isHOVisu = boption(_name="hovisu",_prefix=this->prefix());
    //--------------------------------------------------------------//
    M_haveSourceAdded=false;//true when update
    M_velocityDivIsEqualToZero=true;

    M_dirichletBCnitscheGamma = doption(_name="dirichletbc.nitsche.gamma",_prefix=this->prefix());

    M_useSemiImplicitTimeScheme = boption(_name="use-semi-implicit-time-scheme",_prefix=this->prefix());
    std::string _solver = soption(_name="solver",_prefix=this->prefix());
    if ( _solver != "automatic" )
        this->setSolverName( _solver );
    else
        M_solverName = _solver;
    M_useVelocityExtrapolated = M_useSemiImplicitTimeScheme;

    //--------------------------------------------------------------//
    // fsi options
    M_useFSISemiImplicitScheme = false;
    M_couplingFSIcondition = "dirichlet-neumann";

    //--------------------------------------------------------------//
    // time stepping
    M_timeStepping = soption(_name="time-stepping",_prefix=this->prefix());
    M_timeStepThetaValue = doption(_name="time-stepping.theta.value",_prefix=this->prefix());
    //--------------------------------------------------------------//
    // start solver options
    M_startBySolveNewtonian = boption(_prefix=this->prefix(),_name="start-by-solve-newtonian");
    M_hasSolveNewtonianAtKickOff = false;
    M_startBySolveStokesStationary = boption(_prefix=this->prefix(),_name="start-by-solve-stokes-stationary");
    M_hasSolveStokesStationaryAtKickOff = false;

    //--------------------------------------------------------------//
    // stabilisation options
    M_stabilizationGLS = boption(_name="stabilization-gls",_prefix=this->prefix());
    M_stabilizationGLSType = soption(_name="stabilization-gls.type",_prefix=this->prefix());
    M_stabilizationGLSDoAssembly = true;
    M_stabilizationGLS_checkViscosityDependencyOnCoordinates = boption(_name="stabilization-gls.check-viscosity-dependency-on-coordinates",_prefix=this->prefix());

    M_applyCIPStabOnlyOnBoundaryFaces=false;
    M_doCIPStabConvection = boption(_name="stabilisation-cip-convection",_prefix=this->prefix());
    M_doCIPStabDivergence = boption(_name="stabilisation-cip-divergence",_prefix=this->prefix());
    M_doCIPStabPressure = boption(_name="stabilisation-cip-pressure",_prefix=this->prefix());
    M_stabCIPConvectionGamma = doption(_name="stabilisation-cip-convection-gamma",_prefix=this->prefix());
    M_stabCIPDivergenceGamma = doption(_name="stabilisation-cip-divergence-gamma",_prefix=this->prefix());
    M_stabCIPPressureGamma = doption(_name="stabilisation-cip-pressure-gamma",_prefix=this->prefix());

    // M_doStabDivDiv = boption(_name="stabilisation-div-div",_prefix=this->prefix());
    //M_doCstPressureStab = boption(_name="stabilisation-cstpressure",_prefix=this->prefix());
    M_doStabConvectionEnergy = boption(_name="stabilisation-convection-energy",_prefix=this->prefix());

    M_definePressureCst = boption(_name="define-pressure-cst",_prefix=this->prefix());
    M_definePressureCstMethod = soption(_name="define-pressure-cst.method",_prefix=this->prefix());
    CHECK( M_definePressureCstMethod == "lagrange-multiplier" || M_definePressureCstMethod == "penalisation" ||
           M_definePressureCstMethod == "algebraic" ) << "lagrange-multiplier or penalisation or algebraic";
    M_definePressureCstPenalisationBeta = doption(_name="define-pressure-cst.penalisation-beta",_prefix=this->prefix());
    M_definePressureCstMarkers.clear();
    if ( Environment::vm().count( prefixvm(this->prefix(),"define-pressure-cst.markers").c_str() ) )
    {
        std::vector<std::string> inputMarkers = Environment::vm()[ prefixvm(this->prefix(),"define-pressure-cst.markers").c_str() ].template as<std::vector<std::string> >();
        std::string inputMarkersAsString;
        for ( std::string const& marker : inputMarkers )
            inputMarkersAsString += marker;

        boost::char_separator<char> sep(",");
        boost::char_separator<char> sep2(":");
        boost::tokenizer< boost::char_separator<char> > kvlist( inputMarkersAsString, sep );
        for( const auto& ikvl : kvlist )
        {
            boost::tokenizer< boost::char_separator<char> > kvlist2( ikvl, sep2);
            std::set<std::string> markerList;
            for( const auto& ikvl2 : kvlist2 )
                markerList.insert( ikvl2 );

            if ( !markerList.empty() )
                M_definePressureCstMarkers.push_back( markerList );
        }
    }

    M_dist2WallEnabled = boption(_prefix=this->prefix(),_name="distance-to-wall.enabled");
    if ( Environment::vm().count( prefixvm(this->prefix(),"distance-to-wall.markers").c_str() ) )
    {
        std::vector<std::string> tmp =  Environment::vm()[ prefixvm(this->prefix(),"distance-to-wall.markers").c_str() ].template as<std::vector<std::string> >();
        M_dist2WallMarkers.insert( tmp.begin(), tmp.end() );
    }

    M_useSemiImplicitTurbulenceCoupling = boption(_prefix=this->prefix(),_name="use-semi-implicit-turbulence-coupling");

    // prec
    M_preconditionerAttachPMM = boption(_prefix=this->prefix(),_name="preconditioner.attach-pmm");
    M_pmmNeedUpdate = false;
    M_preconditionerAttachPCD = boption(_prefix=this->prefix(),_name="preconditioner.attach-pcd");

    this->log("FluidMechanics","loadParameterFromOptionsVm", "finish");
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("FluidMechanics","initMesh", "start");
    this->timerTool("Constructor").start();

    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( this->mesh() ) << "mesh generation fail";

    double tElapsed = this->timerTool("Constructor").stop("initMesh");
    this->log("FluidMechanics","initMesh", (boost::format("finish in %1% s") %tElapsed).str() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initMaterialProperties()
{
    this->log("FluidMechanics","initMesh", "start");
    this->timerTool("Constructor").start();

    // auto paramValues = this->modelProperties().parameters().toParameterValues();
    // this->modelProperties().materials().setParameterValues( paramValues );
    if ( !M_materialsProperties )
    {
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    // modif expression of a material properties with non newtonian cases
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
        if ( physicFluidData->dynamicViscosity().isNewtonianLaw() )
            continue;
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            std::string _viscositySymbol = (boost::format("%1%_%2%_mu")%this->keyword() %matName).str();
            std::string newMuExprStr = (boost::format("%1%:%1%")%_viscositySymbol ).str();
            ModelExpression newMuExpr;
            newMuExpr.setExpr( newMuExprStr, this->worldComm(), this->repository().expr() );
            this->materialsProperties()->addProperty( this->materialsProperties()->materialProperties( matName ), "dynamic-viscosity", newMuExpr, true );
        }
    }

    double tElapsed = this->timerTool("Constructor").stop("initMesh");
    this->log("FluidMechanics","initMesh", (boost::format("finish in %1% s") %tElapsed).str() );
}
//---------------------------------------------------------------------------------------------------------//
#if 0
namespace detail
{
template <typename FMtype>
typename FMtype::space_fluid_ptrtype
createFluidFunctionSpaces( FMtype const& FM, std::vector<bool> const& extendedDT, mpl::false_)
{
    return FMtype::space_fluid_type::New( _mesh=FM.mesh(), _worldscomm=FM.worldsComm(),
                                          _extended_doftable=extendedDT );
}
template <typename FMtype>
typename FMtype::space_fluid_ptrtype
createFluidFunctionSpaces( FMtype const& FM, std::vector<bool> const& extendedDT, mpl::true_)
{
    node_type translat( FMtype::nDim );
    translat[0] = doption(_name="periodicity.translate-x",_prefix=FM.prefix());
    if ( FMtype::nDim >=2 )
        translat[1] = doption(_name="periodicity.translate-y",_prefix=FM.prefix());
    if ( FMtype::nDim == 3 )
        translat[2]= doption(_name="periodicity.translate-z",_prefix=FM.prefix());
    std::string marker1 = soption(_name="periodicity.marker1",_prefix=FM.prefix());
    std::string marker2 = soption(_name="periodicity.marker2",_prefix=FM.prefix());
    auto theperiodicity = periodicity( Periodic<>( FM.mesh()->markerName(marker1),FM.mesh()->markerName(marker2), translat), NoPeriodicity() );
    return FMtype::space_fluid_type::New( _mesh=FM.mesh(), _worldscomm=FM.worldsComm(),
                                          _extended_doftable=extendedDT,
                                          _periodicity=theperiodicity );
}

} // namespace detail
#endif

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    this->log("FluidMechanics","initFunctionSpaces","start");
    this->timerTool("Constructor").start();

    // maybe build extended dof table
    std::vector<bool> extendedDT( 2,false );
    bool hasExtendedDofTable = false;
    if ( (this->doCIPStabConvection() || this->doCIPStabDivergence()) && !this->applyCIPStabOnlyOnBoundaryFaces() )
    {
        this->log("FluidMechanics","createFunctionSpaces", "use buildDofTableMPIExtended on velocity" );
        extendedDT[0] = true;
        hasExtendedDofTable = true;
    }
    if ( this->doCIPStabPressure() )
    {
        this->log("FluidMechanics","createFunctionSpaces", "use buildDofTableMPIExtended on pressure" );
        extendedDT[1] = true;
        hasExtendedDofTable = true;
    }

    // fluid spaces : velocity and pressure
    auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
    if ( mom->isDefinedOnWholeMesh( this->physicsAvailableFromCurrentType() ) )
    {
        M_rangeMeshElements = elements(this->mesh());
        M_XhVelocity = space_velocity_type::New( _mesh=this->mesh(),
                                                 _extended_doftable=extendedDT[0] );
        M_XhPressure = space_pressure_type::New( _mesh=this->mesh(),
                                                 _extended_doftable=extendedDT[1] );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), mom->markers( this->physicsAvailableFromCurrentType() ));
        M_XhVelocity = space_velocity_type::New( _mesh=this->mesh(),
                                                 _extended_doftable=extendedDT[0],
                                                 _range=M_rangeMeshElements );
        M_XhPressure = space_pressure_type::New( _mesh=this->mesh(),
                                                 _extended_doftable=extendedDT[1],
                                                 _range=M_rangeMeshElements );
    }


    M_fieldVelocity = M_XhVelocity->elementPtr("velocity");
    M_fieldPressure = M_XhPressure->elementPtr("pressure");

    double tElapsed = this->timerTool("Constructor").stop("initFunctionSpaces");
    this->log("FluidMechanics","createFunctionSpaces", (boost::format("finish in %1% s") %tElapsed).str() );
}

//---------------------------------------------------------------------------------------------------------//


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::createALE()
{
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( !this->markerALEMeshBC( "moving" ).empty() )
    {
        this->log("FluidMechanics","createALE", "start" );
        this->timerTool("Constructor").start();

        M_isMoveDomain=true;

        M_meshALE = meshale( _mesh=this->mesh(),_prefix=this->prefix(),_directory=this->repository() );
        if ( this->functionSpaceVelocity()->dof()->hasMeshSupport() && this->functionSpaceVelocity()->dof()->meshSupport()->isPartialSupport() )
            M_meshALE->setComputationalDomain( this->keyword(), M_rangeMeshElements );
        else
            M_meshALE->setWholeMeshAsComputationalDomain( this->keyword() );
        M_meshALE->init();
        this->log("FluidMechanics","createALE", "create meshale object done" );
        // mesh displacement only on moving
        //M_meshDisplacementOnInterface.reset( new element_mesh_disp_type(M_meshALE->displacement()->functionSpace(),"mesh_disp_on_interface") );
        // mesh velocity used with stab CIP terms (need extended dof table)
        if ( this->doCIPStabConvection() )
            M_fieldMeshVelocityUsedWithStabCIP.reset( new element_velocity_type( this->functionSpaceVelocity() ) );

        double tElapsed = this->timerTool("Constructor").stop("createALE");
        this->log("FluidMechanics","createALE", (boost::format("finish in %1% s") %tElapsed).str() );
    }
#endif

}

//---------------------------------------------------------------------------------------------------------//
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    // clear
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    this->clearMarkerALEMeshBC();
    this->clearMarkerSlipBC();
    this->clearMarkerPressureBC();
    this->M_fluidOutletsBCType.clear();

    // boundary conditions
    this->M_isMoveDomain = false;
    this->M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "velocity", "Dirichlet" );
    for( auto const& d : this->M_bcDirichlet )
    {
        std::pair<bool,std::string> dirichletbcTypeRead = this->modelProperties().boundaryConditions().sparam( "velocity", "Dirichlet", name(d), "method" );
        std::string dirichletbcType = ( dirichletbcTypeRead.first )? dirichletbcTypeRead.second : soption(_name="dirichletbc.type",_prefix=this->prefix());
        CHECK( dirichletbcType=="elimination" || dirichletbcType=="nitsche" || dirichletbcType=="lm" ) << "invalid dirichletbc.type " << dirichletbcType;

        this->setMarkerDirichletBCByNameId( dirichletbcType, name(d), markers(d),ComponentType::NO_COMPONENT );

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Dirichlet", name(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        this->addMarkerALEMeshBC(bcTypeMeshALE,markers(d) );

        std::pair<bool,std::string> bcTypeTurbulenceRead = this->modelProperties().boundaryConditions().sparam( "velocity", "Dirichlet", name(d), "turbulence_bc" );
        if ( bcTypeTurbulenceRead.first )
        {
            if ( bcTypeTurbulenceRead.second == "inlet" )
            {
                typename TurbulenceModelBoundaryConditions::Inlet bcTurbInlet;
                bcTurbInlet.addMarkers( markers(d) );
                M_turbulenceModelBoundaryConditions.addInlet( name(d), bcTurbInlet );
            }
            else if ( bcTypeTurbulenceRead.second == "wall" )
            {
                typename TurbulenceModelBoundaryConditions::Wall bcTurbWall;
                bcTurbWall.addMarkers( markers(d) );
                M_turbulenceModelBoundaryConditions.addWall( name(d), bcTurbWall );
            }
        }
    }
    for ( ComponentType comp : std::vector<ComponentType>( { ComponentType::X, ComponentType::Y, ComponentType::Z } ) )
    {
        std::string compTag = ( comp ==ComponentType::X )? "x" : (comp == ComponentType::Y )? "y" : "z";
        std::string bcDirichletCompField = (boost::format("velocity_%1%")%compTag).str();
        std::string bcDirichletCompKeyword = "Dirichlet";
        this->M_bcDirichletComponents[comp] = this->modelProperties().boundaryConditions().getScalarFields( { { bcDirichletCompField, bcDirichletCompKeyword } } );
        for( auto const& d : this->M_bcDirichletComponents.find(comp)->second )
        {
            std::pair<bool,std::string> dirichletbcTypeRead = this->modelProperties().boundaryConditions().sparam( bcDirichletCompField, bcDirichletCompKeyword, name(d), "method" );
            std::string dirichletbcType = ( dirichletbcTypeRead.first )? dirichletbcTypeRead.second : soption(_name="dirichletbc.type",_prefix=this->prefix());
            CHECK( dirichletbcType=="elimination" || dirichletbcType=="nitsche" || dirichletbcType=="lm" ) << "invalid dirichletbc.type " << dirichletbcType;

            this->setMarkerDirichletBCByNameId( dirichletbcType, name(d), markers(d), comp );

            std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( bcDirichletCompField, bcDirichletCompKeyword, name(d), "alemesh_bc" );
            std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
            this->addMarkerALEMeshBC(bcTypeMeshALE,markers(d));

            std::pair<bool,std::string> bcTypeTurbulenceRead = this->modelProperties().boundaryConditions().sparam( bcDirichletCompField, "Dirichlet", name(d), "turbulence_bc" );
            if ( bcTypeTurbulenceRead.first )
            {
                if ( bcTypeTurbulenceRead.second == "inlet" )
                {
                    typename TurbulenceModelBoundaryConditions::Inlet bcTurbInlet;
                    bcTurbInlet.addMarkers( markers(d) );
                    M_turbulenceModelBoundaryConditions.addInlet( name(d), bcTurbInlet );
                }
                else if ( bcTypeTurbulenceRead.second == "wall" )
                {
                    typename TurbulenceModelBoundaryConditions::Wall bcTurbWall;
                    bcTurbWall.addMarkers( markers(d) );
                    M_turbulenceModelBoundaryConditions.addWall( name(d), bcTurbWall );
                }
            }

        }
    }

    this->M_bcNeumannScalar = this->modelProperties().boundaryConditions().getScalarFields( "velocity", "Neumann_scalar" );
    for( auto const& d : this->M_bcNeumannScalar )
    {
        this->setMarkerNeumannBC(NeumannBCShape::SCALAR,name(d),markers(d));

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Neumann_scalar", name(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        this->addMarkerALEMeshBC(bcTypeMeshALE,markers(d));
    }
    this->M_bcNeumannVectorial = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "velocity", "Neumann_vectorial" );
    for( auto const& d : this->M_bcNeumannVectorial )
    {
        this->setMarkerNeumannBC(NeumannBCShape::VECTORIAL,name(d),markers(d));

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Neumann_vectorial", name(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        this->addMarkerALEMeshBC(bcTypeMeshALE,markers(d));
    }
    this->M_bcNeumannTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<nDim>( "velocity", "Neumann_tensor2" );
    for( auto const& d : this->M_bcNeumannTensor2 )
    {
        this->setMarkerNeumannBC(NeumannBCShape::TENSOR2,name(d),markers(d));

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Neumann_tensor2", name(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        this->addMarkerALEMeshBC(bcTypeMeshALE,markers(d));
    }

    this->M_bcPressure = this->modelProperties().boundaryConditions().getScalarFields( "pressure", "Dirichlet" );
    for( auto const& d : this->M_bcPressure )
    {
        this->setMarkerPressureBC(name(d),markers(d));

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "pressure", "Dirichlet", name(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        this->addMarkerALEMeshBC(bcTypeMeshALE,markers(d));
    }
    for( std::string const& bcMarker : this->modelProperties().boundaryConditions().markers("fluid", "slip") )
    {
        this->addMarkerSlipBC( bcMarker );
        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "slip", bcMarker, "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        this->addMarkerALEMeshBC(bcTypeMeshALE,bcMarker);
    }
    for( std::string const& bcMarker : this->modelProperties().boundaryConditions().markers("fluid", "outlet") )
    {
        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "fluid", "outlet", bcMarker, "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");

        std::string typeOutlet = soption(_name="fluid-outlet.type", _prefix=this->prefix());//"free";
        std::pair<bool,std::string> typeOutletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "outlet", bcMarker, "model" );
        if ( typeOutletRead.first )
        {
            typeOutlet = typeOutletRead.second;
            CHECK( typeOutlet == "free" || typeOutlet == "windkessel" ) << "invalid outlet model " << typeOutlet;
        }
        std::string typeCouplingWindkesselOutlet = soption(_name="fluid-outlet.windkessel.coupling", _prefix=this->prefix());
        std::pair<bool,std::string> typeCouplingWindkesselOutletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "outlet", bcMarker, "windkessel_coupling" );
        if ( typeCouplingWindkesselOutletRead.first )
        {
            typeCouplingWindkesselOutlet = typeCouplingWindkesselOutletRead.second;
            CHECK( typeCouplingWindkesselOutlet == "implicit" || typeCouplingWindkesselOutlet == "explicit" ) << "invalid windkessel coupling type " << typeCouplingWindkesselOutlet;
        }
        std::pair<bool,double> WindkesselRdRead = this->modelProperties().boundaryConditions().dparam( "fluid", "outlet", bcMarker, "windkessel_Rd" );
        std::pair<bool,double> WindkesselRpRead = this->modelProperties().boundaryConditions().dparam( "fluid", "outlet", bcMarker, "windkessel_Rp" );
        std::pair<bool,double> WindkesselCdRead = this->modelProperties().boundaryConditions().dparam( "fluid", "outlet", bcMarker, "windkessel_Cd" );
        double WindkesselRd = ( WindkesselRdRead.first )? WindkesselRdRead.second : 1.;
        double WindkesselRp = ( WindkesselRpRead.first )? WindkesselRpRead.second : 1.;
        double WindkesselCd = ( WindkesselCdRead.first )? WindkesselCdRead.second : 1.;

        std::tuple<std::string,double,double,double> windkesselParam = std::make_tuple(typeCouplingWindkesselOutlet,WindkesselRd,WindkesselRp,WindkesselCd);

        this->M_fluidOutletsBCType.push_back(std::make_tuple(bcMarker,typeOutlet, windkesselParam ));
        this->addMarkerALEMeshBC(bcTypeMeshALE,bcMarker);
    }
    for( std::string const& bcMarker : this->modelProperties().boundaryConditions().markers("fluid", "inlet") )
    {
        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "fluid", "inlet", bcMarker, "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");

        std::string shapeInlet;
        std::pair<bool,std::string> shapeInletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "inlet", bcMarker, "shape" );
        if ( shapeInletRead.first )
        {
            shapeInlet = shapeInletRead.second;
            CHECK( shapeInlet == "constant" || shapeInlet == "parabolic" ) << "invalid inlet shape " << shapeInlet;
        }
        else
            CHECK( false ) << "inlet shape not given";

        std::string constraintInlet;
        std::pair<bool,std::string> constraintInletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "inlet", bcMarker, "constraint" );
        if ( constraintInletRead.first )
        {
            constraintInlet = constraintInletRead.second;
            CHECK( constraintInlet == "velocity_max" || constraintInlet == "flow_rate" ) << "invalid inlet constraint " << constraintInlet;
        }
        else
            CHECK( false ) << "inlet constraint not given";

        std::string fullTypeInlet = (boost::format("%1%_%2%")%constraintInlet %shapeInlet).str();

        std::string exprFluidInlet;
        std::pair<bool,std::string> exprFluidInletRead = this->modelProperties().boundaryConditions().sparam( "fluid", "inlet", bcMarker, "expr" );
        if ( exprFluidInletRead.first )
            exprFluidInlet = exprFluidInletRead.second;
        else
            CHECK( false ) << "inlet expr not given";

        this->M_fluidInletDesc.push_back(std::make_tuple(bcMarker,fullTypeInlet, expr<2>( exprFluidInlet,"",this->worldComm(),this->repository().expr() )) );
        this->addMarkerALEMeshBC(bcTypeMeshALE,bcMarker);

        typename TurbulenceModelBoundaryConditions::Inlet bcTurbInlet;
        bcTurbInlet.addMarkers( bcMarker/*markers(d)*/ );
        M_turbulenceModelBoundaryConditions.addInlet( bcMarker/*name(d)*/, bcTurbInlet );
    }

    M_bcMovingBoundaryImposed = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "fluid", "moving_boundary_imposed" );
    for( auto const& d : M_bcMovingBoundaryImposed )
    {
        for( std::string const& bcMarker : markers(d) )
            this->addMarkerALEMeshBC("moving",bcMarker);

        std::string dirichletbcType = "elimination";
        M_bcMarkersMovingBoundaryImposed.setMarkerDirichletBCByNameId( dirichletbcType, name(d), markers(d),ComponentType::NO_COMPONENT );
    }

    for( std::string const& bcMarker : this->modelProperties().boundaryConditions().markers( { { "velocity", "interface_fsi" }, { "fluid","interface_fsi"} } ) )
    {
        this->addMarkerALEMeshBC("moving",bcMarker);
        M_markersFSI.insert( bcMarker );
    }

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "fluid", "VolumicForces" );



    if ( auto _bcPTree = this->modelProperties().pTree().get_child_optional("BoundaryConditions") )
    {
        if ( auto _fluidPTree = _bcPTree->get_child_optional("fluid") )
        {
            if ( auto _bodyPtree = _fluidPTree->get_child_optional("body") )
            {
                for ( auto const& item : *_bodyPtree )
                {
                    std::string bodyName = item.first;
                    BodyBoundaryCondition bpbc( *this );
                    bpbc.setup( bodyName, item.second, *this );
                    if ( true ) // check if setup is enough
                    {
                        M_bodySetBC.emplace(bpbc.name(), bpbc );
                        for( std::string const& bcMarker : bpbc.markers() )
                            this->addMarkerALEMeshBC("moving",bcMarker);
                    }
                }
            }
        }
    }

    // Dirichlet bc using a lagrange multiplier
    if ( this->hasMarkerDirichletBClm() )
    {
        bool useSubMeshRelation = boption(_name="dirichletbc.lm.use-submesh-relation",_prefix=this->prefix());
        size_type useSubMeshRelationKey = (useSubMeshRelation)? EXTRACTION_KEEP_MESH_RELATION : 0;
        M_meshDirichletLM = createSubmesh(_mesh=this->mesh(),_range=markedfaces(this->mesh(),this->markerDirichletBClm()), _context=useSubMeshRelationKey,_view=true );
        if ( boption(_name="dirichletbc.lm.savemesh",_prefix=this->prefix()) )
        {
            std::string nameMeshDirichletLM = "nameMeshDirichletLM.msh";
            saveGMSHMesh(_mesh=M_meshDirichletLM,_filename=nameMeshDirichletLM);
        }
        M_XhDirichletLM = space_trace_velocity_type::New( _mesh=M_meshDirichletLM );
        M_fieldDirichletLM = M_XhDirichletLM->elementPtr();
    }

    // lagrange multiplier for pressure bc
    if ( this->hasMarkerPressureBC() )
    {
        M_meshLagrangeMultiplierPressureBC = createSubmesh(_mesh=this->mesh(),_range=markedfaces(this->mesh(),this->markerPressureBC()),_view=true );
        M_spaceLagrangeMultiplierPressureBC = space_trace_velocity_component_type::New( _mesh=M_meshLagrangeMultiplierPressureBC );
        M_fieldLagrangeMultiplierPressureBC1 = M_spaceLagrangeMultiplierPressureBC->elementPtr();
        if constexpr ( nDim == 3 )
            M_fieldLagrangeMultiplierPressureBC2 =  M_spaceLagrangeMultiplierPressureBC->elementPtr();
    }

    // init fluid outlet
    // this->initFluidOutlet();  (MOVE in init but should be fixed : TODO)
    // init fluid inlet
    this->initFluidInlet();

    // init bc body
    M_bodySetBC.init( *this );


    this->updateBoundaryConditionsForUse();
}
//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::createPostProcessExporters()
{
    this->log("FluidMechanics","createPostProcessExporters", "start" );

    //bool doExport = boption(_name="exporter.export");
    //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;//(this->isMoveDomain())?ExporterGeometry::EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY:ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
#if 0
    std::string geoExportType="static";//change_coords_only, change, static
#else
    bool useStaticExporter = boption(_name="exporter.use-static-mesh",_prefix=this->prefix());
    std::string geoExportType = useStaticExporter? "static":"change";
#endif

    if constexpr ( nOrderGeo <= 2 /*&& doExport*/ )
    {
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export",
                               //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export")),
                               _geo=geoExportType,
                               _path=this->exporterPath() );
    }

    if ( M_isHOVisu /*&& doExport*/ )
    {
#if 1 //defined(FEELPP_HAS_VTK)
        //M_exporter_ho = export_ho_type::New( this->application()->vm(), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Export_HO"))/*.c_str()*/, M_Xh->worldComm() );

// #if defined( FEELPP_MODELS_HAS_MESHALE )
//         if (M_isMoveDomain) this->meshALE()->revertReferenceMesh();
// #endif
        //auto Xh_create_ho = space_create_ho_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );

        std::shared_ptr<mesh_visu_ho_type> meshVisuHO;
        std::string hovisuSpaceUsed = soption(_name="hovisu.space-used",_prefix=this->prefix());
        bool doLagP1parallel=false;
        if ( hovisuSpaceUsed == "velocity" )
        {
            // with velocity field
            auto Xh_create_ho = M_XhVelocity->compSpace();
            auto opLagP1 = lagrangeP1( _space=Xh_create_ho,
                                       _backend=this->algebraicBackend(),
                                       //_worldscomm=this->localNonCompositeWorldsComm(),
                                       _path=this->rootRepository(),
                                       _prefix=this->prefix(),
                                       _rebuild=!this->doRestart(),
                                       _parallel=doLagP1parallel );
            meshVisuHO = opLagP1->mesh();
        }
        else if ( hovisuSpaceUsed == "pressure" )
        {
            // with pressure velocity field
            auto Xh_create_ho = M_XhPressure;
            auto opLagP1 = lagrangeP1( _space=Xh_create_ho,
                                       _backend=this->algebraicBackend(),
                                       //_worldscomm=this->localNonCompositeWorldsComm(),
                                       _path=this->rootRepository(),
                                       _prefix=this->prefix(),
                                       _rebuild=!this->doRestart(),
                                       _parallel=doLagP1parallel );
            meshVisuHO = opLagP1->mesh();
        }
        else if ( hovisuSpaceUsed == "p1" )
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            if ( M_meshALE )
                meshVisuHO = M_meshALE->referenceMesh();
            else
                meshVisuHO = this->mesh()->createP1mesh();
#else
            meshVisuHO = this->mesh()->createP1mesh();
#endif
        }
        else CHECK( false ) << "invalid hovisu.space-used " << hovisuSpaceUsed;

        M_exporter_ho = exporter( _mesh=meshVisuHO,//opLagP1->mesh(),
                                  //_name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"ExportHO")),
                                  _name="ExportHO",
                                  _geo=geoExportType,
                                  _path=this->exporterPath() );


        M_XhVectorialVisuHO = space_vectorial_visu_ho_type::New(_mesh=meshVisuHO/*opLagP1->mesh()*/, _worldscomm=this->localNonCompositeWorldsComm());
        //M_XhScalarVisuHO = space_scalar_visu_ho_type::New(_mesh=opLagP1->mesh(),_worldscomm=this->localNonCompositeWorldsComm());
        M_XhScalarVisuHO = M_XhVectorialVisuHO->compSpace();

        M_velocityVisuHO.reset( new element_vectorial_visu_ho_type(M_XhVectorialVisuHO,"u_visuHO"));
        M_pressureVisuHO.reset( new element_scalar_visu_ho_type(M_XhScalarVisuHO,"p_visuHO"));
        if (M_isMoveDomain) M_meshdispVisuHO.reset( new element_vectorial_visu_ho_type(M_XhVectorialVisuHO,"meshdisp_visuHO"));

        this->log("FluidMechanics","createPostProcessExporters", "start opInterpolation" );
        boost::mpi::timer timerOpI;

        M_opIvelocity = opInterpolation(_domainSpace=M_XhVelocity,
                                        _imageSpace=M_XhVectorialVisuHO,
                                        _range=elements(M_XhVectorialVisuHO->mesh()),
                                        _backend=this->algebraicBackend(),
                                        _type=InterpolationNonConforme(false,true,false,15) );

        this->log("FluidMechanics","createPostProcessExporters", "step1 done" );

        M_opIpressure = opInterpolation(_domainSpace=M_XhPressure,
                                        _imageSpace=M_XhScalarVisuHO,
                                        _range=elements(M_XhScalarVisuHO->mesh()),
                                        _backend=this->algebraicBackend(),
                                        _type=InterpolationNonConforme(false,true,false,15) );
#if 0
        if ( this->hasPostProcessFieldExported( "normal-stress" ) ||
             this->hasPostProcessFieldExported( "wall-shear-stress" ) )
        {
            M_XhVectorialDiscVisuHO = space_vectorialdisc_visu_ho_type::New(_mesh=meshVisuHO/*opLagP1->mesh()*/,_worldscomm=this->localNonCompositeWorldsComm());
            if ( this->hasPostProcessFieldExported( "normal-stress" ) )
                M_normalStressVisuHO.reset( new element_vectorialdisc_visu_ho_type(M_XhVectorialDiscVisuHO,"normalstress_visuHO") );
            if ( this->hasPostProcessFieldExported( "wall-shear-stress" ) )
                M_fieldWallShearStressVisuHO.reset( new element_vectorialdisc_visu_ho_type(M_XhVectorialDiscVisuHO,"wallshearstress_visuHO") );
            M_opIstress = opInterpolation(_domainSpace=M_XhNormalBoundaryStress,
                                          _imageSpace=M_XhVectorialDiscVisuHO,
                                          _range=elements(M_XhVectorialDiscVisuHO->mesh()),
                                          _backend=this->algebraicBackend(),
                                          _type=InterpolationNonConforme(false,true,false,15) );
        }
#endif
        this->log("FluidMechanics","createPostProcessExporters", "step2 done" );

        if (M_isMoveDomain )
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            M_opImeshdisp = opInterpolation(_domainSpace=M_meshALE->functionSpace(),
                                            _imageSpace=M_XhVectorialVisuHO,
                                            _range=elements(M_XhVectorialVisuHO->mesh()),
                                            _backend=this->algebraicBackend(),
                                            _type=InterpolationNonConforme(false,true,false,15) );
#endif
        }

        double timeElapsedOpI = timerOpI.elapsed();
        this->log("FluidMechanics","createPostProcessExporters", "finish all opInterpolation in " + (boost::format("%1% s") % timeElapsedOpI).str() );

// #if defined( FEELPP_MODELS_HAS_MESHALE )
//         if (M_isMoveDomain) this->meshALE()->revertMovingMesh();
// #endif

#endif
    }

    this->log("FluidMechanics","createPostProcessExporters", "finish" );

} // createPostProcessExporters

//---------------------------------------------------------------------------------------------------------//

// FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
// void
// FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::createFunctionSpacesNormalStress()
// {
//     if ( M_XhNormalBoundaryStress ) return;
//     M_XhNormalBoundaryStress = space_normalstress_type::New( _mesh=M_meshTrace );
//     M_fieldNormalStress.reset(new element_normalstress_type(M_XhNormalBoundaryStress));
//     M_fieldWallShearStress.reset(new element_normalstress_type(M_XhNormalBoundaryStress));
// }

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::createFunctionSpacesSourceAdded()
{
    if ( this->functionSpaceVelocity()->dof()->hasMeshSupport() && this->functionSpaceVelocity()->dof()->meshSupport()->isPartialSupport() ) //M_materialProperties->isDefinedOnWholeMesh() )
        M_XhSourceAdded=space_vectorial_PN_type::New( _mesh=this->mesh(),_worldscomm=this->localNonCompositeWorldsComm(),
                                                       _range=M_rangeMeshElements );
    else
        M_XhSourceAdded=space_vectorial_PN_type::New( _mesh=this->mesh(),_worldscomm=this->localNonCompositeWorldsComm() );
    M_SourceAdded.reset( new element_vectorial_PN_type(M_XhSourceAdded,"SourceAdded"));
}

//---------------------------------------------------------------------------------------------------------//


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initFluidInlet()
{
    if ( !this->hasFluidInlet() ) return;

    M_fluidInletMesh.clear();
    M_fluidInletSpace.clear();
    M_fluidInletVelocity.clear();
    for ( auto const& inletbc : M_fluidInletDesc )
    {
        std::string const& marker = std::get<0>( inletbc );
        std::string const& type = std::get<1>( inletbc );
        auto const& valMaxExpr = std::get<2>( inletbc );
        auto meshinlet = createSubmesh( _mesh=this->mesh(),_range=markedfaces(this->mesh(),marker), _view=true );
        auto spaceinlet = space_fluidinlet_type::New( _mesh=meshinlet,_worldscomm=this->localNonCompositeWorldsComm() );
        auto velinlet = spaceinlet->elementPtr();
        auto velinletInterpolated = functionSpaceVelocity()->compSpace()->elementPtr();
        auto opIfluidinlet = opInterpolation(_domainSpace=spaceinlet,
                                             _imageSpace=this->functionSpaceVelocity()->compSpace(),
                                             _range=markedfaces(this->mesh(),marker),
                                             _backend=this->backend() );
        M_fluidInletMesh[marker] = meshinlet;
        M_fluidInletSpace[marker] = spaceinlet;
        M_fluidInletVelocity[marker] = velinlet;
        M_fluidInletVelocityInterpolated[marker] = std::make_tuple(velinletInterpolated,opIfluidinlet);

        double areainlet = integrate(_range=elements(meshinlet),
                                     _expr=cst(1.)).evaluate()(0,0);
        auto velinletRef = spaceinlet->elementPtr();
        double maxVelRef = 0.;
        if ( type == "velocity_max_constant" || type == "flow_rate_constant" )
        {
            maxVelRef = areainlet;
            velinletRef->on(_range=elements(meshinlet),_expr=cst(areainlet) );
            velinletRef->on(_range=boundaryfaces(meshinlet),_expr=cst(0.) );
        }
        else if ( type == "velocity_max_parabolic" || type == "flow_rate_parabolic" )
        {
            auto l = form1( _test=spaceinlet );
            l = integrate(_range=elements(meshinlet),
                          _expr=cst(areainlet)*id(velinlet));
            auto a = form2( _trial=spaceinlet, _test=spaceinlet);
            a = integrate(_range=elements(meshinlet),
                          _expr=gradt(velinlet)*trans(grad(velinlet)) );
            a+=on(_range=boundaryfaces(meshinlet), _rhs=l, _element=*velinlet, _expr=cst(0.) );

            auto backendinlet = backend_type::build( soption( _name="backend" ), prefixvm(this->prefix(),"fluidinlet"), this->worldCommPtr() );
            backendinlet->solve(_matrix=a.matrixPtr(),_rhs=l.vectorPtr(),_solution=*velinletRef );
            maxVelRef = velinletRef->max();
        }
        double flowRateRef = integrate(_range=markedfaces(this->mesh(),marker),
                                       _expr=inner( idv(velinletRef)*N(),N() ) ).evaluate()(0,0);
        M_fluidInletVelocityRef[marker] = std::make_tuple(velinletRef,maxVelRef,flowRateRef);

    }

    this->updateFluidInletVelocity();
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateFluidInletVelocity()
{
   for ( auto & inletbc : M_fluidInletDesc )
   {
        std::string const& marker = std::get<0>( inletbc );
        std::string const& type = std::get<1>( inletbc );
        auto & exprFluidInlet = std::get<2>( inletbc );
        exprFluidInlet.setParameterValues( this->modelProperties().parameters().toParameterValues() );
        double evalExprFluidInlet = expr( exprFluidInlet, this->symbolsExpr() ).evaluate()(0,0);

        auto itMesh = M_fluidInletMesh.find( marker );
        CHECK( itMesh != M_fluidInletMesh.end() ) << "fluid inlet not init for this marker" << marker;
        auto meshinlet = itMesh->second;

        auto itVelRef = M_fluidInletVelocityRef.find(marker);
        CHECK( itVelRef != M_fluidInletVelocityRef.end() ) << "fluid inlet not init for this marker" << marker;
        auto const& velRef = std::get<0>(itVelRef->second);
        double maxVelRef = std::get<1>(itVelRef->second);
        double flowRateRef = std::get<2>(itVelRef->second);

        if ( type == "velocity_max_constant" || type == "velocity_max_parabolic" )
        {
            M_fluidInletVelocity[marker]->zero();
            M_fluidInletVelocity[marker]->add( evalExprFluidInlet/maxVelRef, *velRef );
            //M_fluidInletVelocity[marker]->on(_range=elements(meshinlet),_expr=cst(evalExprFluidInlet) );
            //M_fluidInletVelocity[marker]->on(_range=boundaryfaces(meshinlet),_expr=cst(0.) );

        }
        else if ( type == "flow_rate_constant" || type == "flow_rate_parabolic" )
        {
            M_fluidInletVelocity[marker]->zero();
            M_fluidInletVelocity[marker]->add( evalExprFluidInlet/flowRateRef, *velRef );
        }

        auto const& velSubmesh = M_fluidInletVelocity.find(marker)->second;
        auto opI = std::get<1>( M_fluidInletVelocityInterpolated[marker] );
        auto & velInterp = std::get<0>( M_fluidInletVelocityInterpolated[marker] );
        opI->apply( *velSubmesh , *velInterp );

#if 0
        double flowRateComputed = integrate(_range=markedfaces(this->mesh(),marker),
                                            _expr=-idv(velInterp)*N() ).evaluate()(0,0);
        double maxVelComputed = velInterp->max();
        if ( this->worldComm().isMasterRank() )
            std::cout << "flowRateComputed : " << flowRateComputed << "\n"
                      << "maxVelComputed : " << maxVelComputed << "\n";
#endif
   }
}


namespace detail
{
template <typename SpaceType>
NullSpace<double> getNullSpace( SpaceType const& space, mpl::int_<2> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( vec(Py(),-Px()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3 } );
    return userNullSpace;
}
template <typename SpaceType>
NullSpace<double> getNullSpace( SpaceType const& space, mpl::int_<3> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( oneZ() );
    auto mode4 = space->element( vec(Py(),-Px(),cst(0.)) );
    auto mode5 = space->element( vec(-Pz(),cst(0.),Px()) );
    auto mode6 = space->element( vec(cst(0.),Pz(),-Py()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3,mode4,mode5,mode6 } );
    return userNullSpace;
}

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    if ( this->isUpdatedForUse() ) return;

    this->log("FluidMechanics","init", "start" );
    this->timerTool("Constructor").start();

    if ( this->physics().empty() )
        this->initPhysics( this->keyword(), this->modelProperties().models() );

    this->initMaterialProperties();

    if ( !this->mesh() )
        this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    // backend
    this->initAlgebraicBackend();

    if ( M_solverName == "automatic" )
    {
        bool isLinear = true;
        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
            if ( physicFluidData->equation() == "Navier-Stokes" || !physicFluidData->dynamicViscosity().isNewtonianLaw() )
            {
                isLinear = false;
                break;
            }
        }
        if ( isLinear || M_useSemiImplicitTimeScheme )
            M_solverName="LinearSystem";
        else
            M_solverName="Newton";
    }

    // functionSpaces and elements
    this->initFunctionSpaces();

    this->initBoundaryConditions();

    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();

    // init fluid outlet
    this->initFluidOutlet(); // (MOVE in initBoundaryConditions but should be fixed : TODO)  because defined time steping inside

    // ALE mode (maybe)
    this->createALE();

    // update definePressureCst respect to the method choosen
    if ( this->definePressureCst() )
        this->updateDefinePressureCst();

    // update marker in mesh (mainly used with CIP stab)
    if ( (this->doCIPStabConvection() || this->doCIPStabDivergence() || this->doCIPStabPressure() ) && !this->applyCIPStabOnlyOnBoundaryFaces() )
        this->updateMarkedZonesInMesh();

    //-------------------------------------------------//
    // init stabilization
    if ( M_stabilizationGLS )
    {
        //static const uint16_type nStabGlsOrderPoly = (nOrderVelocity>1)? nOrderVelocity : 2;
        typedef StabilizationGLSParameter<mesh_type, nOrderVelocity> stab_gls_parameter_velocity_impl_type;
        typedef StabilizationGLSParameter<mesh_type, nOrderPressure> stab_gls_parameter_pressure_impl_type;
        M_stabilizationGLSParameterConvectionDiffusion.reset( new stab_gls_parameter_velocity_impl_type( this->mesh(),prefixvm(this->prefix(),"stabilization-gls.parameter") ) );
        M_stabilizationGLSParameterConvectionDiffusion->init();
        if ( nOrderVelocity == nOrderPressure )
             M_stabilizationGLSParameterPressure = M_stabilizationGLSParameterConvectionDiffusion;
        else
        {
            M_stabilizationGLSParameterPressure.reset( new stab_gls_parameter_pressure_impl_type( this->mesh(),prefixvm(this->prefix(),"stabilization-gls.parameter") ) );
            M_stabilizationGLSParameterPressure->init();
        }
        this->updateStabilizationGLSRange();
    }

    //-------------------------------------------------//
    // distance to wall
    if ( this->hasTurbulenceModel() )
        M_dist2WallEnabled = true;
    if ( M_dist2WallEnabled )
        this->initDist2Wall();

    if ( this->hasTurbulenceModel() )
        this->initTurbulenceModel();
    //-------------------------------------------------//
    // init function defined in json
    this->initUserFunctions();
    // init post-processinig (exporter, measure at point, ...)
    this->initPostProcess();
    //-------------------------------------------------//
    // update ALE mesh for use
    if (this->isMoveDomain())
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        for ( auto const& [bcName,markers] : this->markerALEMeshBC() )
            M_meshALE->addMarkersInBoundaryCondition( bcName, markers );


        std::set<std::string> markersbcMovingBoundaryImposed;
        for( auto const& d : M_bcMovingBoundaryImposed )
        {
            auto bcmarkers = markers(d);
            markersbcMovingBoundaryImposed.insert( bcmarkers.begin(), bcmarkers.end() );
        }
        if ( !markersbcMovingBoundaryImposed.empty() )
            M_meshALE->setDisplacementImposedOnInitialDomainOverFaces( this->keyword(), markersbcMovingBoundaryImposed );


        // if ( !M_bodySetBC.empty() )
        std::set<std::string> markersbcBodyOnlyBoundaryFaces, markersbcBodyWithElements;
        for ( auto const& [bname,bbc] : M_bodySetBC )
        {
            if ( bbc.body().hasMaterialsProperties() )
            {
                auto mom = bbc.body().materialsProperties()->materialsOnMesh( this->mesh() );
                auto markersOnPhysics = mom->markers( bbc.body().modelPhysics()->physicsAvailableFromCurrentType() );
                markersbcBodyWithElements.insert( markersOnPhysics.begin(), markersOnPhysics.end() );
            }
            else
            {
                markersbcBodyOnlyBoundaryFaces.insert( bbc.markers().begin(),  bbc.markers().end() );
            }
        }
        if ( !markersbcBodyOnlyBoundaryFaces.empty() )
            M_meshALE->setDisplacementImposedOnInitialDomainOverFaces( this->keyword(), markersbcBodyOnlyBoundaryFaces );

        if ( !markersbcBodyWithElements.empty() )
            M_meshALE->setDisplacementImposedOnInitialDomainOverElements( this->keyword(), markersbcBodyWithElements );


        if ( this->doRestart() )
            M_meshALE->revertMovingMesh();
        this->log("FluidMechanics","finish init", "meshALE done" );
#endif
    }

    //-------------------------------------------------//
    // bc body (call after meshALE->init() in case of restart)
    M_bodySetBC.updateForUse( *this );

    if ( M_useSemiImplicitTimeScheme )
        M_useVelocityExtrapolated = true;

    // update constant parameters
    this->updateParameterValues();

    // update initial conditions
    this->updateInitialConditions( this->symbolsExpr() );

    //-------------------------------------------------//
    // define start dof index ( lm , windkessel )
    this->initStartBlockIndexFieldsInMatrix();
    //-------------------------------------------------//
    // build solution block vector
    this->buildBlockVector();

    //-------------------------------------------------//
    // InHousePreconditioner : operatorPCD
    this->initInHousePreconditioner();

    //-------------------------------------------------//
    // algebraric data : solver, preconditioner, matrix, vector
    if ( buildModelAlgebraicFactory )
    {
        this->initAlgebraicFactory();
    }

    //-------------------------------------------------//
    //-------------------------------------------------//
    this->setIsUpdatedForUse( true );

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("FluidMechanics","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::applyRemesh( mesh_ptrtype const& newMesh )
{
    RemeshInterpolation remeshInterp;

    this->log("FluidMechanics","applyRemesh", "start" );
    mesh_ptrtype oldMesh = this->mesh();

    // material prop
    this->materialsProperties()->removeMesh( oldMesh );
    this->materialsProperties()->addMesh( newMesh );

    this->setMesh( newMesh );
    this->log("FluidMechanics","applyRemesh", "mesh done" );

    // function space and fields
    space_velocity_ptrtype old_spaceVelocity = M_XhVelocity;
    space_pressure_ptrtype old_spacePressure = M_XhPressure;
    element_velocity_ptrtype old_fieldVelocity = M_fieldVelocity;
    element_pressure_ptrtype old_fieldPressure = M_fieldPressure;
    this->initFunctionSpaces();
    this->log("FluidMechanics","applyRemesh", "functionspace done" );

    // createInterpolationOp
    auto matrixInterpolation_velocity = remeshInterp.computeMatrixInterpolation( old_spaceVelocity, M_XhVelocity, M_rangeMeshElements );
    remeshInterp.registeringBlockIndex( this->keyword(), this->startSubBlockSpaceIndex("velocity"), old_spaceVelocity, M_XhVelocity );
    remeshInterp.interpolate( old_fieldVelocity, M_fieldVelocity );

    auto matrixInterpolation_pressure = remeshInterp.computeMatrixInterpolation( old_spacePressure, M_XhPressure, M_rangeMeshElements );
    remeshInterp.registeringBlockIndex( this->keyword(), this->startSubBlockSpaceIndex("pressure"), old_spacePressure, M_XhPressure );
    remeshInterp.interpolate( old_fieldPressure, M_fieldPressure );

    this->log("FluidMechanics","applyRemesh", "interpolation velocity/pressure done" );

    if ( M_meshALE )
    {
        std::vector<std::tuple<std::string,elements_reference_wrapper_t<mesh_type>>> computationDomains;
        if ( this->functionSpaceVelocity()->dof()->hasMeshSupport() && this->functionSpaceVelocity()->dof()->meshSupport()->isPartialSupport() )
            computationDomains.push_back( std::make_tuple( this->keyword(), M_rangeMeshElements ) );
        M_meshALE->applyRemesh( newMesh, computationDomains );
        this->log("FluidMechanics","applyRemesh", "mesh ale done" );
    }

    // time stepping
    if ( M_bdfVelocity )
        M_bdfVelocity->applyRemesh( this->functionSpaceVelocity(), matrixInterpolation_velocity );
    if ( M_savetsPressure )
        M_savetsPressure->applyRemesh( this->functionSpacePressure(), matrixInterpolation_pressure );
    this->log("FluidMechanics","applyRemesh", "time stepping done" );

    // update definePressureCst respect to the method choosen
    if ( this->definePressureCst() )
        this->updateDefinePressureCst();

    // body bc
    M_bodySetBC.applyRemesh( *this, remeshInterp );
    this->log("FluidMechanics","applyRemesh", "body bc done" );

    // velocity imposed by using a lagrange multiplier
    if ( M_meshDirichletLM )
    {
        M_meshDirichletLM = createSubmesh(_mesh=this->mesh(),_range=markedfaces(this->mesh(),this->markerDirichletBClm())/*, _context=useSubMeshRelationKey*/,_view=true );
        space_trace_velocity_ptrtype old_XhDirichletLM = M_XhDirichletLM;
        element_trace_velocity_ptrtype old_fieldDirichletLM = M_fieldDirichletLM;
        M_XhDirichletLM = space_trace_velocity_type::New( _mesh=M_meshDirichletLM );
        M_fieldDirichletLM = M_XhDirichletLM->elementPtr();

        auto matrixInterpolation_lmvelocitybc = remeshInterp.computeMatrixInterpolation( old_XhDirichletLM, M_XhDirichletLM );
        remeshInterp.registeringBlockIndex( this->keyword(), this->startSubBlockSpaceIndex("dirichletlm"), old_XhDirichletLM, M_XhDirichletLM );
        remeshInterp.interpolate( old_fieldDirichletLM, M_fieldDirichletLM );
    }

    // pressure bc
    if ( M_meshLagrangeMultiplierPressureBC )
    {
        M_meshLagrangeMultiplierPressureBC = createSubmesh(_mesh=this->mesh(),_range=markedfaces(this->mesh(),this->markerPressureBC()),_view=true );
        space_trace_velocity_component_ptrtype old_spaceLagrangeMultiplierPressureBC = M_spaceLagrangeMultiplierPressureBC;
        element_trace_velocity_component_ptrtype old_fieldLagrangeMultiplierPressureBC1 = M_fieldLagrangeMultiplierPressureBC1;
        element_trace_velocity_component_ptrtype old_fieldLagrangeMultiplierPressureBC2 = M_fieldLagrangeMultiplierPressureBC2;
        M_spaceLagrangeMultiplierPressureBC = space_trace_velocity_component_type::New( _mesh=M_meshLagrangeMultiplierPressureBC );
        M_fieldLagrangeMultiplierPressureBC1 = M_spaceLagrangeMultiplierPressureBC->elementPtr();
        if constexpr ( nDim == 3 )
            M_fieldLagrangeMultiplierPressureBC2 = M_spaceLagrangeMultiplierPressureBC->elementPtr();

        auto matrixInterpolation_lmpressurebc = remeshInterp.computeMatrixInterpolation( old_spaceLagrangeMultiplierPressureBC, M_spaceLagrangeMultiplierPressureBC );
        remeshInterp.registeringBlockIndex( this->keyword(), this->startSubBlockSpaceIndex("pressurelm1"), old_spaceLagrangeMultiplierPressureBC, M_spaceLagrangeMultiplierPressureBC );
        remeshInterp.interpolate( old_fieldLagrangeMultiplierPressureBC1, M_fieldLagrangeMultiplierPressureBC1 );
        if constexpr ( nDim == 3 )
        {
            remeshInterp.registeringBlockIndex( this->keyword(), this->startSubBlockSpaceIndex("pressurelm2"), old_spaceLagrangeMultiplierPressureBC, M_spaceLagrangeMultiplierPressureBC );
            remeshInterp.interpolate( old_fieldLagrangeMultiplierPressureBC2, M_fieldLagrangeMultiplierPressureBC2 );
        }
    }

    // fluid inlet bc
    this->initFluidInlet();

    // velocity extrapolated
    if ( M_vectorVelocityExtrapolated )
    {
        vector_ptrtype old_vectorVelocityExtrapolated = M_vectorVelocityExtrapolated;
        M_vectorVelocityExtrapolated = this->algebraicBackend()->newVector( this->functionSpaceVelocity() );
        M_fieldVelocityExtrapolated = this->functionSpaceVelocity()->elementPtr( *M_vectorVelocityExtrapolated, 0 );
        matrixInterpolation_velocity->multVector( *old_vectorVelocityExtrapolated, *M_vectorVelocityExtrapolated );
    }
    if ( M_vectorPreviousVelocityExtrapolated )
    {
        vector_ptrtype old_vectorPreviousVelocityExtrapolated = M_vectorPreviousVelocityExtrapolated;
        M_vectorPreviousVelocityExtrapolated = this->algebraicBackend()->newVector( this->functionSpaceVelocity() );
        matrixInterpolation_velocity->multVector( *old_vectorPreviousVelocityExtrapolated, *M_vectorPreviousVelocityExtrapolated );
    }

    // stabilization GLS familily
    if ( M_stabilizationGLSParameterConvectionDiffusion )
    {
        M_stabilizationGLSParameterConvectionDiffusion->applyRemesh( this->mesh() );
        if ( M_stabilizationGLSParameterPressure != M_stabilizationGLSParameterConvectionDiffusion )
            M_stabilizationGLSParameterPressure->applyRemesh( this->mesh() );
        this->updateStabilizationGLSRange();
    }


    // TODO : post process ??

    // algebraic data/tools
    vector_ptrtype old_vectorPreviousSolution = M_vectorPreviousSolution;
    vector_ptrtype old_timeStepThetaSchemePreviousContrib = M_timeStepThetaSchemePreviousContrib;
    this->removeAllAlgebraicDataAndTools();
    this->initAlgebraicModel();
    this->initAlgebraicFactory();

    if ( old_timeStepThetaSchemePreviousContrib )
    {
        remeshInterp.interpolateBlockVector( this->keyword(), old_timeStepThetaSchemePreviousContrib, M_timeStepThetaSchemePreviousContrib, *this->algebraicBlockVectorSolution() );
        this->log("FluidMechanics","applyRemesh", "timeStepThetaSchemePreviousContrib done" );
    }

    if ( old_vectorPreviousSolution )
    {
        remeshInterp.interpolateBlockVector( this->keyword(), old_vectorPreviousSolution, M_vectorPreviousSolution, *this->algebraicBlockVectorSolution() );
        this->log("FluidMechanics","applyRemesh", "vectorPreviousSolution done" );
    }


    this->log("FluidMechanics","applyRemesh", "finish" );

}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initAlgebraicModel()
{
    // backend
    this->initAlgebraicBackend();

    // define start dof index ( lm , windkessel )
    this->initStartBlockIndexFieldsInMatrix();

    // build solution block vector
    this->buildBlockVector();

    // dof eliminitation ids
    this->updateAlgebraicDofEliminationIds();

    // InHousePreconditioner : operatorPCD
    this->initInHousePreconditioner();
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateAlgebraicDofEliminationIds()
{
    this->updateBoundaryConditionsForUse();
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
    this->setAlgebraicFactory( algebraicFactory );

    if ( !M_bodySetBC.empty() )
    {
        M_bodySetBC.initAlgebraicFactory( *this, this->algebraicFactory() );
    }

    if ( boption(_name="use-velocity-near-null-space",_prefix=this->prefix() ) )
    {
        std::string nearNullSpacePrefix = this->prefix();
        if ( Environment::vm().count(prefixvm(this->prefix(),"use-velocity-near-null-space.prefix").c_str()) )
            nearNullSpacePrefix = soption( _name="use-velocity-near-null-space.prefix", _prefix=this->prefix() );

        NullSpace<double> userNullSpace = detail::getNullSpace(this->functionSpaceVelocity(), mpl::int_<nDim>() ) ;
        algebraicFactory->attachNearNullSpace( 0,userNullSpace, nearNullSpacePrefix ); // for block velocity in fieldsplit
    }

    bool attachMassMatrix = boption(_prefix=this->prefix(),_name="preconditioner.attach-mass-matrix");
    if ( attachMassMatrix )
    {
        auto massbf = form2( _trial=this->functionSpaceVelocity(), _test=this->functionSpaceVelocity());
        //auto themassMatrixGraph = stencil( _trial=this->functionSpaceVelocity(), _test=this->functionSpaceVelocity() );
        auto const& u = this->fieldVelocity();
        massbf += integrate( _range=M_rangeMeshElements, _expr=inner( idt(u),id(u) ) );

        massbf.matrixPtr()->close();
        if ( this->algebraicFactory() )
            this->algebraicFactory()->attachAuxiliarySparseMatrix( "mass-matrix", massbf.matrixPtr() );
    }

    if ( this->hasOperatorPCD() )
        this->algebraicFactory()->attachOperatorPCD("pcd", this->operatorPCD());

    if ( this->timeStepping() == "Theta" )
    {
        M_timeStepThetaSchemePreviousContrib = this->backend()->newVector(this->algebraicBlockVectorSolution()->vectorMonolithic()->mapPtr() );
        algebraicFactory->addVectorResidualAssembly( M_timeStepThetaSchemePreviousContrib, 1.0, "Theta-Time-Stepping-Previous-Contrib", true );
        algebraicFactory->addVectorLinearRhsAssembly( M_timeStepThetaSchemePreviousContrib, -1.0, "Theta-Time-Stepping-Previous-Contrib", false );
        if ( M_stabilizationGLS )
            algebraicFactory->dataInfos().addVectorInfo( "time-stepping.previous-solution", M_vectorPreviousSolution/*this->backend()->newVector( M_blockVectorSolution.vectorMonolithic()->mapPtr() )*/ );
    }


}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateMarkedZonesInMesh()
{
    this->log("FluidMechanics","updateMarkedZonesInMesh", "start" );

    MarkedMeshTool<mesh_type> markMesh( this->mesh() );

    if ( Environment::vm().count( prefixvm(this->prefix(),"marked-zones.markedfaces" ) ) )
    {
        std::vector<std::string> mymarkedfaces = Environment::vm()[prefixvm(this->prefix(),"marked-zones.markedfaces").c_str()].template as<std::vector<std::string> >();
        markMesh.setFaceMarker( mymarkedfaces );
        markMesh.updateFaceMarker3FromFaceMarker();
        this->applyCIPStabOnlyOnBoundaryFaces( true );
    }
    if ( Environment::vm().count( prefixvm(this->prefix(),"marked-zones.elements-from-markedfaces" ) ) )
    {
        std::vector<std::string> mymarkedfaces = Environment::vm()[prefixvm(this->prefix(),"marked-zones.elements-from-markedfaces").c_str()].template as<std::vector<std::string> >();
        markMesh.setFaceMarker( mymarkedfaces );
        markMesh.updateFaceMarker3FromEltConnectedToFaceMarker();
        this->applyCIPStabOnlyOnBoundaryFaces( false );
    }
    if ( Environment::vm().count( prefixvm(this->prefix(),"marked-zones.expressions" ) ) )
    {
        std::vector<std::string> myexpressions = Environment::vm()[prefixvm(this->prefix(),"marked-zones.expressions").c_str()].template as<std::vector<std::string> >();
        for ( std::string const& mystringexpr : myexpressions )
        {
            auto myexpr = expr( mystringexpr );
            markMesh.updateFaceMarker3FromExpr(myexpr,false);
            this->applyCIPStabOnlyOnBoundaryFaces( false );
        }
        if ( myexpressions.size() >0 )
            markMesh.updateForUseFaceMarker3();
    }

    if ( boption(_name="marked-zones.internal-faces",_prefix=this->prefix() ) )
    {
        markMesh.updateFaceMarker3FromInternalFaces();
        this->applyCIPStabOnlyOnBoundaryFaces( false );
    }

    if ( this->verbose() )
    {
        markMesh.verbose();
    }

    if ( false )
    {
        markMesh.saveSubMeshFromMarked3Faces();
        markMesh.exportP0EltMarkerFromFaceMarker();
    }

    this->log("FluidMechanics","updateMarkedZonesInMesh", "finish" );
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateStabilizationGLSRange()
{
    if ( Environment::vm().count( prefixvm(this->prefix(),"stabilization-gls.convection-diffusion.location.expressions" ) ) )
    {
        std::string locationExpression = soption(_prefix=this->prefix(),_name="stabilization-gls.convection-diffusion.location.expressions");
        auto rangeStab = elements(this->mesh(),expr(locationExpression));
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physicsAvailableFromCurrentType() ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            M_stabilizationGLSEltRangeConvectionDiffusion[matName] = intersect( rangeStab, range );
        }
    }
    else
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physicsAvailableFromCurrentType() ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            M_stabilizationGLSEltRangeConvectionDiffusion[matName] = range;
        }
    }
    for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physicsAvailableFromCurrentType() ) )
    {
        auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
        M_stabilizationGLSEltRangePressure[matName] = range;
    }

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    this->log("FluidMechanics","initTimeStep", "start" );
    this->timerTool("Constructor").start();

    if ( this->isStationaryModel() ) // force BDF with Stokes
        M_timeStepping = "BDF";

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
#if 0
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
         suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    fs::path saveTsDir = fs::path(this->rootRepository())/fs::path( prefixvm(this->prefix(),prefixvm(this->subPrefix(),"ts")) );

    double ti = this->timeInitial();
    double tf = this->timeFinal();
    double dt = this->timeStep();
#endif
    int bdfOrder = 1;
    if ( M_timeStepping == "BDF" )
        bdfOrder = ioption(_prefix=this->prefix(),_name="bdf.order");
    int nConsecutiveSave = std::max( 3, bdfOrder ); // at least 3 is required when restart with theta scheme
#if 0
    M_bdfVelocity = bdf( _space=this->functionSpaceVelocity(),
                         _name="velocity"+suffixName,
                         _prefix=this->prefix(),
                         _order=bdfOrder,
                         // don't use the fluid.bdf {initial,final,step}time but the general bdf info, the order will be from fluid.bdf
                         _initial_time=ti, _final_time=tf, _time_step=dt,
                         _restart=this->doRestart(),
                         _restart_path=this->restartPath(),
                         _restart_at_last_save=this->restartAtLastSave(),
                         _save=this->tsSaveInFile(), _format=myFileFormat, _freq=this->tsSaveFreq(),
                         _n_consecutive_save=nConsecutiveSave );
    M_bdfVelocity->setfileFormat( myFileFormat );
    M_bdfVelocity->setPathSave( ( saveTsDir/"velocity" ).string() );

    M_savetsPressure = bdf( _space=this->functionSpacePressure(),
                            _name="pressure"+suffixName,
                            _prefix=this->prefix(),
                            _order=1,
                            _initial_time=ti, _final_time=tf, _time_step=dt,
                            _restart=this->doRestart(),
                            _restart_path=this->restartPath(),
                            _restart_at_last_save=this->restartAtLastSave(),
                            _save=this->tsSaveInFile(), _format=myFileFormat, _freq=this->tsSaveFreq(),
                            _n_consecutive_save=nConsecutiveSave );
    M_savetsPressure->setfileFormat( myFileFormat );
    M_savetsPressure->setPathSave( ( saveTsDir/"pressure" ).string() );
#else
    M_bdfVelocity = this->createBdf( this->functionSpaceVelocity(),"velocity", bdfOrder, nConsecutiveSave, myFileFormat );
    M_savetsPressure = this->createBdf( this->functionSpacePressure(),"pressure", 1, nConsecutiveSave, myFileFormat );
#endif

    // velocity extrapolated
    M_vectorVelocityExtrapolated = this->algebraicBackend()->newVector( this->functionSpaceVelocity() );
    M_vectorPreviousVelocityExtrapolated = this->algebraicBackend()->newVector( this->functionSpaceVelocity() );
    M_fieldVelocityExtrapolated = this->functionSpaceVelocity()->elementPtr( *M_vectorVelocityExtrapolated, 0 );


    double tir = M_bdfVelocity->timeInitial();
    if ( this->doRestart() )
    {
        // start time step
        tir = M_bdfVelocity->restart();
        M_savetsPressure->restart();
        // load a previous solution as current solution
        *M_fieldVelocity = M_bdfVelocity->unknown(0);
        *M_fieldPressure = M_savetsPressure->unknown(0);
    }

    M_bodySetBC.initTimeStep( *this, bdfOrder, nConsecutiveSave, myFileFormat );

    if ( this->doRestart() )
        this->setTimeInitial( tir );

    this->updateTime( tir );

#if 0
    // start or restart time step scheme
    if ( !this->doRestart() )
    {
        // up current time
        this->updateTime( M_bdfVelocity->timeInitial() );
    }
    else
    {
        // start time step
        double tir = M_bdfVelocity->restart();
        M_savetsPressure->restart();
        // load a previous solution as current solution
        *M_fieldVelocity = M_bdfVelocity->unknown(0);
        *M_fieldPressure = M_savetsPressure->unknown(0);
        // up initial time
        this->setTimeInitial( tir );
        // up current time
        this->updateTime( tir );

        this->log("FluidMechanics","initTimeStep", "restart bdf/exporter done" );
    }
#endif


    double tElapsed = this->timerTool("Constructor").stop("initTimeStep");
    this->log("FluidMechanics","initTimeStep", (boost::format("finish in %1% s") %tElapsed).str() );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initFluidOutlet()
{
    this->log("FluidMechanics","initFluidOutlet", "start" );

    // create submesh, functionspace and interpolation operator
    if ( this->hasFluidOutletWindkesselImplicit() )
    {
        // list usefull to create the outlets submesh
        std::list<std::string> markerNameBFOutletForSubmesh;
        for (int k=0;k<this->nFluidOutlet();++k)
        {
            if ( std::get<1>( M_fluidOutletsBCType[k] ) == "windkessel" &&  std::get<0>( std::get<2>( M_fluidOutletsBCType[k] ) ) == "implicit" )
                markerNameBFOutletForSubmesh.push_back( std::get<0>( M_fluidOutletsBCType[k] ) );
        }

        M_fluidOutletWindkesselMesh = createSubmesh( _mesh=this->mesh(), _range=markedfaces(this->mesh(),markerNameBFOutletForSubmesh), _view=true );
        M_fluidOutletWindkesselSpace = space_fluidoutlet_windkessel_type::New( _mesh=M_fluidOutletWindkesselMesh,
                                                                               _worldscomm=makeWorldsComm(2,this->worldCommPtr()) );
    }

    // clean
    M_fluidOutletWindkesselPressureDistal.clear();
    M_fluidOutletWindkesselPressureProximal.clear();
    M_fluidOutletWindkesselPressureDistal_old.clear();

    // windkessel outlet
    if ( this->hasFluidOutletWindkessel() )
    {
        for (int k=0;k<this->nFluidOutlet();++k)
        {
            if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" ) continue;

            // init containers
            M_fluidOutletWindkesselPressureDistal[k] = 0;
            M_fluidOutletWindkesselPressureProximal[k] = 0;
            M_fluidOutletWindkesselPressureDistal_old[k].resize( Feel::BDF_MAX_ORDER, 0 );
        }

        std::string nameFile = this->rootRepository() + "/" + prefixvm(this->prefix(),"fluidoutletbc.windkessel.data");

        if (!this->doRestart())
        {
            if (this->worldComm().isMasterRank())
            {
                std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::trunc);
                file.precision( 8 );
                file.setf( std::ios::scientific );
                file.width( 15 );
                file.setf( std::ios::left );
                file << int(0);
                file.width( 20 );
                file << this->timeInitial();

                for (int k=0;k<this->nFluidOutlet();++k)
                {
                    if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" ) continue;
                    // write value on disk
                    file.width( 20 );
                    file << M_fluidOutletWindkesselPressureDistal.find(k)->second;
                    file.width( 20 );
                    file << M_fluidOutletWindkesselPressureProximal.find(k)->second;
                }
                file << "\n";
                file.close();
            }
        }
        else
        {
            if (this->worldComm().isMasterRank())
            {
                std::ifstream fileI(nameFile.c_str(), std::ios::in);
                int cptIter=0; double timeIter=0;double valPresDistal=0,valPresProximal=0;
                int askedIter = M_bdfVelocity->iteration() - 1;
                bool find=false; std::ostringstream buffer;
                buffer.precision( 8 );
                buffer.setf( std::ios::scientific );

                while ( !fileI.eof() && !find )
                {
                    fileI >> cptIter >> timeIter;
                    buffer.width( 15 );
                    buffer.setf( std::ios::left );
                    buffer << cptIter;
                    buffer.width( 20 );
                    buffer << timeIter;

                    for (int k=0;k<this->nFluidOutlet();++k)
                    {
                        if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" ) continue;
                        fileI >> valPresDistal >> valPresProximal;
                        buffer.width( 20 );
                        buffer << valPresDistal;
                        buffer.width( 20 );
                        buffer << valPresProximal;

                        for (int l=0 ; l< M_fluidOutletWindkesselPressureDistal_old.find(k)->second.size() ; ++l)
                            if (cptIter == askedIter - l) M_fluidOutletWindkesselPressureDistal_old[k][l] = valPresDistal;

                        if (cptIter == askedIter)
                        {
                            M_fluidOutletWindkesselPressureDistal[k] = valPresDistal;
                            M_fluidOutletWindkesselPressureProximal[k] = valPresProximal;
                        }
                    }
                    buffer << "\n";

                    if (cptIter == askedIter) find=true;
                }
                fileI.close();
                std::ofstream fileW(nameFile.c_str(), std::ios::out | std::ios::trunc);
                fileW << buffer.str();
                fileW.close();

                //std::cout << cptIter <<" " << timeIter << " " << valPresDistal << std::endl;
                //M_fluidOutletWindkesselPressureDistal_old[k][0] = valPresDistal;
            }

            // broadcast windkessel data
            if ( this->worldComm().globalSize() > 1 )
            {
                auto dataToBroadcast = boost::make_tuple(M_fluidOutletWindkesselPressureDistal,
                                                         M_fluidOutletWindkesselPressureProximal,
                                                         M_fluidOutletWindkesselPressureDistal_old );
                mpi::broadcast( this->worldComm().globalComm(), dataToBroadcast, this->worldComm().masterRank() );
                if ( !this->worldComm().isMasterRank() )
                {
                    M_fluidOutletWindkesselPressureDistal = boost::get<0>( dataToBroadcast );
                    M_fluidOutletWindkesselPressureProximal = boost::get<1>( dataToBroadcast );
                    M_fluidOutletWindkesselPressureDistal_old = boost::get<2>( dataToBroadcast );
                }
                this->log("FluidMechanics","init", "restart windkessel broadcast done" );
            }
        }
        if (this->verbose())
        {
            std::ostringstream bufferLog;
            bufferLog << "\n";
            for (int k=0;k<this->nFluidOutlet();++k)
            {
                if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" ) continue;
                bufferLog << " (" << M_fluidOutletWindkesselPressureDistal[k] <<","<< M_fluidOutletWindkesselPressureProximal[k] << ")";
                bufferLog << " [";
                const int sizeOld = M_fluidOutletWindkesselPressureDistal_old.find(k)->second.size();
                bool hasDoneFirstElt = false;
                for (int l=0;l<sizeOld;++l)
                {
                    if ( hasDoneFirstElt ) bufferLog << " , ";
                    bufferLog << M_fluidOutletWindkesselPressureDistal_old[k][l];
                    hasDoneFirstElt = true;
                }
                bufferLog << "]\n";
            }
            this->log("FluidMechanics","windkessel initialisation", bufferLog.str() );
        }

        this->log("FluidMechanics","init", "restart windkessel done" );

    } // if (this->hasFluidOutletWindkessel())

    this->log("FluidMechanics","initFluidOutlet", "finish" );

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initUserFunctions()
{
    for ( auto const& modelfunc : this->modelProperties().functions() )
    {
        auto const& funcData = modelfunc.second;
        std::string funcName = funcData.name();

        if ( funcData.isScalar() )
        {
            if ( this->hasFieldUserScalar( funcName ) )
                continue;
            M_fieldsUserScalar[funcName] = this->functionSpaceVelocity()->compSpace()->elementPtr();
        }
        else if ( funcData.isVectorial2() )
        {
            if ( nDim != 2 ) continue;
            if ( this->hasFieldUserVectorial( funcName ) )
                continue;
            M_fieldsUserVectorial[funcName] = this->functionSpaceVelocity()->elementPtr();
        }
        else if ( funcData.isVectorial3() )
        {
            if ( nDim != 3 ) continue;
            if ( this->hasFieldUserVectorial( funcName ) )
                continue;
            M_fieldsUserVectorial[funcName] = this->functionSpaceVelocity()->elementPtr();
        }
    }

    // update custom field given by registerCustomField
    for ( auto & [name,uptr] : M_fieldsUserScalar )
        if ( !uptr )
            uptr = this->functionSpaceVelocity()->compSpace()->elementPtr();
    for ( auto & [name,uptr] : M_fieldsUserVectorial )
        if ( !uptr )
            uptr = this->functionSpaceVelocity()->elementPtr();

    this->updateUserFunctions();
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateUserFunctions( bool onlyExprWithTimeSymbol )
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
            M_fieldsUserScalar[funcName]->on(_range=M_rangeMeshElements,_expr=funcData.expressionScalar() );
        }
        else if ( funcData.isVectorial2() )
        {
            if constexpr( nDim == 2 )
            {
                CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
                M_fieldsUserVectorial[funcName]->on(_range=M_rangeMeshElements,_expr=funcData.expressionVectorial2() );
            }
            else CHECK( false ) << "TODO";
        }
        else if ( funcData.isVectorial3() )
        {
            if constexpr( nDim == 3 )
            {
                CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
                M_fieldsUserVectorial[funcName]->on(_range=M_rangeMeshElements,_expr=funcData.expressionVectorial3() );
            }
            else CHECK( false ) << "TODO";
        }
    }
}

#if 0
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::set<std::string>
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::postProcessFieldExported( std::set<std::string> const& ifields, std::string const& prefix ) const
{
    std::set<std::string> res;
    for ( auto const& o : ifields )
    {
        if ( o == prefixvm(prefix,"velocity") || o == prefixvm(prefix,"all") )
            res.insert( "velocity" );
        if ( o == prefixvm(prefix,"pressure") || o == prefixvm(prefix,"all") )
            res.insert( "pressure" );
        if ( o == prefixvm(prefix,"vorticity") || o == prefixvm(prefix,"all") )
            res.insert( "vorticity" );
        if ( o == prefixvm(prefix,"density") || o == prefixvm(prefix,"all") )
            res.insert( "density" );
        if ( o == prefixvm(prefix,"viscosity") || o == prefixvm(prefix,"all") )
            res.insert( "viscosity" );
        if ( o == prefixvm(prefix,"pid") || o == prefixvm(prefix,"all") )
            res.insert( "pid" );

        if ( o == prefixvm(prefix,"pressurebc") || o == prefixvm(prefix,"all") )
            res.insert( "pressurebc" );

        if ( this->isMoveDomain() )
        {
            if ( o == prefixvm(prefix,"displacement") || o == prefixvm(prefix,"all") )
                res.insert( "displacement" );
            if ( o == prefixvm(prefix,"alemesh") || o == prefixvm(prefix,"all") )
                res.insert( "alemesh" );
        }

        // add user functions
        if ( this->hasFieldUserScalar( o ) || this->hasFieldUserVectorial( o ) )
            res.insert( o );
    }
    return res;
}
#endif

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->setPostProcessExportsAllFieldsAvailable( {"velocity","pressure","vorticity","displacement","alemesh"} );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    this->setPostProcessExportsAllFieldsAvailable( "trace_mesh", {"trace.normal-stress","trace.wall-shear-stress" /*, "trace.body.translational-velocity", "trace.body.angular-velocity"*/ } );
    for ( auto const& [bname,bbc] : M_bodySetBC )
    {
        this->addPostProcessExportsAllFieldsAvailable( "trace_mesh", { (boost::format("body.%1%.translational-velocity")%bname).str(), (boost::format("body.%1%.angular-velocity")%bname).str() } );
        this->addPostProcessMeasuresQuantitiesAllNamesAvailable( { (boost::format("body_%1%.mass_center")%bname).str(),
                                                                   (boost::format("body_%1%.rigid_rotation_angles")%bname).str(),
                                                                   (boost::format("body_%1%.moment_of_inertia")%bname).str(),
                                                                   (boost::format("body_%1%.moment_of_inertia_body_frame")%bname).str(),
                                                                   (boost::format("body_%1%.fluid_forces")%bname).str(), (boost::format("body_%1%.fluid_torques")%bname).str() } );
    }
    this->setPostProcessExportsPidName( "trace_mesh", "trace.pid" );
    this->setPostProcessSaveAllFieldsAvailable( {"velocity","pressure","vorticity","displacement"} );
    super_type::initPostProcess();

    // init exporters
    if ( boption(_name="exporter.export") )
    {
        if ( !this->postProcessExportsFields().empty() )
        {
            this->createPostProcessExporters();
            // restart exporters if restart is activated
            if ( this->doRestart() && this->restartPath().empty() )
            {
                // if restart and same directory, update the exporter for new value, else nothing (create a new exporter)
                if ( M_exporter && M_exporter->doExport() )
                    M_exporter->restart( this->timeInitial() );
                if ( M_exporter_ho && M_exporter_ho->doExport() )
                    M_exporter_ho->restart( this->timeInitial() );
            }
        }

         if ( !this->postProcessExportsFields( "trace_mesh" ).empty() && nOrderGeo <= 2  )
        {
#if 1
            auto velocityMeshSupport = this->functionSpaceVelocity()->template meshSupport<0>();
            auto rangeTrace = velocityMeshSupport->rangeBoundaryFaces(); // not very nice, need to store the meshsupport
            M_meshTrace = createSubmesh( _mesh=velocityMeshSupport/*this->mesh()*/, _range=rangeTrace, _context=size_type(EXTRACTION_KEEP_MESH_RELATION|EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT),_view=true );
#else
            auto rangeTrace = M_bodySetBC.begin()->second.rangeMarkedFacesOnFluid();
            M_meshTrace = M_bodySetBC.begin()->second.mesh();
#endif
            this->updateRangeDistributionByMaterialName( "trace_mesh", rangeTrace );
            std::string geoExportType = "static";//change_coords_only, change, static
            M_exporterTrace = exporter( _mesh=M_meshTrace,
                                        _name="Export_trace",
                                        _geo=geoExportType,
                                        _worldcomm=this->worldComm(),
                                        _path=this->exporterPath() );
            if ( this->doRestart() && this->restartPath().empty() )
            {
                // if restart and same directory, update the exporter for new value, else nothing (create a new exporter)
                if ( M_exporterTrace && M_exporterTrace->doExport() )
                    M_exporterTrace->restart( this->timeInitial() );
            }
        }
    }


    // forces (lift, drag) and flow rate measures
    pt::ptree ptree = this->modelProperties().postProcess().pTree( this->keyword() );
    std::string ppTypeMeasures = "Measures";
    for( auto const& ptreeLevel0 : ptree )
    {
        std::string ptreeLevel0Name = ptreeLevel0.first;
        if ( ptreeLevel0Name != ppTypeMeasures ) continue;
        for( auto const& ptreeLevel1 : ptreeLevel0.second )
        {
            std::string ptreeLevel1Name = ptreeLevel1.first;
            if ( ptreeLevel1Name == "Forces" )
            {
                // get list of marker
                std::set<std::string> markerSet;
                std::string markerUnique = ptreeLevel1.second.template get_value<std::string>();
                if ( markerUnique.empty() )
                {
                    for (auto const& ptreeMarker : ptreeLevel1.second )
                    {
                        std::string marker = ptreeMarker.second.template get_value<std::string>();
                        markerSet.insert( marker );
                    }
                }
                else
                {
                    markerSet.insert( markerUnique );
                }
                // save forces measure for each marker
                for ( std::string const& marker : markerSet )
                {
                    ModelMeasuresForces myPpForces;
                    myPpForces.addMarker( marker );
                    myPpForces.setName( marker );
                    std::string name = myPpForces.name();
                    M_postProcessMeasuresForces.push_back( myPpForces );
                }
            }
            else if ( ptreeLevel1Name == "FlowRate" )
            {
                for( auto const& ptreeLevel2 : ptreeLevel1.second )
                {
                    ModelMeasuresFlowRate myPpFlowRate;
                    std::string name = ptreeLevel2.first;
                    myPpFlowRate.setup( ptreeLevel2.second, name );
                    M_postProcessMeasuresFlowRate.push_back( myPpFlowRate );
                }
            }
            else if ( ptreeLevel1Name == "Pressure" )
            {
                // this->modelProperties().postProcess().operator[](ppTypeMeasures).push_back( "Pressure" );
                M_postProcessMeasuresFields["pressure"] = "";
            }
            else if ( ptreeLevel1Name == "VelocityDivergence" )
            {
                M_postProcessMeasuresFields["velocity-divergence"] = "";
            }
        }
    }

    // point measures
    auto fieldNamesWithSpaceVelocity = std::make_pair( std::set<std::string>({"velocity"}), this->functionSpaceVelocity() );
    auto fieldNamesWithSpacePressure = std::make_pair( std::set<std::string>({"pressure"}), this->functionSpacePressure() );
    auto fieldNamesWithSpaces = hana::make_tuple( fieldNamesWithSpaceVelocity, fieldNamesWithSpacePressure );
    M_measurePointsEvaluation = std::make_shared<measure_points_evaluation_type>( fieldNamesWithSpaces );
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( this->keyword() ) )
    {
       M_measurePointsEvaluation->init( evalPoints );
    }


    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasuresIO().restart( "time", this->timeInitial() );
        else
            this->postProcessMeasuresIO().setMeasure( "time", this->timeInitial() ); //just for have time in first column
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::size_type
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initStartBlockIndexFieldsInMatrix()
{
    size_type currentStartIndex = 0;
    this->setStartSubBlockSpaceIndex( "velocity", currentStartIndex++ );
    this->setStartSubBlockSpaceIndex( "pressure", currentStartIndex++ );
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        this->setStartSubBlockSpaceIndex( "define-pressure-cst-lm", currentStartIndex );
        currentStartIndex += M_XhMeanPressureLM.size();
    }
    if (this->hasMarkerDirichletBClm())
    {
        this->setStartSubBlockSpaceIndex( "dirichletlm", currentStartIndex++ );
    }
    if ( this->hasMarkerPressureBC() )
    {
        this->setStartSubBlockSpaceIndex( "pressurelm1", currentStartIndex++ );
        if ( nDim == 3 )
            this->setStartSubBlockSpaceIndex( "pressurelm2", currentStartIndex++ );
    }
    if ( this->hasFluidOutletWindkesselImplicit() )
    {
        this->setStartSubBlockSpaceIndex( "windkessel", currentStartIndex++ );
    }

    for ( auto const& [bpname,bbc] : M_bodySetBC )
    {
        this->setStartSubBlockSpaceIndex( "body-bc."+bbc.name()+".translational-velocity", currentStartIndex++ );
        if ( !bbc.isInNBodyArticulated() )
            this->setStartSubBlockSpaceIndex( "body-bc."+bbc.name()+".angular-velocity", currentStartIndex++ );
    }
    for ( auto const& nba : M_bodySetBC.nbodyArticulated() )
    {
        this->setStartSubBlockSpaceIndex( "body-bc."+nba.name()+".angular-velocity", currentStartIndex );
        // defined alias in each body
        for ( auto bbc : nba.bodyList() )
            this->setStartSubBlockSpaceIndex( "body-bc."+bbc->name()+".angular-velocity", currentStartIndex );
        ++currentStartIndex;

        if ( nba.articulationMethod() != "lm" )
            continue;
        for ( auto const& ba : nba.articulations() )
            this->setStartSubBlockSpaceIndex( "body-bc.articulation-lm."+ba.name()+".translational-velocity", currentStartIndex++ );
    }


    return currentStartIndex;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::buildBlockVector()
{
    this->initBlockVector();
    this->algebraicBlockVectorSolution()->buildVector( this->backend() );

    M_usePreviousSolution = false;
    if ( timeStepping() == "Theta" && this->stabilizationGLS() )
        M_usePreviousSolution = true;
    if  ( !M_usePreviousSolution && this->hasTurbulenceModel() )
    {
        for ( auto const& cfpde : M_turbulenceModelType->coefficientFormPDEs() )
        {
            if ( cfpde->timeStepping() == "Theta" && cfpde->applyStabilization() )
            {
                M_usePreviousSolution = true;
                break;
            }
        }
    }

    if ( M_usePreviousSolution )
        M_vectorPreviousSolution = this->backend()->newVector( this->algebraicBlockVectorSolution()->vectorMonolithic()->mapPtr() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
int
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initBlockVector()
{
    int nBlock = this->nBlockMatrixGraph();
    auto bvs = this->initAlgebraicBlockVectorSolution( nBlock );
    int cptBlock = 0;
    bvs->operator()(cptBlock++) = this->fieldVelocityPtr();
    bvs->operator()(cptBlock++) = this->fieldPressurePtr();
    // impose mean pressure by lagrange multiplier
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        for ( int k=0;k<M_XhMeanPressureLM.size();++k )
            bvs->operator()(cptBlock++) = this->backend()->newVector( M_XhMeanPressureLM[k] );
    }
    // lagrange multiplier for Dirichlet BC
    if (this->hasMarkerDirichletBClm())
    {
        bvs->operator()(cptBlock) = M_fieldDirichletLM;//this->backend()->newVector( this->XhDirichletLM() );
        ++cptBlock;
    }
    if ( this->hasMarkerPressureBC() )
    {
        bvs->operator()(cptBlock++) = M_fieldLagrangeMultiplierPressureBC1;
        if ( nDim == 3 )
            bvs->operator()(cptBlock++) = M_fieldLagrangeMultiplierPressureBC2;
    }
    // windkessel outel with implicit scheme
    if ( this->hasFluidOutletWindkesselImplicit() )
    {
        for (int k=0;k<this->nFluidOutletWindkesselImplicit();++k)
        {
            bvs->operator()(cptBlock) = this->backend()->newVector( M_fluidOutletWindkesselSpace );
            ++cptBlock;
        }
    }

    if ( !M_bodySetBC.empty() )
    {
        for ( auto const& [bpname,bpbc] : M_bodySetBC )
        {
            bvs->operator()(cptBlock++) = bpbc.fieldTranslationalVelocityPtr();
            if ( !bpbc.isInNBodyArticulated() )
                bvs->operator()(cptBlock++) = bpbc.fieldAngularVelocityPtr();
        }
        for ( auto const& nba : M_bodySetBC.nbodyArticulated() )
        {
            bvs->operator()(cptBlock++) = nba.fieldAngularVelocityPtr();
            if ( nba.articulationMethod() != "lm" )
                continue;
            for ( auto const& ba : nba.articulations() )
                bvs->operator()(cptBlock++) = ba.vectorLagrangeMultiplierTranslationalVelocity();
        }
    }

    return cptBlock;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initInHousePreconditioner()
{


    if ( M_preconditionerAttachPCD )
    {
        typedef Feel::Alternatives::OperatorPCD<space_velocity_type,space_pressure_type> op_pcd_type;
        auto opPCD = std::make_shared<op_pcd_type>( this->functionSpaceVelocity(), this->functionSpacePressure(),
                                                    this->backend(), this->prefix(), true);

        for( auto const& d : M_bcDirichlet )
        {
            std::set<std::string> themarkers;
            for ( std::string const& type : std::vector<std::string>( { "elimination", "nitsche", "lm" } ) )
            {
                auto ret = detail::distributeMarkerListOnSubEntity( this->mesh(),this->markerDirichletBCByNameId( type, name(d) ) );
                themarkers.insert( std::get<0>( ret ).begin(), std::get<0>( ret ).end() );
            }
            opPCD->addRangeDirichletBC( name(d), markedfaces( this->mesh(), themarkers ) );
        }
        for( auto const& d : M_bcMovingBoundaryImposed )
        {
            std::set<std::string> themarkers;
            for ( std::string const& type : std::vector<std::string>( { "elimination", "nitsche", "lm" } ) )
            {
                auto ret = detail::distributeMarkerListOnSubEntity( this->mesh(), M_bcMarkersMovingBoundaryImposed.markerDirichletBCByNameId( type, name(d) ) );
                themarkers.insert( std::get<0>( ret ).begin(), std::get<0>( ret ).end() );
            }
            opPCD->addRangeDirichletBC( name(d), markedfaces( this->mesh(), themarkers ) );
        }
        for ( auto const& inletbc : M_fluidInletDesc )
        {
            std::string const& themarker = std::get<0>( inletbc );
            opPCD->addRangeDirichletBC( themarker, markedfaces( this->mesh(), themarker ) ); // warning marker is the name
        }

        std::set<std::string> markersNeumann;
        for( auto const& d : M_bcNeumannScalar )
        {
            auto themarkers = this->markerNeumannBC(NeumannBCShape::SCALAR,name(d));
            markersNeumann.insert( themarkers.begin(), themarkers.end() );
        }
        for( auto const& d : M_bcNeumannVectorial )
        {
            auto themarkers = this->markerNeumannBC(NeumannBCShape::VECTORIAL,name(d));
            markersNeumann.insert( themarkers.begin(), themarkers.end() );
        }
        for( auto const& d : M_bcNeumannTensor2 )
        {
            auto themarkers = this->markerNeumannBC(NeumannBCShape::TENSOR2,name(d));
            markersNeumann.insert( themarkers.begin(), themarkers.end() );
        }
        for ( auto const& bcOutlet : M_fluidOutletsBCType )
        {
            markersNeumann.insert( std::get<0>(bcOutlet) );
        }
        opPCD->addRangeNeumannBC( "FluidNeumann", markedfaces( this->mesh(), markersNeumann ) );

        for ( auto const& f : M_addUpdateInHousePreconditionerPCD )
            f.second.first( *opPCD );

        opPCD->initialize();
        M_operatorPCD = opPCD;
    }

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initDist2Wall()
{
    space_dist2wall_ptrtype M_spaceDist2Wall;
    M_spaceDist2Wall = space_dist2wall_type::New(_mesh=this->mesh() );
    M_fieldDist2Wall = M_spaceDist2Wall->elementPtr();

    auto thefms = fms( M_spaceDist2Wall );

    //auto phio = Xh->element();
    //phio = vf::project(Xh, elements(mesh), h() );
    M_fieldDist2Wall->on(_range=elements(this->mesh()),_expr=h() );

    auto rangeWall = M_dist2WallMarkers.empty()? boundaryfaces(this->mesh()) : markedfaces(this->mesh(), M_dist2WallMarkers );

    (*M_fieldDist2Wall) +=vf::project(_space=M_spaceDist2Wall,
                                      _range=rangeWall,
                                      _expr= -idv(M_fieldDist2Wall) - h()/100. );
    *M_fieldDist2Wall = thefms->march(*M_fieldDist2Wall);
    M_fieldDist2Wall->on(_range=rangeWall,_expr=cst(0.),_close=true);
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initTurbulenceModel()
{
    M_turbulenceModelType.reset( new turbulence_model_type( prefixvm(this->prefix(),"turbulence"),
                                                            "cfpdes",
                                                            this->worldCommPtr(), "", this->repository() ) );

    bool isSpalartAllmarasTurbulenceModel = this->hasTurbulenceModel( "Spalart-Allmaras" );
    
    std::string eqkeyword;
    if ( isSpalartAllmarasTurbulenceModel )
    {
        /*std::string*/ eqkeyword = prefixvm( this->keyword(), "turbulence_SA", "_" );
        M_turbulenceModelType->addGenericPDE( typename FeelModels::ModelGenericPDE<nDim>::infos_type( eqkeyword, "solution", "nu", "Pch1" ) );
    }
    else
    {
        std::string eqkeyword_k = prefixvm( this->keyword(), "turbulence_k", "_" );
        std::string eqkeyword_epsilon = prefixvm( this->keyword(), "turbulence_epsilon", "_" );
        M_turbulenceModelType->addGenericPDE( typename FeelModels::ModelGenericPDE<nDim>::infos_type( eqkeyword_k, "solution_k", "k", "Pch1" ) );
        M_turbulenceModelType->addGenericPDE( typename FeelModels::ModelGenericPDE<nDim>::infos_type( eqkeyword_epsilon, "solution_epsilon", "epsilon", "Pch1" ) );
    }

    std::string symb_dist2wall = prefixvm( this->keyword(),"dist2wall","_" );
    std::string symb_curl_magintude = prefixvm( this->keyword(), "curl_U_magnitude", "_");
    std::string symb_velocity_magintude = prefixvm( this->keyword(), "U_magnitude", "_");
    std::string symb_strain_rate_magnitude = prefixvm( this->keyword(), "strain_rate_magnitude", "_");
    std::string symb_velocity_x =  prefixvm( this->keyword(),"U","_") +"_0";
    std::string symb_velocity_y =  prefixvm( this->keyword(),"U","_") +"_1";
    std::string symb_velocity_z =  prefixvm( this->keyword(),"U","_") +"_2";
    std::string exprstr_velocity = nDim==2? (boost::format("{%1%,%2%}:%1%:%2%")%symb_velocity_x %symb_velocity_y).str() :
        (boost::format("{%1%,%2%,%3%}:%1%:%2%:%3%")%symb_velocity_x %symb_velocity_y %symb_velocity_z).str() ;

    auto myCvrtSeqToStr = []( auto const& cont ) {
        std::string markerListStr;
        for ( std::string const& mark : cont )
        {
            if ( !markerListStr.empty() )
                markerListStr += ",";
            markerListStr += "\"" + mark + "\"";
        }
        return markerListStr;
    };


    std::ostringstream ostr;
    //std::vector<std::string> strInitialConditions;
    std::map<std::string,std::vector<std::string>> strInitialConditions;
    ostr << "{";

    ostr << "\"Materials\":{";
    bool isFirstMaterial = true;
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
        if ( !physicFluidData->turbulence().isEnabled() )
            continue;
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            //auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto & matProps = this->materialsProperties()->materialProperties( matName );

            std::string symbDensity = "materials_" + matName + "_rho";
            std::string symbDynViscosity = "materials_" + matName + "_mu";
            std::string symbTurbulentDynViscosity = "materials_" + matName + "_mu_t";

            if ( !isFirstMaterial )
            {
                ostr  << ",";
                isFirstMaterial = false;
            }
            ostr << "\""<<matName<<"\":{"
                 << "\"markers\":[";
            bool isFirstMark = true;
            for ( std::string const& m : matProps.markers() )
            {
                if ( !isFirstMark )
                    ostr << ",";
                ostr << "\"" << m << "\"";
                isFirstMark = false;
            }
            ostr << "],";

            if ( physicFluidData->turbulence().model() == "Spalart-Allmaras" )
            {
                // std::cout << "HOLA eqname : " << std::get<0>( M_turbulenceModelType->pdes()[0] ).equationName() <<  " versus " << eqkeyword << " and "
                //           << std::get<0>( M_turbulenceModelType->pdes()[0] ).unknownSymbol() << std::endl;
                std::string eqkeyword = std::get<0>( M_turbulenceModelType->pdes()[0] ).equationName();

                std::string symb_sol_SA = prefixvm( eqkeyword, "nu" , "_" );
                std::string symbbase_grad_sol_SA = prefixvm( eqkeyword, "grad_nu", "_" );
                std::string exprstrnosymb_inner_grad_sol_SA = nDim==2? (boost::format("%1%_0^2+%1%_1^2")%symbbase_grad_sol_SA ).str() :
                    (boost::format("(%1%_0^2+%1%_1^2+%1%_2^2)")%symbbase_grad_sol_SA ).str();
                std::string depsymb_inner_grad_sol_SA = nDim==2? (boost::format("%1%_0:%1%_1")%symbbase_grad_sol_SA ).str() :
                    (boost::format("%1%_0:%1%_1:%1%_2")%symbbase_grad_sol_SA ).str();


                std::string symb_c_b1 = physicFluidData->addParameter( prefixvm( eqkeyword, "c_b1", "_" ), 0.1355 );
                std::string symb_c_b2 = physicFluidData->addParameter( prefixvm( eqkeyword, "c_b2", "_" ), 0.622 );
                std::string symb_c_v1 = physicFluidData->addParameter( prefixvm( eqkeyword, "c_v1", "_" ), 7.1 );
                std::string symb_c_v2 = physicFluidData->addParameter( prefixvm( eqkeyword, "c_v2", "_" ), 0.7 );
                std::string symb_c_v3 = physicFluidData->addParameter( prefixvm( eqkeyword, "c_v3", "_" ), 0.9 );
                std::string symb_c_t3 = physicFluidData->addParameter( prefixvm( eqkeyword, "c_t3", "_" ), 1.2 );
                std::string symb_c_t4 = physicFluidData->addParameter( prefixvm( eqkeyword, "c_t4", "_" ), 0.5 );
                std::string symb_c_w2 = physicFluidData->addParameter( prefixvm( eqkeyword, "c_w2", "_" ), 0.3 );
                std::string symb_c_w3 = physicFluidData->addParameter( prefixvm( eqkeyword, "c_w3", "_" ), 2 );
                std::string symb_c_kappa = physicFluidData->addParameter( prefixvm( eqkeyword, "c_kappa", "_" ), 0.41 );
                std::string symb_c_sigma = physicFluidData->addParameter( prefixvm( eqkeyword, "c_sigma", "_" ), 2./3. );
                std::string symb_c_w1 = physicFluidData->addParameter( prefixvm( eqkeyword, "c_w1", "_" ), (boost::format("%1%/(%2%^2) + (1+%3%)/%4% :%1%:%2%:%3%:%4%")%symb_c_b1 %symb_c_kappa %symb_c_b2 %symb_c_sigma ).str(),
                                                                       this->worldComm(), this->repository().expr() );

                std::string symb_khi = physicFluidData->addParameter( prefixvm( eqkeyword, "khi", "_" ), (boost::format("%1%*%2%/%3%:%1%:%2%:%3%")%symb_sol_SA %symbDensity %symbDynViscosity ).str(),
                                                                      this->worldComm(), this->repository().expr() );
                std::string symb_f_v1 = physicFluidData->addParameter( prefixvm( eqkeyword, "f_v1", "_" ), (boost::format("%1%^3/(%1%^3+%2%^3):%1%:%2%")%symb_khi %symb_c_v1 ).str(), this->worldComm(), this->repository().expr() );
                std::string symb_f_v2 = physicFluidData->addParameter( prefixvm( eqkeyword, "f_v2", "_" ), (boost::format("1 - %1%/(1+%1%*%2%) :%1%:%2%")%symb_khi %symb_f_v1 ).str(), this->worldComm(), this->repository().expr() );
                std::string symb_f_t2 = physicFluidData->addParameter( prefixvm( eqkeyword, "f_t2", "_" ), (boost::format("%1%*exp(-%2%*%3%^2) :%1%:%2%:%3%")%symb_c_t3 %symb_c_t4 %symb_khi ).str(),
                                                                       this->worldComm(), this->repository().expr() );

                std::string symb_S_bar =  physicFluidData->addParameter( prefixvm( eqkeyword, "S_bar", "_" ), (boost::format( "(%1%/(%2%^2*%3%^2))*%4% :%1%:%2%:%3%:%4%" )%symb_sol_SA  %symb_c_kappa %symb_dist2wall % symb_f_v2 ).str(),
                                                                         this->worldComm(), this->repository().expr() );
#if 0
                // classic
                std::string symb_S = physicFluidData->addParameter( prefixvm( eqkeyword, "S", "_" ),
                                                                    (boost::format("%1% + (%2%/(%3%^2*%4%^2))*%5% :%1%:%2%:%3%:%4%:%5%")%symb_curl_magintude %symb_sol_SA  %symb_c_kappa %symb_dist2wall %symb_f_v2 ).str(),
                                                                    this->worldComm(), this->repository().expr() );
#elif 0
                // variant note(c)
                std::string symb_S = physicFluidData->addParameter( prefixvm( eqkeyword, "S", "_" ), (boost::format("%1% + (1-(%2% < (-%3%*%1%)))*%2% + (%2% < (-%3%*%1%))*( %1%*(%3%^2*%1%+%4%*%2%)/((%4%-2*%3%)*%1% - %2%) ):%1%:%2%:%3%:%4%")%symb_curl_magintude %symb_S_bar %symb_c_v2 %symb_c_v3 ).str(), this->worldComm(), this->repository().expr() );
#else
                // variant note1b +
                std::string symb_c_rot = physicFluidData->addParameter( prefixvm( eqkeyword, "c_rot", "_" ), 2.0 );
                std::string symb_S = physicFluidData->addParameter( prefixvm( eqkeyword, "S", "_" ),
                                                                    (boost::format("max( %1% + %2%*min(0,%3%-%1%) + %4%, 0.3*%1% ):%1%:%2%:%3%:%4%")%symb_curl_magintude %symb_c_rot %symb_strain_rate_magnitude %symb_S_bar  ).str(),
                                                                    this->worldComm(), this->repository().expr() );
#endif

                std::string symb_r = physicFluidData->addParameter( prefixvm( eqkeyword, "r", "_" ),
                                                                    //(boost::format("(%2%>0)*min(%1%/(%2%*%3%^2*%4%^2),10) + (1-(%2%>0))*10:%1%:%2%:%3%:%4%")%symb_sol_SA %symb_S %symb_c_kappa %symb_dist2wall ).str(),
                                                                    //(boost::format("min(%1%/(%2%*%3%^2*%4%^2 + (1-(%2%>0))*1e-5 ),10):%1%:%2%:%3%:%4%")%symb_sol_SA %symb_S %symb_c_kappa %symb_dist2wall ).str(),
                                                                    (boost::format("min(%1%/(%2%*%3%^2*%4%^2 + (1-(%2%>1e-7))*1e-5 ),10):%1%:%2%:%3%:%4%")%symb_sol_SA %symb_S %symb_c_kappa %symb_dist2wall ).str(),
                                                                    this->worldComm(), this->repository().expr() );
                std::string symb_g = physicFluidData->addParameter( prefixvm( eqkeyword, "g", "_" ), (boost::format("%1% + %2%*(%1%^6 - %1%):%1%:%2%")%symb_r %symb_c_w2 ).str(), this->worldComm(), this->repository().expr() );
                std::string symb_f_w = physicFluidData->addParameter( prefixvm( eqkeyword, "f_w", "_" ), (boost::format("%1%*( (1+%2%^6)/(%1%^6+%2%^6) )^(1/6):%1%:%2%")%symb_g %symb_c_w3 ).str(), this->worldComm(), this->repository().expr() );

                std::string exprstr_diffusion = (boost::format("(1/%1%)*(%2%/%3%+%4%):%1%:%2%:%3%:%4%:") %symb_c_sigma %symbDynViscosity %symbDensity %symb_sol_SA ).str();
#if 1
                std::string exprstr_reaction = (boost::format("-%1%*(1-%2%)*%3% + (%4%*%5% - %1%*%2%/(%6%^2))*(%7%/(%8%^2)) :%1%:%2%:%3%:%4%:%5%:%6%:%7%:%8%")%symb_c_b1 %symb_f_t2 %symb_S %symb_c_w1 %symb_f_w %symb_c_kappa %symb_sol_SA %symb_dist2wall ).str();
                std::string exprstr_source = (boost::format("(1/%1%)*%2%*%3%:%1%:%2%:%4%") %symb_c_sigma %symb_c_b2 %exprstrnosymb_inner_grad_sol_SA %depsymb_inner_grad_sol_SA ).str();
#else
                std::string exprstr_source = (boost::format("(1/%1%)*%2%*%3% -%11%*( -%5%*(1-%6%)*%7% + (%8%*%9% - %5%*%6%/(%10%^2))*(%11%/(%12%^2))  )  :%1%:%2%:%4%:%5%:%6%:%7%:%8%:%9%:%10%:%11%:%12%") %symb_c_sigma %symb_c_b2 %exprstrnosymb_inner_grad_sol_SA %depsymb_inner_grad_sol_SA      %symb_c_b1 %symb_f_t2 %symb_S %symb_c_w1 %symb_f_w %symb_c_kappa %symb_sol_SA %symb_dist2wall).str();
#endif
                std::string exprstr_convection = nDim==2? (boost::format("{%1% -(%3%/%4%)*%5%_0 ,%2% -(%3%/%4%)*%5%_1  }:%1%:%2%:%3%:%4%:%5%_0:%5%_1")%symb_velocity_x %symb_velocity_y %symb_c_b2 %symb_c_sigma %symbbase_grad_sol_SA).str() :
                    (boost::format("{%1%,%2%,%3%}:%1%:%2%:%3%")%symb_velocity_x %symb_velocity_y %symb_velocity_z).str() ; // NOT FINISH IN 3D

                ostr << "\""<<eqkeyword << "_beta\":\"" << /*exprstr_convection*/exprstr_velocity << "\","
                     << "\""<<eqkeyword << "_c\":\"" << exprstr_diffusion << "\","
                     << "\""<<eqkeyword << "_a\":\"" << exprstr_reaction << "\","
                     << "\""<<eqkeyword << "_f\":\"" << exprstr_source << "\"";
                if ( !this->isStationaryModel() )
                    ostr << ",\""<<eqkeyword << "_d\":\"1\"";
                //ostr << "}"; // end mymat1

                std::string mutExprStr = (boost::format("(%2%>0)*%1%*%2%*%3%:%1%:%2%:%3%")%symbDensity %symb_sol_SA %symb_f_v1 ).str();
                ModelExpression mutExpr;
                mutExpr.setExpr( mutExprStr, this->worldComm(), this->repository().expr() );
                this->materialsProperties()->addProperty( matProps, "turbulent-dynamic-viscosity", mutExpr, true );

                // initials coondition
                std::ostringstream ostr_ic;
                ostr_ic << "\"markers\":[" << myCvrtSeqToStr( matProps.markers() ) <<  "],"
                        << "\"expr\":\"" << (boost::format("%1%/%2%:%1%:%2%")%symbDynViscosity %symbDensity ).str() << "\"";
                strInitialConditions["solution"].push_back( ostr_ic.str() );
            }
            else if ( physicFluidData->turbulence().model() == "k-epsilon" )
            {
                std::string eqkeyword_k = std::get<0>( M_turbulenceModelType->pdes()[0] ).equationName();
                std::string unknownName_k = std::get<0>( M_turbulenceModelType->pdes()[0] ).unknownName();
                std::string symb_sol_k_base = std::get<0>( M_turbulenceModelType->pdes()[0] ).unknownSymbol();
                std::string symb_sol_k = prefixvm( eqkeyword_k, symb_sol_k_base , "_" );
                std::string symb_sol_k_previous = prefixvm( eqkeyword_k, symb_sol_k_base + "_previous" , "_" );
                std::string eqkeyword_epsilon = std::get<0>( M_turbulenceModelType->pdes()[1] ).equationName();
                std::string unknownName_epsilon = std::get<0>( M_turbulenceModelType->pdes()[1] ).unknownName();
                std::string symb_sol_epsilon_base = std::get<0>( M_turbulenceModelType->pdes()[1] ).unknownSymbol();
                std::string symb_sol_epsilon = prefixvm( eqkeyword_epsilon, symb_sol_epsilon_base , "_" );
                std::string symb_sol_epsilon_previous = prefixvm( eqkeyword_epsilon, symb_sol_epsilon_base + "_previous" , "_" );

                std::string prefix_symbol_physic_nomat = prefixvm( this->keyword(), "turbulence_k_epsilon", "_" );
                std::string prefix_symbol_physic_mat = prefixvm( this->keyword(), prefixvm( matName, "turbulence_k_epsilon" , "_" ), "_" );

                std::string symb_c_1epsilon = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_nomat, "c_1epsilon", "_" ), 1.44 );
                std::string symb_c_2epsilon = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_nomat, "c_2epsilon", "_" ), 1.92 );
                std::string symb_c_mu = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_nomat, "c_mu", "_" ), 0.09 );
                std::string symb_sigma_k = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_nomat, "sigma_k", "_" ), 1. );
                std::string symb_sigma_epsilon = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_nomat, "sigma_epsilon", "_" ), 1.3 );

                std::string symb_sol_k_positive = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "sol_k_positive", "_" ),
                                                                        (boost::format("max(%1%,0):%1%") %symb_sol_k ).str(),
                                                                        this->worldComm(), this->repository().expr() );
                std::string symb_sol_epsilon_positive = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "sol_epsilon_positive", "_" ),
                                                                        (boost::format("max(%1%,0):%1%") %symb_sol_epsilon ).str(),
                                                                        this->worldComm(), this->repository().expr() );

                std::string symb_sol_k_previous_positive = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "sol_k_previous_positive", "_" ),
                                                                        (boost::format("max(%1%,0):%1%") %symb_sol_k_previous ).str(),
                                                                        this->worldComm(), this->repository().expr() );
                std::string symb_sol_epsilon_previous_positive = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "sol_epsilon_previous_positive", "_" ),
                                                                                                (boost::format("max(%1%,0):%1%") %symb_sol_epsilon_previous ).str(),
                                                                                                this->worldComm(), this->repository().expr() );


                std::string symb_l_max = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_nomat, "l_max", "_" ), 0.0635/2. );
                std::string symb_l_star = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "l_star", "_" ),
                                                                         //(boost::format("(%1%*%2%^(3/2) < %3%*%4%)*(%1%*%2%^(3/2))/%3% + (1-(%1%*%2%^(3/2) < %3%*%4%))*%4%  :%1%:%2%:%3%:%4%")%symb_c_mu %symb_sol_k_positive %symb_sol_epsilon_positive %symb_l_max ).str(),
                                                                         //(boost::format("(%1%*%2%^(3/2) < %3%*%4%)*(%1%*%2%^(3/2))/(max(%3%,1e-16)) + (1-(%1%*%2%^(3/2) < %3%*%4%))*%4%  :%1%:%2%:%3%:%4%")%symb_c_mu %symb_sol_k_positive %symb_sol_epsilon_positive %symb_l_max ).str(),
                                                   (boost::format("min( %1%*%2%^(3/2)/(max(%3%,1e-16)),%4% ):%1%:%2%:%3%:%4%")%symb_c_mu %symb_sol_k_positive %symb_sol_epsilon_positive %symb_l_max ).str(),
                                                                         this->worldComm(), this->repository().expr() );


                std::string symb_toto_chi = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "toto_chi", "_" ),
                                                                           (boost::format("(0.09*%1%^(1.5) < %2%*0.0635) :%1%:%2%") %symb_sol_k_positive %symb_sol_epsilon_positive ).str(),
                                                                           this->worldComm(), this->repository().expr() );

                std::string symb_l_starBIS = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "l_starBIS", "_" ),
                                             //(boost::format("(%1%*%2%^(3/2) < %3%*%4%)*(%1%*%2%^(3/2))/%3% + (1-(%1%*%2%^(3/2) < %3%*%4%))*%4%  :%1%:%2%:%3%:%4%")%symb_c_mu %symb_sol_k_previous_positive %symb_sol_epsilon_previous_positive %symb_l_max ).str(),
                                                                            // (boost::format("(%1%*%2%^(3/2) < %3%*%4%)*(%1%*%2%^(3/2))/(max(%3%,1e-16)) + (1-(%1%*%2%^(3/2) < %3%*%4%))*%4%  :%1%:%2%:%3%:%4%")%symb_c_mu %symb_sol_k_previous_positive %symb_sol_epsilon_previous_positive %symb_l_max ).str(),
                                     (boost::format("min( %1%*%2%^(3/2)/(max(%3%,1e-16)),%4% ):%1%:%2%:%3%:%4%")%symb_c_mu %symb_sol_k_previous_positive %symb_sol_epsilon_previous_positive %symb_l_max ).str(),
                                                                         this->worldComm(), this->repository().expr() );

                std::string symbTurbulentDynViscosityBIS =  physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "mu_tBIS", "_" ),
                                                                                           //(boost::format("max( 1e-8, %1%*%2%*sqrt(%3%) ):%1%:%2%:%3%")%symbDensity %symb_l_starBIS %symb_sol_k_previous_positive ).str(),
                                                                                           (boost::format("max( 1e-8, %1%*%2%*sqrt(max(1e-9,%3%)) ):%1%:%2%:%3%")%symbDensity %symb_l_starBIS %symb_sol_k_previous_positive ).str(),
                                                                                           //(boost::format("max( 1e-8, %1%*%2% ):%1%:%2%")%symbDensity %symb_l_starBIS ).str(),
                                                                                           //(boost::format("max( 1e-8, %1%*sqrt(max(1e-9,%2%)) ):%1%:%2%")%symbDensity  %symb_sol_k_previous_positive ).str(),
                                                                                           this->worldComm(), this->repository().expr() );

#if 0
                std::string symbTurbulentDynViscosityBIS2 = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "mu_tBIS", "_" ),
                                                                                           (boost::format("max( 1e-8, %1%*sqrt(max(1e-9,%2%)) ):%1%:%2%")%symbDensity  %symb_sol_k_previous_positive ).str(),
                                                                                           this->worldComm(), this->repository().expr() );
#endif
                std::string symb_gamma = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "gamma", "_" ),
                                                                        //(boost::format("max(%1%*%2%/%3%,0):%1%:%2%:%3%")%symb_c_mu %symb_sol_k_previous_positive %symbTurbulentDynViscosityBIS2 ).str(),
                                                                        //(boost::format("max(%1%*%2%/%3%,1e-9):%1%:%2%:%3%")%symb_c_mu %symb_sol_k_previous_positive %symbTurbulentDynViscosityBIS ).str(),
                                                                        //(boost::format("(%2%>0)*%1%*%2%/%3% :%1%:%2%:%3%")%symb_c_mu %symb_sol_k_previous_positive %symbTurbulentDynViscosityBIS ).str(),
                                                                        (boost::format("(%2%>0)*%1%*%2%/%3% :%1%:%2%:%3%")%symb_c_mu %symb_sol_k_positive %symbTurbulentDynViscosityBIS ).str(), //AIE
                                                                        //(boost::format("max(%1%*%2%/%3%,0):%1%:%2%:%3%")%symb_c_mu %symb_sol_k_positive %symbTurbulentDynViscosityBIS ).str(),
                                                                        this->worldComm(), this->repository().expr() );

                std::string exprstr_convection = nDim==2? (boost::format("{%1%*%2%,%1%*%3%}:%1%:%2%:%3%")%symbDensity %symb_velocity_x %symb_velocity_y).str() :
                    (boost::format("{%1%*%2%,%1%*%3%,%1%*%4%}:%1%:%2%:%3%:%4%")%symbDensity %symb_velocity_x %symb_velocity_y %symb_velocity_z).str() ;
                std::string exprstr_diffusion_k = (boost::format("%1%+%2%/%3%:%1%:%2%:%3%")%symbDynViscosity %symbTurbulentDynViscosityBIS %symb_sigma_k ).str() ;//AIE
#if 0
                std::string exprstr_source_k = (boost::format("%1%*%2%^2-%3%*%4%:%1%:%2%:%3%:%4%")%symbTurbulentDynViscosity %symb_strain_rate_magnitude %symbDensity %symb_sol_epsilon).str();
#else
                std::string exprstr_source_k = (boost::format("%1%*%2%^2:%1%:%2%")%symbTurbulentDynViscosityBIS %symb_strain_rate_magnitude).str();
                std::string exprstr_reaction_k = (boost::format("%1%*%2%:%1%:%2%") %symbDensity %symb_gamma).str();
#endif

                ostr << "\""<<eqkeyword_k << "_beta\":\"" << exprstr_convection << "\","
                     << "\""<<eqkeyword_k << "_c\":\"" << exprstr_diffusion_k << "\","
                     << "\""<<eqkeyword_k << "_a\":\"" << exprstr_reaction_k << "\","
                    << "\""<<eqkeyword_k << "_f\":\"" << exprstr_source_k << "\"";
                if ( !this->isStationaryModel() )
                    ostr << ",\""<<eqkeyword_k << "_d\":\"" << (boost::format("%1%:%1%")%symbDensity).str() << "\"";

                std::string exprstr_diffusion_epsilon = (boost::format("%1%+%2%/%3%:%1%:%2%:%3%")%symbDynViscosity %symbTurbulentDynViscosityBIS %symb_sigma_epsilon ).str() ; //AIE
#if 0
                std::string exprstr_reaction_epsilon = (boost::format("-(%1%/%2%)*( %3%*%4%^2 ) + %5%*%6%*%7%/%2% :%1%:%2%:%3%:%4%:%5%:%6%:%7%") %symb_c_1epsilon %symb_sol_k %symbTurbulentDynViscosity %symb_strain_rate_magnitude %symb_c_2epsilon %symbDensity %symb_sol_epsilon ).str();
#else
                std::string exprstr_reaction_epsilon = (boost::format("%1%*%2%*%3%:%1%:%2%:%3%")%symb_c_2epsilon %symbDensity% symb_gamma ).str();


                // WARNING!!!!!
                std::string exprstr_source_epsilon = (boost::format("%1%*%2%*( %3%*%4%^2 ) :%1%:%2%:%3%:%4%")%symb_gamma %symb_c_1epsilon %symbTurbulentDynViscosityBIS %symb_strain_rate_magnitude ).str();
                //std::string exprstr_source_epsilon = (boost::format("%1%*%2%*%3%:%1%:%2%:%3%")%symb_c_1epsilon %symb_sol_k_positive %symb_strain_rate_magnitude ).str();
#endif
                ostr << ",";
                ostr << "\""<<eqkeyword_epsilon << "_beta\":\"" << exprstr_convection << "\","
                     << "\""<<eqkeyword_epsilon << "_c\":\"" << exprstr_diffusion_epsilon << "\","
                     << "\""<<eqkeyword_epsilon << "_a\":\"" << exprstr_reaction_epsilon << "\","
                     << "\""<<eqkeyword_epsilon << "_f\":\"" << exprstr_source_epsilon << "\"";
                if ( !this->isStationaryModel() )
                    ostr << ",\""<<eqkeyword_epsilon << "_d\":\"" << (boost::format("%1%:%1%")%symbDensity).str() << "\"";

#if 1
                //std::string mutExprStr = (boost::format("%1%*%2%*%3%^2/%4%:%1%:%2%:%3%:%4%")%symbDensity %symb_c_mu %symb_sol_k %symb_sol_epsilon ).str();
                //std::string mutExprStr = (boost::format("%1%*%2%*%3%^2/( %4%*%4%+1e-5 ) :%1%:%2%:%3%:%4%")%symbDensity %symb_c_mu %symb_sol_k %symb_sol_epsilon ).str();
                //std::string mutExprStr = (boost::format("max( 1e-8, %1%*%2%*%3%^2/( %4%*%4%+1e-5 ) ) :%1%:%2%:%3%:%4%")%symbDensity %symb_c_mu %symb_sol_k %symb_sol_epsilon ).str();
                //std::string mutExprStr = (boost::format("max( 1e-8, %1%*%2%*sqrt(%3%) ):%1%:%2%:%3%")%symbDensity %symb_c_mu %symb_sol_k_positive  ).str();
                //std::string mutExprStr = (boost::format("max( 1e-8, %1%*%2%*%3% ):%1%:%2%:%3%")%symbDensity %symb_l_star %symb_sol_k/*_positive*/ ).str(); // TODO define 1e-8 as mut min
                //std::string mutExprStr = (boost::format("%1%*%2%*%3%^2/max(1e-8,%4%) :%1%:%2%:%3%:%4%")%symbDensity %symb_c_mu %symb_sol_k %symb_sol_epsilon ).str();
                //std::string mutExprStr = (boost::format("%1%*%2%*%3%^2/(%4%*(%4%>1e-8)+ 1e-8*(1-(%4%>1e-8))) :%1%:%2%:%3%:%4%")%symbDensity %symb_c_mu %symb_sol_k %symb_sol_epsilon ).str();
                std::string mutExprStr = (boost::format("max( 1e-8, %1%*%2%*sqrt(max(1e-9,%3%) ) ):%1%:%2%:%3%")%symbDensity %symb_l_star %symb_sol_k_positive ).str(); // TODO define 1e-8 as mut min
#else
                std::string mutExprStr = (boost::format("max( 1e-8, %1%*%2%*sqrt(%3%) ):%1%:%2%:%3%")%symbDensity %symb_l_star %symb_sol_k_positive ).str(); // TODO define 1e-8 as mut min
#endif
                ModelExpression mutExpr;
                mutExpr.setExpr( mutExprStr, this->worldComm(), this->repository().expr() );
                this->materialsProperties()->addProperty( matProps, "turbulent-dynamic-viscosity", mutExpr, true );

                std::string tkeStr = (boost::format("%1%:%1%")%symb_sol_k).str();
                ModelExpression tkeExpr;
                tkeExpr.setExpr( tkeStr, this->worldComm(), this->repository().expr() );
                this->materialsProperties()->addProperty( matProps, "turbulent-kinetic-energy", tkeExpr, true );

                // initials coondition
                double mixing_length_limit = 0.035*0.0635/2.; //???
                std::ostringstream ostr_ic_k;
                ostr_ic_k << "\"markers\":[" << myCvrtSeqToStr( matProps.markers() ) <<  "],"
                          << "\"expr\":\"" << (boost::format(" (1*%1%/( %2%*0.1*%3% ) )^2 :%1%:%2%")%symbDynViscosity %symbDensity %mixing_length_limit ).str() << "\"";
                strInitialConditions[unknownName_k].push_back( ostr_ic_k.str() );
                std::ostringstream ostr_ic_epsilon;
                ostr_ic_epsilon << "\"markers\":[" << myCvrtSeqToStr( matProps.markers() ) <<  "],"
                                << "\"expr\":\"" << (boost::format("( (1*%1%/( %2%*0.1*%3% ) )^3 )*%4%/( %2%*0.1*%3% ) :%1%:%2%:%4%")%symbDynViscosity %symbDensity %mixing_length_limit %symb_c_mu ).str() << "\"";
                strInitialConditions[unknownName_epsilon].push_back( ostr_ic_epsilon.str() );

                // wall function
                std::string symb_y_plus_star = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_nomat, "y_plus_star", "_" ), 11.06 );
                std::string symb_u_tau = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "u_tau", "_" ),
                                                                        //(boost::format("%1%^(1/4)*sqrt(%2%):%1%:%2%") %symb_c_mu %symb_sol_k_previous_positive ).str(),
                                                                        (boost::format("max(%1%^(1/4)*sqrt(%2%), %3%/%4%):%1%:%2%:%3%:%4%") %symb_c_mu %symb_sol_k_previous_positive %symb_velocity_magintude %symb_y_plus_star ).str(),
                                                                        this->worldComm(), this->repository().expr() );

                std::string symb_u_tauBC = physicFluidData->addParameter( prefixvm( prefix_symbol_physic_mat, "u_tauBC", "_" ),
                                                                          (boost::format("%1%^(1/4)*sqrt(%2%):%1%:%2%") %symb_c_mu %symb_sol_k_previous_positive ).str(),
                                                                          this->worldComm(), this->repository().expr() );


            } // k-epsilon
            ostr << "}"; // end mymat1

        } // foreach mat
    }
    ostr << "}"; //end Materials


    ostr << ",";
    ostr << "\"InitialConditions\":{";
    bool writeFirstField = true;
    for ( auto const& [ nameField, icExprSubSections ] : strInitialConditions )
    {
        if ( !writeFirstField )
            ostr << ",";
        ostr << "\""<< nameField << "\":{";
        ostr << "\"Expression\":{";
        for ( int k=0;k<icExprSubSections.size();++k )
            ostr << "\"" << (boost::format("myic_%1%")%k ).str() << "\":{" << icExprSubSections[k] << "}";
        ostr << "}"; // end Expression
        ostr << "}"; // end solution
        writeFirstField = false;
    }
    ostr << "}"; //end InitialConditions


    ostr << ",";
    ostr << "\"BoundaryConditions\":{";
    if ( isSpalartAllmarasTurbulenceModel )
    {
        ostr << "\"" << eqkeyword << "\":{";
        ostr << "\"Dirichlet\":{";
        bool writeFirstBc=true;
        for ( auto const& [bcName,bcInlet] : M_turbulenceModelBoundaryConditions.inlet() )
        {
            if ( !writeFirstBc )
                ostr << ",";
            ostr << "\""<< bcName << "\":{"
                 << "\"markers\":[" << myCvrtSeqToStr( bcInlet.markers() ) << "],"
                 << "\"expr\":\"3*materials_mu/materials_rho:materials_mu:materials_rho\""
                 << "}";
            writeFirstBc = false;
        }
        for ( auto const& [bcName,bcWall] : M_turbulenceModelBoundaryConditions.wall() )
        {
            if ( !writeFirstBc )
                ostr << ",";
            ostr << "\""<< bcName << "\":{"
                 << "\"markers\":[" << myCvrtSeqToStr( bcWall.markers() ) << "],"
                 << "\"expr\":\"0\""
                 << "}";
            writeFirstBc = false;
        }
        ostr << "}"; // end Dirichlet
        ostr << "}"; // end eqkeyword
    }
    else // k-epsilon
    {
        std::string eqkeyword_k = std::get<0>( M_turbulenceModelType->pdes()[0] ).equationName();
        ostr << "\"" << eqkeyword_k << "\":{";
        ostr << "\"Dirichlet\":{";
        bool writeFirstBc=true;
        for ( auto const& [bcName,bcInlet] : M_turbulenceModelBoundaryConditions.inlet() )
        {
            double turbulenceIntensity = 0.05;
            //double umax = 15.6;
            double L = 0.0635/2.;
            double c_mu = 0.09;
            if ( !writeFirstBc )
                ostr << ",";
            ostr << "\""<< bcName << "\":{"
                 << "\"markers\":[" << myCvrtSeqToStr( bcInlet.markers() ) << "],"
                 << "\"expr\":\""<<(boost::format("(3/2)*(%1%*%2%)^2:%1%") %symb_velocity_magintude %turbulenceIntensity).str() << "\""
                 << "}";
            writeFirstBc = false;
        }
        ostr << "}"; // end Dirichlet
        ostr << "}"; // end eqkeyword_k
        ostr << ",";
        std::string eqkeyword_epsilon = std::get<0>( M_turbulenceModelType->pdes()[1] ).equationName();
        ostr << "\"" << eqkeyword_epsilon << "\":{";
        ostr << "\"Dirichlet\":{";
        writeFirstBc=true;
        for ( auto const& [bcName,bcInlet] : M_turbulenceModelBoundaryConditions.inlet() )
        {
            double turbulenceIntensity = 0.05;
            double umax = 15.6;
            double L = 0.0635/2.;
            double c_mu = 0.09;
            if ( !writeFirstBc )
                ostr << ",";
            ostr << "\""<< bcName << "\":{"
                 << "\"markers\":[" << myCvrtSeqToStr( bcInlet.markers() ) << "],"
                 << "\"expr\":\""<<(boost::format("%1%^(3/4)* ( ( (3/2)*(%2%*%3%)^2 )^(3/2) )/(0.035*%4% ) :%2%")%c_mu %symb_velocity_magintude %turbulenceIntensity %L).str() << "\""
                 << "}";
            writeFirstBc = false;
        }
        for ( auto const& [bcName,bcWall] : M_turbulenceModelBoundaryConditions.wall() )
        {
            // TODO get symbol from physics and mat
            std::string symb_u_tau = "physics_fluid_fluid_fluid_Omega_turbulence_k_epsilon_u_tau";
            std::string symbDynViscosity = "materials_Omega_mu";
            std::string symb_y_plus_star = "physics_fluid_fluid_fluid_turbulence_k_epsilon_y_plus_star";
            if ( !writeFirstBc )
                ostr << ",";
            ostr << "\""<< bcName << "\":{"
                 << "\"markers\":[" << myCvrtSeqToStr( bcWall.markers() ) << "],"
                 << "\"expr\":\""<<(boost::format("%1%^4/(0.41*%2%*%3%):%1%:%2%:%3%")%symb_u_tau %symb_y_plus_star %symbDynViscosity).str() << "\""
                 << "}";
        }

        ostr << "}"; // end Dirichlet
        ostr << "}"; // end eqkeyword_epsilon

    }
    ostr << "}"; //end BoundaryConditions

    ostr << "}";

    if ( this->worldComm().isMasterRank() )
        std::cout << "\n" << ostr.str()  << "\n\n";


    if ( !M_turbulenceModelType->hasModelProperties() )
    {
        pt::ptree pt;
        std::istringstream istr( ostr.str() );
        pt::read_json(istr, pt);
        if ( this->worldComm().isMasterRank() )
            pt::write_json( (fs::path( M_turbulenceModelType->repository().root() )/"turbulence.json").string(), pt);
        auto myModelProp = std::make_shared<ModelProperties>( pt, M_turbulenceModelType->repository().expr(), M_turbulenceModelType->worldCommPtr() );
        M_turbulenceModelType->setModelProperties( myModelProp );
    }
    M_turbulenceModelType->setMesh( this->mesh() );
    M_turbulenceModelType->setManageParameterValues( false );
    M_turbulenceModelType->init();

    M_turbulenceModelType->algebraicFactory()->setFunctionLinearAssembly( boost::bind( static_cast<void(self_type::*)(DataUpdateLinear&) const>(&self_type::updateLinear_Turbulence),
                                                                                       boost::ref( *this ), _1 ) );
    M_turbulenceModelType->algebraicFactory()->setFunctionLinearDofElimination( boost::bind( static_cast<void(self_type::*)(DataUpdateLinear&) const>(&self_type::updateLinearDofElimination_Turbulence),
                                                                                             boost::ref( *this ), _1 ) );
    M_turbulenceModelType->algebraicFactory()->setFunctionNewtonInitialGuess( boost::bind( static_cast<void(self_type::*)(DataNewtonInitialGuess&) const>(&self_type::updateNewtonInitialGuess_Turbulence),
                                                                                           boost::ref( *this ), _1 ) );
    M_turbulenceModelType->algebraicFactory()->setFunctionResidualAssembly( boost::bind( static_cast<void(self_type::*)(DataUpdateResidual&) const>(&self_type::updateResidual_Turbulence),
                                                                                         boost::ref( *this ), _1 ) );
    M_turbulenceModelType->algebraicFactory()->setFunctionJacobianAssembly( boost::bind( static_cast<void(self_type::*)(DataUpdateJacobian&) const>(&self_type::updateJacobian_Turbulence),
                                                                                         boost::ref( *this ), _1 ) );

    M_turbulenceModelType->printAndSaveInfo();

    // if ( this->worldComm().isMasterRank() )
    //     std::cout << "holalla turb \n "<< M_turbulenceModelType->modelFields().symbolsExpr().names() << std::endl;


}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::Body::setup( pt::ptree const& p, ModelMaterials const& mats, mesh_ptrtype mesh )
{
    M_mesh = mesh;
    std::set<std::string> matNames;
    if ( auto ptMatNames = p.get_child_optional("names") )
    {
        if ( ptMatNames->empty() ) // value case
            matNames.insert( ptMatNames->get_value<std::string>() );
        else // array case
        {
            for ( auto const& item : *ptMatNames )
            {
                CHECK( item.first.empty() ) << "should be an array, not a subtree";
                matNames.insert( item.second.template get_value<std::string>() );
            }
        }
    }

    ModelMarkers onlyMarkers;
    if ( auto ptmarkers = p.get_child_optional("markers") )
        onlyMarkers.setPTree(*ptmarkers/*, indexes*/);

    M_materialsProperties.reset( new materialsproperties_type( M_modelPhysics ) );
    M_materialsProperties->updateForUse( mats, matNames, onlyMarkers );
    M_materialsProperties->addMesh( M_mesh );

    // init displacement space
    auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
    auto M_rangeMeshElements = markedelements(this->mesh(), mom->markers( M_modelPhysics->physicsAvailableFromCurrentType() ) );
    M_spaceDisplacement = space_displacement_type::New(_mesh=M_mesh,_range=M_rangeMeshElements);
    M_fieldDisplacement = M_spaceDisplacement->elementPtr();

    this->updateForUse();

}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::Body::applyRemesh( mesh_ptrtype const& newMesh )
{
    if ( M_mesh )
    {
        mesh_ptrtype oldMesh = this->mesh();

        // material prop
        this->materialsProperties()->removeMesh( oldMesh );
        this->materialsProperties()->addMesh( newMesh );

        M_mesh = newMesh;

        // function space and fields
        space_displacement_ptrtype old_spaceDisplacement = M_spaceDisplacement;
        element_displacement_ptrtype old_fieldDisplacement = M_fieldDisplacement;
        element_displacement_ptrtype old_fieldElasticDisplacement = M_fieldElasticDisplacement;

        auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
        auto M_rangeMeshElements = markedelements(this->mesh(), mom->markers( M_modelPhysics->physicsAvailableFromCurrentType() ) );
        M_spaceDisplacement = space_displacement_type::New(_mesh=M_mesh,_range=M_rangeMeshElements);
        M_fieldDisplacement = M_spaceDisplacement->elementPtr();

        // createInterpolationOp
        auto opI_displacement = opInterpolation(_domainSpace=old_spaceDisplacement,
                                                _imageSpace=M_spaceDisplacement,
                                                _range=M_rangeMeshElements
                                                );

        auto matrixInterpolation_displacement = opI_displacement->matPtr();
        matrixInterpolation_displacement->multVector( *old_fieldDisplacement, *M_fieldDisplacement );


        if ( old_fieldElasticDisplacement )
        {
            M_fieldElasticDisplacement = M_spaceDisplacement->elementPtr();
            matrixInterpolation_displacement->multVector( *old_fieldElasticDisplacement, *M_fieldElasticDisplacement );
        }

        if ( M_fieldElasticVelocity )
        {
            space_velocity_ptrtype old_spaceElasticVelocity = M_spaceElasticVelocity;
            element_velocity_ptrtype old_fieldElasticVelocity = M_fieldElasticVelocity;
            M_fieldElasticVelocity.reset();
            this->initElasticVelocity();

            auto opI_elasticVelocity = opInterpolation(_domainSpace=old_spaceElasticVelocity,
                                                       _imageSpace=M_spaceElasticVelocity,
                                                       _range=M_rangeMeshElements );

            auto matrixInterpolation_elasticVelocity = opI_elasticVelocity->matPtr();
            matrixInterpolation_elasticVelocity->multVector( *old_fieldElasticVelocity, *M_fieldElasticVelocity );

        }

    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::Body::updateForUse()
{
    CHECK( M_materialsProperties ) << "no materialsProperties defined";

    auto mom = M_materialsProperties->materialsOnMesh(M_mesh);
    M_mass = 0;
    M_massCenter = eigen_vector_type<nRealDim>::Zero();
    for ( auto const& rangeData : mom->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& density = M_materialsProperties->density( matName );
        auto const& densityExpr = density.exprScalar();

        M_mass += integrate(_range=range,_expr=densityExpr).evaluate()(0,0);
        M_massCenter += integrate(_range=range,_expr=densityExpr*P()).evaluate();
    }
    M_massCenter /= M_mass;

    this->computeMomentOfInertia_bodyFrame( this->massCenterExpr(), this->rigidRotationMatrix(), M_momentOfInertia_bodyFrame );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::BodyBoundaryCondition( self_type const& fluidToolbox )
    :
    M_gravityForceEnabled( false )
{}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::setup( std::string const& bodyName, pt::ptree const& pt, self_type const& fluidToolbox )
{
    M_name = bodyName;
    if ( auto ptmarkers = pt.get_child_optional("markers") )
        M_markers.setPTree(*ptmarkers/*, indexes*/);
    else
        M_markers.insert( bodyName );

    if ( auto ptmaterials = pt.get_child_optional("materials") )
    {
        auto bodyPhysics = std::make_shared<ModelPhysics<nRealDim>>( "body", fluidToolbox );
        if ( bodyPhysics->physics().empty() )
            bodyPhysics->initPhysics( "body", ModelModels{}/*fluidToolbox.modelProperties().models()*/ );
        //if ( M_body->physics().empty() )
        //M_body->initPhysics( "body", ModelModels{}/*fluidToolbox.modelProperties().models()*/ );
        M_body.reset( new Body( bodyPhysics ) );
        M_body->setup( *ptmaterials, fluidToolbox.modelProperties().materials(), fluidToolbox.mesh() );
    }
    else
    {
        M_body.reset( new Body );

        ModelExpression massExpr, momentOfInertiaExpr, initialMassCenterExpr;
        massExpr.setExpr( "mass", pt, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
        if ( massExpr.template hasExpr<1,1>() )
            M_body->setMass( massExpr.template expr<1,1>().evaluate()(0,0) );
        momentOfInertiaExpr.setExpr( "moment-of-inertia", pt, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
        if constexpr ( nDim == 2 )
        {
            if ( momentOfInertiaExpr.template hasExpr<1,1>() )
                M_body->setMomentOfInertia_bodyFrame( momentOfInertiaExpr.template expr<1,1>().evaluate()(0,0) );
        }
        else
        {
            if ( momentOfInertiaExpr.template hasExpr<nDim,nDim>() )
                M_body->setMomentOfInertia_bodyFrame( momentOfInertiaExpr.template expr<nDim,nDim>().evaluate() );
        }
        initialMassCenterExpr.setExpr( "mass-center", pt, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
        if ( initialMassCenterExpr.template hasExpr<nRealDim,1>() )
        {
            auto initMassCenter = initialMassCenterExpr.template expr<nRealDim,1>();
            M_massCenterRef = initMassCenter.evaluate();
            M_body->setMassCenter( M_massCenterRef );
        }
    }

    M_translationalVelocityExpr.setExpr( "translational-velocity", pt, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
    M_angularVelocityExpr.setExpr( "angular-velocity", pt, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );

    if ( auto ptElasticVelocity = pt.get_child_optional("elastic-velocity") )
    {
        if ( ptElasticVelocity->empty() )
        {
            std::tuple< ModelExpression, std::set<std::string>> dataExpr;
            std::get<0>( dataExpr ).setExpr( "elastic-velocity", pt, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
            if ( std::get<0>( dataExpr ).template hasExpr<nDim,1>() )
                M_elasticVelocityExprBC.emplace( "", dataExpr );
        }
        else
        {
            for ( auto const& item : *ptElasticVelocity )
            {
                std::string bcElasticVelocityName = item.first;
                std::tuple< ModelExpression, std::set<std::string>> dataExpr;
                std::get<0>( dataExpr ).setExpr( "expr", item.second, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
                if ( !std::get<0>( dataExpr ).template hasExpr<nDim,1>() )
                    continue;
                ModelMarkers bcElasticVelocityMarkers;
                if ( auto ptmarkers = item.second.get_child_optional("markers") )
                    bcElasticVelocityMarkers.setPTree(*ptmarkers/*, indexes*/);
                std::get<1>( dataExpr ) = bcElasticVelocityMarkers;
                M_elasticVelocityExprBC.emplace( bcElasticVelocityName, dataExpr );
            }
        }
    }

    if ( auto ptElasticDisplacement = pt.get_child_optional("elastic-displacement") )
    {
        if ( ptElasticDisplacement->empty() )
        {
            std::tuple< ModelExpression, std::set<std::string>> dataExpr;
            std::get<0>( dataExpr ).setExpr( "elastic-displacement", pt, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
            if ( std::get<0>( dataExpr ).template hasExpr<nDim,1>() )
                M_elasticDisplacementExprBC.emplace( "", dataExpr );
        }
        else
        {
            for ( auto const& item : *ptElasticDisplacement )
            {
                std::string bcElasticDisplacementName = item.first;
                std::tuple< ModelExpression, std::set<std::string>> dataExpr;
                std::get<0>( dataExpr ).setExpr( "expr", item.second, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
                if ( !std::get<0>( dataExpr ).template hasExpr<nDim,1>() )
                    continue;
                ModelMarkers bcElasticDisplacementMarkers;
                if ( auto ptmarkers = item.second.get_child_optional("markers") )
                    bcElasticDisplacementMarkers.setPTree(*ptmarkers/*, indexes*/);
                std::get<1>( dataExpr ) = bcElasticDisplacementMarkers;
                M_elasticDisplacementExprBC.emplace( bcElasticDisplacementName, dataExpr );
            }
        }
    }


    if ( auto ptArticulation = pt.get_child_optional("articulation") )
    {
        if ( auto ptBodyName = ptArticulation->get_optional<std::string>( "body" ) )
        {
            M_articulationTranslationalVelocityExpr[*ptBodyName].setExpr( "translational-velocity", *ptArticulation, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
            if ( !M_articulationTranslationalVelocityExpr[*ptBodyName].template hasExpr<1,1>() )
                CHECK( false ) << "required a scalar expr";
        }
       else
            CHECK( false ) << "require body";
        //M_articulationBodiesUsed->begin()->second M_articulationTranslationalVelocityExpr.setExpr( "translational-velocity", *ptArticulation, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );

    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::init( self_type const& fluidToolbox )
{
    auto rangeBodyBoundary = markedfaces(fluidToolbox.mesh(), std::set<std::string>(M_markers) );
    M_rangeMarkedFacesOnFluid = rangeBodyBoundary;
    M_mesh = createSubmesh(_mesh=/*fluidToolbox.mesh()*/fluidToolbox.functionSpaceVelocity()->template meshSupport<0>(),_range=rangeBodyBoundary,_view=true );

    this->initFunctionSpaces();

    // maybe enable gravity (TODO : select only physic with gravity)
    for ( auto const& [physicName,physicData] : fluidToolbox.physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
        if ( physicFluidData->gravityForceEnabled() )
        {
            M_gravityForceEnabled = true;
            break;
        }
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::initFunctionSpaces()
{
    M_spaceTranslationalVelocity = space_trace_p0c_vectorial_type::New( _mesh=M_mesh );
    M_fieldTranslationalVelocity = M_spaceTranslationalVelocity->elementPtr();

    if ( !this->isInNBodyArticulated() )
    {
        if constexpr ( nDim == 2 )
            M_spaceAngularVelocity = space_trace_angular_velocity_type::New( _mesh=M_mesh );
        else
            M_spaceAngularVelocity = M_spaceTranslationalVelocity;
        M_fieldAngularVelocity = M_spaceAngularVelocity->elementPtr();
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::initElasticBehavior()
{
    this->body().initElasticDisplacement();
    this->body().initElasticVelocity();
    if ( !M_fieldElasticVelocity )
    {
        M_spaceElasticVelocity = space_trace_velocity_type::New( _mesh=M_mesh );
        M_fieldElasticVelocity = M_spaceElasticVelocity->elementPtr();



        if ( M_mesh->worldComm().size() == 1 )
        {
            // TODO : NOT WORK IN // with some special partitioning
            // idea to fix without change opInterlation : introduce an intermediate functionspace
            // defined on mesh with elements connected to the interface only. And then create 2 opInterp
            auto opInterpolationElasticVelocity = opInterpolation(_domainSpace=this->body().fieldElasticVelocity().functionSpace(),
                                                                  _imageSpace=M_spaceElasticVelocity,
                                                                  _range=elements(M_mesh),
                                                                  _type=InterpolationConforme()/*,
                                                                                                _backend=M_fluidModel->backend()*/ );
            M_matrixInterpolationElasticVelocity = opInterpolationElasticVelocity->matPtr();

        }
        else
        {
#if 0
            typename trace_mesh_type::smd_ptrtype smd;
            smd = M_mesh->subMeshData();
            M_mesh->setSubMeshData( typename trace_mesh_type::smd_ptrtype{} );
#endif
            M_spaceElasticVelocityTouchingBodyInterface = space_velocity_type::New( _mesh=M_body->mesh(),_range=elements( M_body->mesh(), this->rangeMarkedFacesOnFluid()/*, ElementsType::MESH_FACES*/ ) );
            //M_spaceElasticVelocityTouchingBodyInterface = space_velocity_type::New( _mesh=M_body->mesh() );
            M_fieldElasticVelocityTouchingBodyInterface = M_spaceElasticVelocityTouchingBodyInterface->elementPtr();
            auto opInterpolationElasticVelocity_tmp1 = opInterpolation(_domainSpace=this->body().fieldElasticVelocity().functionSpace(),
                                                                       _imageSpace=M_spaceElasticVelocityTouchingBodyInterface,
                                                                       //_range= this->rangeMarkedFacesOnFluid(),
                                                                       //_range=elements(support(this->body().fieldElasticVelocity().functionSpace())),
                                                                       _range=intersect(elements(support(this->body().fieldElasticVelocity().functionSpace())),elements(support(M_spaceElasticVelocityTouchingBodyInterface))),
                                                                       _type=InterpolationConforme() );

            auto opInterpolationElasticVelocity_tmp2 = opInterpolation(_domainSpace=M_spaceElasticVelocityTouchingBodyInterface,
                                                                       _imageSpace=M_spaceElasticVelocity,
                                                                       _range=elements(M_mesh),
                                                                       _type=InterpolationConforme() );

            M_matrixInterpolationElasticVelocity_tmp1 = opInterpolationElasticVelocity_tmp1->matPtr();
            M_matrixInterpolationElasticVelocity_tmp2 = opInterpolationElasticVelocity_tmp2->matPtr();
#if 0
            if ( smd )
                M_mesh->setSubMeshData( smd );
#endif
        }


    }

}

namespace utility_constant_functionspace
{

std::shared_ptr<datamap_t<>>
aggregateParallelSupport( std::vector<std::shared_ptr<datamap_t<>>> const& datamaps, worldcomm_ptr_t const& thecomm )
{
    // create a new DataMap that allow to couple constant spaces (map1 and map2) which are not defined on same processor.
    // active dofs are assigned on the proc which contains the active dofs of map1.
    // ghost dofs are defined on all proc that contains active/ghost in map1 and map2.

    CHECK( !datamaps.empty() ) << "no data map";

    // TODO check some properties nDof ..

    auto const& mapRef = datamaps.front();
    auto mapNew = std::make_shared<datamap_t<>>( thecomm );
    // define number of dof (with and without ghost)
    mapNew->setNDof( mapRef->nDof() );
    rank_type themasterRank = 0;
    for (rank_type p=0;p<thecomm->localSize();++p)
    {
        //auto mapUsed =
        auto itMapUsed = std::find_if(datamaps.begin(), datamaps.end(), [&p](auto const& x){return x->nLocalDofWithGhost(p) > 0 ;});
        if ( itMapUsed != datamaps.end() )
        {
            auto const& mapUsed = *itMapUsed;
            mapNew->setNLocalDofWithGhost( p, mapUsed->nLocalDofWithGhost(p) );
            mapNew->setLastDof( p, mapUsed->nLocalDofWithGhost(p) - 1 );
        }
        if ( mapRef->nLocalDofWithoutGhost(p) > 0 )
        {
            themasterRank = p;
            mapNew->setNLocalDofWithoutGhost( p, mapRef->nLocalDofWithoutGhost(p) );
            mapNew->setLastDofGlobalCluster( p, mapRef->nLocalDofWithoutGhost(p) - 1 );
        }
    }
    // define mapping GlobalProcessToGlobalCluster
    auto itMapUsed = std::find_if(datamaps.begin(), datamaps.end(), [](auto const& x){return x->nLocalDofWithGhost() > 0 ;});
    if ( itMapUsed != datamaps.end() )
    {
        auto const& mapUsed = *itMapUsed;
        mapNew->setMapGlobalProcessToGlobalCluster( mapUsed->mapGlobalProcessToGlobalCluster() );
    }

    // add all partition as neighbor (if has localdof)
    if ( mapNew->nLocalDofWithGhost() > 0 )
    {
        for ( rank_type proc=0; proc<mapNew->worldComm().localSize(); ++proc )
            if ( proc!=mapNew->worldComm().localRank() && mapNew->nLocalDofWithGhost(proc) > 0 )
                mapNew->addNeighborSubdomain( proc );
    }
    // update activeDofSharedOnCluster
    if ( themasterRank == mapNew->worldComm().localRank() )
    {
        for ( size_type k = 0; k < mapNew->nLocalDofWithGhost() ; ++k )
        {
            for ( rank_type proc=0; proc<mapNew->worldComm().localSize(); ++proc )
            {
                if ( proc == themasterRank ) continue;
                if ( mapNew->nLocalDofWithGhost(proc) == 0 ) continue;
                mapNew->setActiveDofSharedOnCluster(k, { proc });
            }
        }
    }
    mapNew->initNumberOfDofIdToContainerId( 1 );
    mapNew->initDofIdToContainerIdIdentity( 0,mapNew->nLocalDofWithGhost() );
    mapNew->buildIndexSplit();
    // mapNew->buildIndexSplitWithComponents( nRealDim );

    return mapNew;
}

std::shared_ptr<GraphCSR>
aggregateGraph( std::vector<std::shared_ptr<GraphCSR>> const& graphs, std::shared_ptr<datamap_t<>> mapCol )
{
    auto sparsity_graph = std::make_shared<GraphCSR>( graphs.front()->mapRowPtr(), mapCol );
    for ( auto const& g : graphs )
    {
        for ( auto it = g->begin(), en=g->end(); it!=en ;++it )
        {
            auto/*row_type*/ const& irow = it->second;
            size_type ig = it->first;
            auto& row = sparsity_graph->row( ig );
            row.get<0>() = irow.get<0>(); // rank
            row.get<1>() = irow.get<1>(); // rank
            row.get<2>().insert( irow.get<2>().begin(), irow.get<2>().end() );
        }
    }
    sparsity_graph->close();
    return sparsity_graph;
}

} // namespace utility_constant_functionspace

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::updateForUse( self_type const& fluidToolbox )
{
    if ( !M_mesh )
        this->init( fluidToolbox );

    if ( M_body->hasMaterialsProperties() )
    {
        M_body->updateForUse();
    }
    else if ( fluidToolbox.isMoveDomain() )
    {
        auto disp = mean(_range=M_rangeMarkedFacesOnFluid,_expr=idv(fluidToolbox.meshALE()->displacement()) );
        //M_massCenter = M_massCenterRef + disp;
        M_body->setMassCenter( M_massCenterRef + disp );
    }

    //if ( fluidToolbox.worldComm().isMasterRank() )
    //    std::cout << "M_massCenter=\n " << M_body->massCenter() << std::endl;

    auto XhV = fluidToolbox.functionSpaceVelocity();

    XhV->rebuildDofPoints(); // TODO REMOVE this line, interpolation operator must not use dofPoint!!!

    // matrix interpolation of translational velocity
    if ( !M_matrixPTilde_translational )
    {
        auto opI_partTranslationalVelocity = opInterpolation( _domainSpace=this->spaceTranslationalVelocity() ,_imageSpace=XhV,_range=M_rangeMarkedFacesOnFluid );
        M_matrixPTilde_translational = opI_partTranslationalVelocity->matPtr();
    }

    if ( !this->isInNBodyArticulated() )
        this->updateMatrixPTilde_angular( fluidToolbox, M_matrixPTilde_angular );



    if ( M_gravityForceEnabled )
    {
        // TODO !!!
        auto const& matNames = fluidToolbox.materialsProperties()->physicToMaterials( fluidToolbox.physicsAvailableFromCurrentType() );
        CHECK( matNames.size() == 1 ) << "support only one";
        std::string matName = *matNames.begin();
        auto const& densityExpr = fluidToolbox.materialsProperties()->density( matName ).template expr<1,1>();
        double rho = densityExpr.evaluate()(0,0);
#if 0
        CHECK( fluidToolbox.materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
        std::string matName = fluidToolbox.materialProperties()->rangeMeshElementsByMaterial().begin()->first;
        double rho = fluidToolbox.materialProperties()->cstDensity( matName );
        //auto const& rho = fluidToolbox.materialProperties()->fieldRho();
#endif
        double massBody = M_body->mass();
        double massOfFluid = M_body->evaluateMassFromDensity( cst( rho ) );
        // if ( Environment::isMasterRank() )
        // {
        //     std::cout << "massBody = " << massBody << std::endl;
        //     std::cout << "massOfFluid = " << massOfFluid << std::endl;
        // }

        for ( auto const& [physicName,physicData] : fluidToolbox.physicsFromCurrentType() )
        {
            auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
            if ( physicFluidData->gravityForceEnabled() )
            {
                M_gravityForceWithMass = (massBody-massOfFluid)*physicFluidData->gravityForceExpr().evaluate();
            }
        }
    }

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::applyRemesh( self_type const& fluidToolbox/*, mesh_ptrtype const& newMesh*/ ,RemeshInterpolation & remeshInterp )
{
    if ( M_body )
        M_body->applyRemesh( fluidToolbox.mesh() );

    auto rangeBodyBoundary = markedfaces(fluidToolbox.mesh(), std::set<std::string>(M_markers) );
    M_rangeMarkedFacesOnFluid = rangeBodyBoundary;
    M_mesh = createSubmesh(_mesh=/*fluidToolbox.mesh()*/fluidToolbox.functionSpaceVelocity()->template meshSupport<0>(),_range=rangeBodyBoundary,_view=true );

#if 0
    typename Body::translational_velocity_type translationalVelocity = Body::translational_velocity_type::Zero();
    typename Body::angular_velocity_type angularVelocity = Body::angular_velocity_type::Zero();
    if ( M_fieldTranslationalVelocity )
        translationalVelocity = idv(M_fieldTranslationalVelocity).evaluate();
    if ( M_fieldAngularVelocity )
        angularVelocity = idv(M_fieldAngularVelocity).evaluate();
#endif
    space_trace_p0c_vectorial_ptrtype old_spaceTranslationalVelocity = M_spaceTranslationalVelocity;
    space_trace_angular_velocity_ptrtype old_spaceAngularVelocity = M_spaceAngularVelocity;
    element_trace_p0c_vectorial_ptrtype old_fieldTranslationalVelocity = M_fieldTranslationalVelocity;
    element_trace_angular_velocity_ptrtype old_fieldAngularVelocity = M_fieldAngularVelocity;
    this->initFunctionSpaces();
#if 0
    if ( M_fieldTranslationalVelocity )
    {
        for (int i=0;i<M_fieldTranslationalVelocity->functionSpace()->nLocalDofWithGhost();++i)
            M_fieldTranslationalVelocity->set( i, translationalVelocity(i) );
    }

    if ( M_fieldAngularVelocity )
    {
        for (int i=0;i<M_fieldAngularVelocity->functionSpace()->nLocalDofWithGhost();++i)
            M_fieldAngularVelocity->set( i, angularVelocity(i) );
    }
#endif

    if ( old_spaceTranslationalVelocity )
    {
        auto matrixInterpolation_translationalVelocity = fluidToolbox.backend()->newIdentityMatrix( old_spaceTranslationalVelocity->dof(), M_spaceTranslationalVelocity->dof() );
        remeshInterp.setMatrixInterpolation( old_spaceTranslationalVelocity, M_spaceTranslationalVelocity, matrixInterpolation_translationalVelocity );
        remeshInterp.registeringBlockIndex( fluidToolbox.keyword(), fluidToolbox.startSubBlockSpaceIndex( "body-bc."+this->name()+".translational-velocity" ), old_spaceTranslationalVelocity, M_spaceTranslationalVelocity );
        remeshInterp.interpolate( old_fieldTranslationalVelocity, M_fieldTranslationalVelocity );
        if ( M_bdfTranslationalVelocity )
            M_bdfTranslationalVelocity->applyRemesh( M_spaceTranslationalVelocity, matrixInterpolation_translationalVelocity );
    }

    if ( old_spaceAngularVelocity )
    {
        auto matrixInterpolation_angularVelocity = fluidToolbox.backend()->newIdentityMatrix( old_spaceAngularVelocity->dof(), M_spaceAngularVelocity->dof() );
        remeshInterp.setMatrixInterpolation( old_spaceAngularVelocity, M_spaceAngularVelocity, matrixInterpolation_angularVelocity );
        //if ( !this->isInNBodyArticulated() ) // maybe
        remeshInterp.registeringBlockIndex( fluidToolbox.keyword(), fluidToolbox.startSubBlockSpaceIndex( "body-bc."+this->name()+".angular-velocity" ), old_spaceAngularVelocity, M_spaceAngularVelocity );
        remeshInterp.interpolate( old_fieldAngularVelocity, M_fieldAngularVelocity );
        if ( M_bdfAngularVelocity )
            M_bdfAngularVelocity->applyRemesh( M_spaceAngularVelocity, matrixInterpolation_angularVelocity );
    }

    if ( M_fieldElasticVelocity )
    {
        M_fieldElasticVelocity.reset();
        this->initElasticBehavior();
        // TODO : check if we need to init ElasticVelocity from ElasticVelocity of the body
    }

    auto XhV = fluidToolbox.functionSpaceVelocity();
    XhV->rebuildDofPoints(); // TODO REMOVE this line, interpolation operator must not use dofPoint!!!

    // matrix interpolation of translational velocity
    if ( M_matrixPTilde_translational )
    {
        auto opI_partTranslationalVelocity = opInterpolation( _domainSpace=this->spaceTranslationalVelocity() ,_imageSpace=XhV,_range=M_rangeMarkedFacesOnFluid );
        M_matrixPTilde_translational = opI_partTranslationalVelocity->matPtr();
    }

    if ( M_matrixPTilde_angular )
    {
        M_matrixPTilde_angular.reset();
        this->updateMatrixPTilde_angular( fluidToolbox, M_matrixPTilde_angular );
    }

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::updateInformationObject( nl::json & p ) const
{
    p["mesh"] = M_mesh->journalSection().to_string();
    p["face markers on fluid mesh"] = M_markers;
    p["gravity force enabled"] = M_gravityForceEnabled;
    if ( this->hasTranslationalVelocityExpr() )
        p["translational velocity expr"] = M_translationalVelocityExpr.exprToString();
    if ( this->hasAngularVelocityExpr() )
        p["angular velocity expr"] =  M_angularVelocityExpr.exprToString();
    p["has elastic velocity"] = this->hasElasticVelocity();
    if ( this->hasElasticVelocityFromExpr() )
    {
        // TODO
    }

    nl::json subPt;
    subPt["Translational Velocity"] = this->spaceTranslationalVelocity()->journalSection().to_string();
    subPt["Angular Velocity"] = this->spaceAngularVelocity()->journalSection().to_string();
    p["Function Spaces"] = subPt;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp, std::map<std::string,uint16_type> & jsonPtrFunctionSpacesToLevel ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );

    Feel::Table tabInfoOthers;
    TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoOthers, jsonInfo, tabInfoProp );

    if ( jsonInfo.contains( "face markers on fluid mesh" ) )
        tabInfoOthers.add_row({ "face markers on fluid mesh", TabulateInformationTools::FromJSON::createTableFromArray( jsonInfo.at( "face markers on fluid mesh" ), true ) });

    tabInfoOthers.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    tabInfo->add( "", TabulateInformations::New( tabInfoOthers,tabInfoProp ) );

    Feel::Table tabInfoAdvanced;
    if ( jsonInfo.contains("Function Spaces") )
    {
        auto const& jsonInfoFunctionSpaces = jsonInfo.at("Function Spaces");
        Feel::Table tabInfoFunctionSpaces;
        for ( std::string const& spaceName : std::vector<std::string>({"Translational Velocity","Angular Velocity"}) )
        {
            std::string jsonPtrStr = jsonInfoFunctionSpaces.at( spaceName ).template get<std::string>();
            nl::json::json_pointer jsonPointerSpace( jsonPtrStr );
            if ( JournalManager::journalData().contains( jsonPointerSpace ) )
            {
                tabInfoFunctionSpaces.add_row( {spaceName, jsonPtrStr } );

                uint16_type levelFunctionSpace = tabInfoProp.verboseLevel()+1;
                auto itFindJsonPtr = jsonPtrFunctionSpacesToLevel.find( jsonPtrStr );
                if ( itFindJsonPtr != jsonPtrFunctionSpacesToLevel.end() )
                    levelFunctionSpace = std::min( itFindJsonPtr->second, levelFunctionSpace );
                jsonPtrFunctionSpacesToLevel[jsonPtrStr] = levelFunctionSpace;
            }
        }
        tabInfoAdvanced.add_row({"Function Spaces",tabInfoFunctionSpaces});
    };
    tabInfo->add( "", TabulateInformations::New( tabInfoAdvanced,tabInfoProp.newByIncreasingVerboseLevel() ) );

    return tabInfo;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::updateMatrixPTilde_angular( self_type const& fluidToolbox, sparse_matrix_ptrtype & mat, size_type startBlockIndexVelocity, size_type startBlockIndexAngularVelocity ) const
{
    auto XhV = fluidToolbox.functionSpaceVelocity();
    auto const& w = *this->fieldAngularVelocityPtr();
    // matrix interpolation with angular velocity expr (depends on mesh position and mass center -> rebuild at each call of updateForUse)
    auto massCenter = this->massCenterExpr();

    bool buildNewMatrix = mat? false : true;

    OperatorInterpolationMatrixSetup matSetup( mat, buildNewMatrix? Feel::DIFFERENT_NONZERO_PATTERN : Feel::SAME_NONZERO_PATTERN,
                                               startBlockIndexVelocity, startBlockIndexAngularVelocity );
    if constexpr (nDim == 2 )
    {
        auto opI_AngularVelocity = opInterpolation( _domainSpace=this->spaceAngularVelocity(),_imageSpace=XhV,_range=M_rangeMarkedFacesOnFluid,
                                                    _type= makeExprInterpolation( id(w)*vec(-Py()+massCenter(1,0),Px()-massCenter(0,0) ), nonconforming_t() ),
                                                    _matrix=matSetup );
        if ( buildNewMatrix )
            mat = opI_AngularVelocity->matPtr();
        // M_matrixPTilde_angular = opI_AngularVelocity->matPtr();
    }
    else
    {
        auto r = vec(Px()-massCenter(0,0),Py()-massCenter(1,0),Pz()-massCenter(2,0) );
        auto opI_AngularVelocity = opInterpolation( _domainSpace=this->spaceAngularVelocity(),
                                                    _imageSpace=XhV,
                                                    _range=M_rangeMarkedFacesOnFluid,
                                                    _type= makeExprInterpolation( cross(id(w),r), nonconforming_t() ),
                                                    _matrix=matSetup );
        if ( buildNewMatrix )
            mat = opI_AngularVelocity->matPtr();
        // M_matrixPTilde_angular = opI_AngularVelocity->matPtr();
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::NBodyArticulated::updateMatrixPTilde_angular( self_type const& fluidToolbox, sparse_matrix_ptrtype & mat, size_type startBlockIndexVelocity, size_type startBlockIndexAngularVelocity ) const
{
    auto XhV = fluidToolbox.functionSpaceVelocity();
    auto const& w = *this->fieldAngularVelocityPtr();
    // matrix interpolation with angular velocity expr (depends on mesh position and mass center -> rebuild at each call of updateForUse)
    auto massCenter = this->massCenterExpr();

    bool buildNewMatrix = mat? false : true;

    OperatorInterpolationMatrixSetup matSetup( mat, buildNewMatrix? Feel::DIFFERENT_NONZERO_PATTERN : Feel::SAME_NONZERO_PATTERN,
                                               startBlockIndexVelocity, startBlockIndexAngularVelocity );
    if constexpr (nDim == 2 )
    {
        auto opI_AngularVelocity = opInterpolation( _domainSpace=this->spaceAngularVelocity() ,_imageSpace=XhV,_range=M_rangeMarkedFacesOnFluid,
                                                    _type= makeExprInterpolation( id(w)*vec(-Py()+massCenter(1,0),Px()-massCenter(0,0) ), nonconforming_t() ),
                                                    _matrix=matSetup );
        if ( buildNewMatrix )
            mat = opI_AngularVelocity->matPtr();
        // M_matrixPTilde_angular = opI_AngularVelocity->matPtr();
    }
    else
    {
        auto r = vec(Px()-massCenter(0,0),Py()-massCenter(1,0),Pz()-massCenter(2,0) );
        auto opI_AngularVelocity = opInterpolation( _domainSpace=this->spaceAngularVelocity(),
                                                    _imageSpace=XhV,
                                                    _range=M_rangeMarkedFacesOnFluid,
                                                    _type= makeExprInterpolation( cross(id(w),r), nonconforming_t() ),
                                                    _matrix=matSetup );
        if ( buildNewMatrix )
            mat = opI_AngularVelocity->matPtr();
        // M_matrixPTilde_angular = opI_AngularVelocity->matPtr();
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::NBodyArticulated::updateMatrixPTilde_angular( self_type const& fluidToolbox )
{
    this->updateMatrixPTilde_angular( fluidToolbox, M_matrixPTilde_angular );
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::string
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyArticulation::name() const { return this->body1().name() + "_" + this->body2().name(); }

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
bool
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyArticulation::has( BodyBoundaryCondition const& bbc ) const { return (bbc.name() == this->body1().name()) || (bbc.name() == this->body2().name()); }

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyArticulation::initLagrangeMultiplier( self_type const& fluidToolbox )
{
    M_dataMapLagrangeMultiplierTranslationalVelocity =
        utility_constant_functionspace::aggregateParallelSupport( { this->body1().spaceTranslationalVelocity()->mapPtr(), this->body2().spaceTranslationalVelocity()->mapPtr() },
                                                                  fluidToolbox.worldCommPtr() );
    M_vectorLagrangeMultiplierTranslationalVelocity = fluidToolbox.backend()->newVector( M_dataMapLagrangeMultiplierTranslationalVelocity );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyArticulation::applyRemesh( self_type const& fluidToolbox, RemeshInterpolation & remeshInterp )
{
    if ( M_vectorLagrangeMultiplierTranslationalVelocity )
    {
        vector_ptrtype old_vectorLagrangeMultiplierTranslationalVelocity = M_vectorLagrangeMultiplierTranslationalVelocity;
        this->initLagrangeMultiplier( fluidToolbox );

        auto matInterp = fluidToolbox.backend()->newIdentityMatrix( old_vectorLagrangeMultiplierTranslationalVelocity->mapPtr(), M_vectorLagrangeMultiplierTranslationalVelocity->mapPtr() );
        std::string subBlockSpaceName = "body-bc.articulation-lm."+this->name()+".translational-velocity";
        std::string interpName = prefixvm(fluidToolbox.keyword(),subBlockSpaceName);
        remeshInterp.setMatrixInterpolation( interpName, old_vectorLagrangeMultiplierTranslationalVelocity->mapPtr(), M_vectorLagrangeMultiplierTranslationalVelocity->mapPtr(), matInterp );
        remeshInterp.registeringBlockIndex( fluidToolbox.keyword(), fluidToolbox.startSubBlockSpaceIndex(subBlockSpaceName), interpName );
        remeshInterp.interpolate( interpName, old_vectorLagrangeMultiplierTranslationalVelocity, M_vectorLagrangeMultiplierTranslationalVelocity );
    }

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
eigen_vector_type<FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::nRealDim>
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyArticulation::unitDirBetweenMassCenters() const
{
    auto mc1 = this->body1().body().massCenter();
    auto mc2 = this->body2().body().massCenter();
    return this->unitDirBetweenMassCenters( mc1, mc2 );
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::string
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::NBodyArticulated::name() const { return "articulation_"+this->masterBodyBC().name(); }


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::NBodyArticulated::init( self_type const& fluidToolbox )
{
    std::map<std::string,std::tuple<int,BodyBoundaryCondition const*>> count;
    for ( auto const& ba : M_articulations )
    {
        for ( BodyBoundaryCondition const* bbcPtr : { &ba.body1(), &ba.body2() } )
        {
            auto const& bbc = *bbcPtr;
            if ( count.find( bbc.name() ) == count.end() )
                count[bbc.name()] = std::make_tuple(1,bbcPtr);
            else
                ++std::get<0>( count[bbc.name()] );
        }
    }
    auto itMax = std::max_element(count.begin(), count.end(), [](auto const& e1,auto const& e2) { return std::get<0>(e1.second) < std::get<0>(e2.second); } );
    CHECK ( itMax != count.end() ) << "something wrong";
    M_masterBodyBC = std::get<1>( itMax->second );
    CHECK( M_masterBodyBC->name() == itMax->first ) << "something wrong";

    // std::cout << "M_pmatrixMasterBodyName = "<< M_pmatrixMasterBodyName << std::endl;
    // std::cout << "this->bodyList().size()="<< this->bodyList().size() << std::endl;
    // std::cout << "this->bodyList(false).size()="<< this->bodyList(false).size() << std::endl;


    std::set<std::string> M_markers;
    for ( auto bbc : this->bodyList() )
        M_markers.insert( bbc->markers().begin(), bbc->markers().end() );
    M_rangeMarkedFacesOnFluid = markedfaces(fluidToolbox.mesh(), M_markers );


    if ( this->articulationMethod() == "lm" )
    {
        for ( auto & ba : M_articulations )
            ba.initLagrangeMultiplier( fluidToolbox );
    }
    else
    {
        // create datamap shared on all processess where a body bc is defined (active dofs use the master body)
        std::vector<std::shared_ptr<datamap_t<>>> datamaps = { this->masterBodyBC().spaceTranslationalVelocity()->mapPtr() };
        for ( auto const& bbcPtr : this->bodyList(false) )
            datamaps.push_back( bbcPtr->spaceTranslationalVelocity()->mapPtr() );
        M_dataMapPMatrixTranslationalVelocity = utility_constant_functionspace::aggregateParallelSupport( datamaps, fluidToolbox.worldCommPtr() );
    }


    auto M_mesh = createSubmesh(_mesh=/*fluidToolbox.mesh()*/fluidToolbox.functionSpaceVelocity()->template meshSupport<0>(),_range=M_rangeMarkedFacesOnFluid,_view=true );
    M_spaceAngularVelocity = space_trace_angular_velocity_type::New( _mesh=M_mesh );
    M_fieldAngularVelocity = M_spaceAngularVelocity->elementPtr();
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::NBodyArticulated::applyRemesh( self_type const& fluidToolbox, RemeshInterpolation & remeshInterp )
{
    std::set<std::string> M_markers;
    for ( auto bbc : this->bodyList() )
        M_markers.insert( bbc->markers().begin(), bbc->markers().end() );
    M_rangeMarkedFacesOnFluid = markedfaces(fluidToolbox.mesh(), M_markers );

    for ( auto & ba : M_articulations )
        ba.applyRemesh( fluidToolbox, remeshInterp );

    if ( M_dataMapPMatrixTranslationalVelocity )
    {
        // create datamap shared on all processess where a body bc is defined (active dofs use the master body)
        std::vector<std::shared_ptr<datamap_t<>>> datamaps = { this->masterBodyBC().spaceTranslationalVelocity()->mapPtr() };
        for ( auto const& bbcPtr : this->bodyList(false) )
            datamaps.push_back( bbcPtr->spaceTranslationalVelocity()->mapPtr() );
        M_dataMapPMatrixTranslationalVelocity = utility_constant_functionspace::aggregateParallelSupport( datamaps, fluidToolbox.worldCommPtr() );
    }


    auto M_mesh = createSubmesh(_mesh=/*fluidToolbox.mesh()*/fluidToolbox.functionSpaceVelocity()->template meshSupport<0>(),_range=M_rangeMarkedFacesOnFluid,_view=true );

    space_trace_angular_velocity_ptrtype old_spaceAngularVelocity = M_spaceAngularVelocity;
    element_trace_angular_velocity_ptrtype old_fieldAngularVelocity = M_fieldAngularVelocity;
    M_spaceAngularVelocity = space_trace_angular_velocity_type::New( _mesh=M_mesh );
    M_fieldAngularVelocity = M_spaceAngularVelocity->elementPtr();

    auto matrixInterpolation_angularVelocity = fluidToolbox.backend()->newIdentityMatrix( old_spaceAngularVelocity->dof(), M_spaceAngularVelocity->dof() );
    remeshInterp.setMatrixInterpolation( old_spaceAngularVelocity, M_spaceAngularVelocity, matrixInterpolation_angularVelocity );
    //auto matrixInterpolation_angularVelocity = remeshInterp.computeMatrixInterpolation( old_spaceAngularVelocity,M_spaceAngularVelocity );
    remeshInterp.registeringBlockIndex( fluidToolbox.keyword(), fluidToolbox.startSubBlockSpaceIndex( "body-bc."+this->name()+".angular-velocity" ), old_spaceAngularVelocity, M_spaceAngularVelocity );
    remeshInterp.interpolate( old_fieldAngularVelocity, M_fieldAngularVelocity );
    if ( M_bdfAngularVelocity )
        M_bdfAngularVelocity->applyRemesh( M_spaceAngularVelocity, matrixInterpolation_angularVelocity );

    if ( M_matrixPTilde_angular )
    {
        M_matrixPTilde_angular.reset();
        this->updateMatrixPTilde_angular( fluidToolbox, M_matrixPTilde_angular );
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::vector<typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition const*>
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::NBodyArticulated::bodyList( bool withMaster ) const
{
    std::string const& masterBodyBCName = this->masterBodyBC().name();
    std::vector<BodyBoundaryCondition const*> res;
    for ( auto const& ba : M_articulations )
    {
        auto const& bbc1 = ba.body1();
        if ( withMaster || (bbc1.name() != masterBodyBCName) )
        {
            if ( std::find_if( res.begin(), res.end(),
                               [&bbc1]( BodyBoundaryCondition const* e ) { return e->name() == bbc1.name(); } ) ==res.end() )
                res.push_back( &bbc1 );
        }
        auto const& bbc2 = ba.body2();
        if ( withMaster || (bbc2.name() != masterBodyBCName) )
        {
            if ( std::find_if( res.begin(), res.end(),
                               [&bbc2]( BodyBoundaryCondition const* e ) { return e->name() == bbc2.name(); } ) == res.end() )
                res.push_back( &bbc2 );
        }
    }
    return res;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::NBodyArticulated::updateForUse()
{
    auto allbbc = this->bodyList();
    M_mass = 0;
    M_massCenter = eigen_vector_type<nRealDim>::Zero();
    for ( auto const& bbc : allbbc )
    {
        M_mass += bbc->body().mass();
        M_massCenter += bbc->body().mass()*bbc->body().massCenter();
    }
    M_massCenter /= M_mass;

    M_momentOfInertia_bodyFrame = moment_of_inertia_type::Zero();
    for ( auto const& bbc : allbbc )
        bbc->body().computeMomentOfInertia_bodyFrame( this->massCenterExpr(), this->rigidRotationMatrix(), M_momentOfInertia_bodyFrame, true );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
eigen_vector_type<FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::nRealDim>
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::NBodyArticulated::evaluateRelativeRigidTranslation( BodyBoundaryCondition const& bbc, BodyBoundaryCondition const& bbcMaster ) const
{
    bool findBA = false;
    double coeff = 1;
    for ( BodyArticulation const& ba : this->articulations() )
    {
        if ( ba.body1().name() == bbc.name() && ba.body2().name() == bbcMaster.name() )
        {
            coeff = 1;
            findBA = true;
        }
        else if ( ba.body2().name() == bbc.name() && ba.body1().name() == bbcMaster.name() )
        {
            coeff = -1;
            findBA = true;
        }

        if ( !findBA )
            continue;

        auto [newMass1,newMassCenter1] = ba.body1().body().computeMassAndMassCenterFromDisplacementField( ba.body1().body().fieldDisplacement() );
        auto [newMass2,newMassCenter2] = ba.body2().body().computeMassAndMassCenterFromDisplacementField( ba.body2().body().fieldDisplacement() );

        return coeff*ba.relativeTranslationVector(newMassCenter1,newMassCenter2);
    }

    CHECK( false ) << "not found BodyArticulation related";
    return eigen_vector_type<nRealDim>::Zero();
}



FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodySetBoundaryCondition::updateForUse( self_type const& fluidToolbox )
{
    for ( auto & [name,bpbc] : *this )
        bpbc.updateForUse( fluidToolbox );

    for ( auto & nba : M_nbodyArticulated )
        nba.updateForUse();
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodySetBoundaryCondition::init( self_type const& fluidToolbox )
{
    M_internal_elasticVelocity_is_v0 = boption(_prefix=fluidToolbox.prefix(),_name="body.elastic-behavior.use-old-version");
    std::vector<BodyArticulation> articulations;
    for ( auto const& [name,bbc] : *this )
    {
        if ( bbc.articulationTranslationalVelocityExpr().empty() )
            continue;
        std::string const& bbcName = bbc.articulationTranslationalVelocityExpr().begin()->first; // WARNING : we guess that we have only one body! TODO 
        auto itFind = this->find( bbcName );
        CHECK( itFind != this->end() ) << "body not found";

        BodyArticulation ba( &bbc, &(itFind->second) );
        //ba.setTranslationalVelocityExpr( bbc.articulationTranslationalVelocityModelExpr() );
        ba.setTranslationalVelocityExpr( bbc.articulationTranslationalVelocityExpr().begin()->second );
        articulations.push_back( std::move( ba ) );
    }

    if ( articulations.size() > 0 )
    {
        std::set<int> baDone;
        while ( true )
        {
            int first=-1;
            for (int k=0;k<articulations.size();++k )
                if ( baDone.find(k) == baDone.end() )
                {
                    first = k;
                    break;
                }
            if ( first < 0 )
                break;
            NBodyArticulated nba( fluidToolbox );
            nba.addArticulation( articulations[first] );
            baDone.insert( first );
            while ( true )
            {
                int nAddedCurrently = 0;
                for (int k=0;k<articulations.size();++k )
                {
                    if ( baDone.find( k ) != baDone.end() )
                        continue;

                    if ( nba.canBeConnectedTo( articulations[k] ) )
                    {
                        nba.addArticulation( articulations[k] );
                        baDone.insert( k );
                        ++nAddedCurrently;
                    }
                }
                if ( nAddedCurrently == 0 )
                    break;
            }
            M_nbodyArticulated.push_back( std::move( nba ) );
        }
    }

    // attach nba to bbc
    for ( auto & nba : M_nbodyArticulated )
    {
        for ( auto & [name,bbc] : *this )
        {
            if ( nba.has( bbc ) )
                bbc.attachToNBodyArticulated( nba );
        }
    }

    // std::cout << "M_nbodyArticulated.size()="<<M_nbodyArticulated.size()<<std::endl;
    // for ( auto const& nba : M_nbodyArticulated )
    //     std::cout << "nba.articulations.size()="<<nba.articulations().size()<<std::endl;

    for ( auto & [name,bbc] : *this )
        bbc.init( fluidToolbox );

    for ( auto & nba : M_nbodyArticulated )
        nba.init( fluidToolbox );
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodySetBoundaryCondition::initAlgebraicFactory( self_type const& fluidToolbox, model_algebraic_factory_ptrtype algebraicFactory )
{
    if ( this->empty() )
        return;

    size_type startBlockIndexVelocity = fluidToolbox.startSubBlockSpaceIndex("velocity");

    std::set<size_type> blockIndexNotDiagIdendity;
    std::set<std::string> bodyInPMatrixArticulation;
    for ( auto const& nba : this->nbodyArticulated() )
    {
        if ( nba.articulationMethod() != "p-matrix" )
            continue;

        for ( auto const& bbcPtr : nba.bodyList() )
        {
            blockIndexNotDiagIdendity.insert( fluidToolbox.rowStartInVector() + fluidToolbox.startSubBlockSpaceIndex("body-bc."+bbcPtr->name()+".translational-velocity") );
            bodyInPMatrixArticulation.insert(bbcPtr->name());
        }
    }

    int nBlock = fluidToolbox.nBlockMatrixGraph();
    BlocksBaseSparseMatrix<double> myblockMat(nBlock,nBlock);
    for (int i=0;i<nBlock;++i)
    {
        if ( blockIndexNotDiagIdendity.find( i ) == blockIndexNotDiagIdendity.end() )
            myblockMat(i,i) = fluidToolbox.backend()->newIdentityMatrix( fluidToolbox.algebraicBlockVectorSolution()->operator()(i)->mapPtr(),
                                                                         fluidToolbox.algebraicBlockVectorSolution()->operator()(i)->mapPtr() );
    }


    for ( auto & [bpname,bbc] : *this )
    {
        size_type startBlockIndexTranslationalVelocity = fluidToolbox.startSubBlockSpaceIndex("body-bc."+bbc.name()+".translational-velocity");

        if ( bodyInPMatrixArticulation.find( bbc.name() ) == bodyInPMatrixArticulation.end() )
            myblockMat(startBlockIndexVelocity,startBlockIndexTranslationalVelocity) = bbc.matrixPTilde_translational();

        if ( !bbc.isInNBodyArticulated() )
        {
            size_type startBlockIndexAngularVelocity = fluidToolbox.startSubBlockSpaceIndex("body-bc."+bbc.name()+".angular-velocity");
            myblockMat(startBlockIndexVelocity,startBlockIndexAngularVelocity) = bbc.matrixPTilde_angular();
        }

        auto dofsBody = fluidToolbox.functionSpaceVelocity()->dofs( bbc.rangeMarkedFacesOnFluid() );
        auto matFI_Id = myblockMat(startBlockIndexVelocity,startBlockIndexVelocity);
        for ( auto dofid : dofsBody )
        {
            matFI_Id->set( dofid,dofid, 0.);
        }
        matFI_Id->close();
    }


    // p-matrix articulation
    //for ( auto const& nba : this->nbodyArticulated() )
    for ( auto & nba : M_nbodyArticulated )
    {
        nba.updateMatrixPTilde_angular( fluidToolbox );
        size_type startBlockIndexAngularVelocity = fluidToolbox.startSubBlockSpaceIndex("body-bc."+nba.name()+".angular-velocity");
        myblockMat(startBlockIndexVelocity,startBlockIndexAngularVelocity) = nba.matrixPTilde_angular();

        if ( nba.articulationMethod() != "p-matrix" )
            continue;

        std::vector<std::shared_ptr<GraphCSR>> graphs;
        for ( auto const& bbcPtr : nba.bodyList() )
            graphs.push_back( bbcPtr->matrixPTilde_translational()->graph() );
        auto newGraph =  utility_constant_functionspace::aggregateGraph( graphs, nba.dataMapPMatrixTranslationalVelocity() );

        sparse_matrix_ptrtype matrixPMatrixPTilde_TranslationalVelocity = fluidToolbox.backend()->newMatrix( nba.dataMapPMatrixTranslationalVelocity(), fluidToolbox.functionSpaceVelocity()->mapPtr(),newGraph );

        for ( auto const& bbcPtr : nba.bodyList() )
        {
            auto XhV = fluidToolbox.functionSpaceVelocity();
            OperatorInterpolationMatrixSetup matSetup( matrixPMatrixPTilde_TranslationalVelocity, Feel::SAME_NONZERO_PATTERN );
            auto opI_partTranslationalVelocity = opInterpolation( _domainSpace=bbcPtr->spaceTranslationalVelocity() ,_imageSpace=XhV,_range=bbcPtr->rangeMarkedFacesOnFluid(),_matrix=matSetup );
        }

        size_type startBlockIndexTranslationalVelocityMasterBody = fluidToolbox.startSubBlockSpaceIndex("body-bc."+nba.masterBodyBC().name()+".translational-velocity");
        myblockMat(startBlockIndexVelocity,startBlockIndexTranslationalVelocityMasterBody) = matrixPMatrixPTilde_TranslationalVelocity;

        myblockMat(startBlockIndexTranslationalVelocityMasterBody,startBlockIndexTranslationalVelocityMasterBody) = fluidToolbox.backend()->newIdentityMatrix( nba.dataMapPMatrixTranslationalVelocity(),
                                                                                                                                                               nba.masterBodyBC().spaceTranslationalVelocity()->mapPtr() );

        for ( auto const& bbcPtr : nba.bodyList(false) )
        {
            size_type startBlockIndexTranslationalVelocityOtherBody = fluidToolbox.startSubBlockSpaceIndex("body-bc."+bbcPtr->name()+".translational-velocity");
            myblockMat(startBlockIndexTranslationalVelocityOtherBody,startBlockIndexTranslationalVelocityMasterBody) = fluidToolbox.backend()->newIdentityMatrix( nba.dataMapPMatrixTranslationalVelocity(), bbcPtr->spaceTranslationalVelocity()->mapPtr() );
            myblockMat(startBlockIndexTranslationalVelocityOtherBody,startBlockIndexTranslationalVelocityOtherBody) = fluidToolbox.backend()->newZeroMatrix(  bbcPtr->spaceTranslationalVelocity()->mapPtr() , bbcPtr->spaceTranslationalVelocity()->mapPtr() );
        }
    }

    auto matP = fluidToolbox.backend()->newBlockMatrix(_block=myblockMat, _copy_values=true);

    auto matQ = fluidToolbox.backend()->newIdentityMatrix( matP->mapColPtr(), matP->mapRowPtr() );

    algebraicFactory->initSolverPtAP( matP,matQ );

    std::set<size_type> dofEliminationIdsPtAP;
    bool hasDofEliminationIdsPtAP = false;
    for ( auto const& nba : this->nbodyArticulated() )
    {
        if ( nba.articulationMethod() != "p-matrix" )
            continue;
        hasDofEliminationIdsPtAP = true;
        for ( auto const& bbcPtr : nba.bodyList(false) )
        {
            if ( bbcPtr->spaceTranslationalVelocity()->nLocalDofWithGhost() > 0 )
            {
                size_type startBlockIndexTranslationalVelocityOtherBody = fluidToolbox.startSubBlockSpaceIndex("body-bc."+bbcPtr->name()+".translational-velocity");
                auto const& basisToContainerGpTranslationalVelocity = matP->mapCol().dofIdToContainerId( startBlockIndexTranslationalVelocityOtherBody );
                for( size_type _dofId : basisToContainerGpTranslationalVelocity )
                    dofEliminationIdsPtAP.insert( _dofId );
            }
        }
    }
    if ( hasDofEliminationIdsPtAP )
        algebraicFactory->solverPtAP_setDofEliminationIds( dofEliminationIdsPtAP );

    if ( ( this->hasElasticVelocity() && false ) ||  this->hasArticulationWithMethodPMatrix() )
        algebraicFactory->initExplictPartOfSolution();
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodySetBoundaryCondition::updateAlgebraicFactoryForUse( self_type const& fluidToolbox, model_algebraic_factory_ptrtype algebraicFactory )
{
    if ( this->empty() )
        return;

    // not very nice, we need to update direclty P, not rebuild
    if ( !algebraicFactory->hasInitSolverPtAP() )
        this->initAlgebraicFactory( fluidToolbox, algebraicFactory );
    else
    {
        auto matP = algebraicFactory->solverPtAP_matrixP();
        size_type startBlockIndexVelocity = fluidToolbox.startSubBlockSpaceIndex("velocity");
        for ( auto & [bpname,bbc] : *this )
        {
            if ( true )
            {
                auto dofsBody = fluidToolbox.functionSpaceVelocity()->dofs( bbc.rangeMarkedFacesOnFluid() );
                auto const& dofIdToContainer_velocity_row = matP->mapRow().dofIdToContainerId( startBlockIndexVelocity );
                auto const& dofIdToContainer_velocity_col = matP->mapCol().dofIdToContainerId( startBlockIndexVelocity );
                double val = bbc.hasElasticVelocity() && !this->internal_elasticVelocity_is_v0() ? 1.0 : 0.0;
                for ( auto dofid :  dofsBody )
                    matP->set( dofIdToContainer_velocity_row[dofid],dofIdToContainer_velocity_col[dofid], val );
            }

            if ( bbc.isInNBodyArticulated() )
                continue;

            size_type startBlockIndexAngularVelocity = fluidToolbox.startSubBlockSpaceIndex("body-bc."+bbc.name()+".angular-velocity");
            bbc.updateMatrixPTilde_angular( fluidToolbox, matP, startBlockIndexVelocity, startBlockIndexAngularVelocity );
        }

        for ( auto const& nba : this->nbodyArticulated() )
        {
            size_type startBlockIndexAngularVelocity = fluidToolbox.startSubBlockSpaceIndex("body-bc."+nba.name()+".angular-velocity");
            nba.updateMatrixPTilde_angular( fluidToolbox, matP, startBlockIndexVelocity, startBlockIndexAngularVelocity );
        }
    }

    bool applyCloseInExplictPartOfSolution = false;
    size_type startBlockIndexVelocity = fluidToolbox.startSubBlockSpaceIndex("velocity");

    // update explicit part of solution if use articulation with p-matrix
    for ( auto const& nba : this->nbodyArticulated() )
    {
        if ( nba.articulationMethod() != "p-matrix" )
            continue;

        for ( auto const& ba : nba.articulations() )
        {
            auto const& bbc1 = ba.body1();
            auto const& bbc2 = ba.body2();
            CHECK( bbc1.name() == nba.masterBodyBC().name() || bbc2.name() == nba.masterBodyBC().name() ) << "Case not handle : too complex articulation";
            auto const& bbc = bbc1.name() == nba.masterBodyBC().name()? bbc2 : bbc1;

            applyCloseInExplictPartOfSolution = true;
            size_type startBlockIndexTranslationalVelocity = fluidToolbox.startSubBlockSpaceIndex("body-bc."+bbc.name()+".translational-velocity");

            auto articulationTranslationalVelocityExpr = ba.translationalVelocityExpr( fluidToolbox.symbolsExpr() );
            auto uExplictiPart = fluidToolbox.functionSpaceVelocity()->element( algebraicFactory->explictPartOfSolution(), fluidToolbox.rowStartInVector()+startBlockIndexVelocity);
            uExplictiPart.on(_range=bbc.rangeMarkedFacesOnFluid(),_expr=articulationTranslationalVelocityExpr,_close=true ); // TODO sync all body in one call
            if ( bbc.spaceTranslationalVelocity()->nLocalDofWithGhost() > 0 )
            {
                auto articulationTranslationalVelocityExprEvaluated = articulationTranslationalVelocityExpr.evaluate(false);
                auto const& basisToContainerGpTranslationalVelocity = algebraicFactory->explictPartOfSolution()->map().dofIdToContainerId( fluidToolbox.rowStartInVector()+startBlockIndexTranslationalVelocity );
                for (int d=0;d<nDim;++d)
                    algebraicFactory->explictPartOfSolution()->set(basisToContainerGpTranslationalVelocity[d], articulationTranslationalVelocityExprEvaluated(d) );
            }
        }
    }

    // update explicit part of solution if we have an elastic velocity
    if ( this->hasElasticVelocity() && this->internal_elasticVelocity_is_v0() )
    {
        if ( !algebraicFactory->explictPartOfSolution() )
            algebraicFactory->initExplictPartOfSolution();
        auto uExplictiPart = fluidToolbox.functionSpaceVelocity()->element( algebraicFactory->explictPartOfSolution(), fluidToolbox.rowStartInVector()+startBlockIndexVelocity);
        for ( auto const& [bname,bbc] : *this )
            uExplictiPart.on(_range=bbc.rangeMarkedFacesOnFluid(),_expr=idv(bbc.fieldElasticVelocityPtr()),_close=true ); // TODO sync all body in one call
    }

    if ( applyCloseInExplictPartOfSolution )
        algebraicFactory->explictPartOfSolution()->close();

}


} // namespace FeelModels
} // namespace Feel
