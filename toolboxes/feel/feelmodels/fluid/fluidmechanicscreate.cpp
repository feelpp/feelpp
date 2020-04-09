/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelfilters/savegmshmesh.hpp>
//#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
//#include <feel/feelvf/inv.hpp>
#include <feel/feelpde/operatorpcd.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>
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
    ModelPhysics<nDim>( "fluid" ),
    M_materialProperties( new material_properties_type( prefix ) ),
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
        M_mesh = __mesh;
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
    M_isHOVisu = nOrderGeo > 1;
    if ( Environment::vm().count(prefixvm(this->prefix(),"hovisu").c_str()) )
        M_isHOVisu = boption(_name="hovisu",_prefix=this->prefix());
    //--------------------------------------------------------------//
    M_haveSourceAdded=false;//true when update
    M_velocityDivIsEqualToZero=true;

    M_dirichletBCnitscheGamma = doption(_name="dirichletbc.nitsche.gamma",_prefix=this->prefix());

#if 0
    if ( Environment::vm().count(prefixvm(this->prefix(),"model").c_str()) )
        this->setModelName( soption(_name="model",_prefix=this->prefix()) );
#endif
    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        this->setSolverName( soption(_name="solver",_prefix=this->prefix()) );

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

    M_applyCIPStabOnlyOnBoundaryFaces=false;
    M_doCIPStabConvection = boption(_name="stabilisation-cip-convection",_prefix=this->prefix());
    M_doCIPStabDivergence = boption(_name="stabilisation-cip-divergence",_prefix=this->prefix());
    M_doCIPStabPressure = boption(_name="stabilisation-cip-pressure",_prefix=this->prefix());
    M_stabCIPConvectionGamma = doption(_name="stabilisation-cip-convection-gamma",_prefix=this->prefix());
    M_stabCIPDivergenceGamma = doption(_name="stabilisation-cip-divergence-gamma",_prefix=this->prefix());
    M_stabCIPPressureGamma = doption(_name="stabilisation-cip-pressure-gamma",_prefix=this->prefix());

    M_doStabDivDiv = boption(_name="stabilisation-div-div",_prefix=this->prefix());
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

    //--------------------------------------------------------------//
    // gravity
    std::string gravityStr;
    if ( Environment::vm().count(prefixvm(this->prefix(),"gravity-force").c_str()) )
        gravityStr = soption(_name="gravity-force",_prefix=this->prefix());
    else if (nDim == 2 )
        gravityStr = "{0,-9.80665}";
    else if (nDim == 3 )
        gravityStr = "{0,0,-9.80665}";
    M_gravityForce = expr<nDim,1,2>( gravityStr,"",this->worldComm(),this->repository().expr() );
    M_useGravityForce = boption(_name="use-gravity-force",_prefix=this->prefix());

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

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

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
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::createFunctionSpaces()
{
    this->log("FluidMechanics","createFunctionSpaces","start");
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

    // update rho, mu, nu,...
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().materials().setParameterValues( paramValues );
    M_materialProperties->updateForUse( this->mesh(), this->modelProperties().materials(), hasExtendedDofTable );

    // fluid mix space : velocity and pressure
    if ( M_materialProperties->isDefinedOnWholeMesh() )
    {
        M_rangeMeshElements = elements(this->mesh());
        M_XhVelocity = space_velocity_type::New( _mesh=M_mesh,
                                                 _extended_doftable=extendedDT[0] );
        M_XhPressure = space_pressure_type::New( _mesh=M_mesh,
                                                 _extended_doftable=extendedDT[1] );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), M_materialProperties->markers());
        M_XhVelocity = space_velocity_type::New( _mesh=M_mesh,
                                                 _extended_doftable=extendedDT[0],
                                                 _range=M_rangeMeshElements );
        M_XhPressure = space_pressure_type::New( _mesh=M_mesh,
                                                 _extended_doftable=extendedDT[1],
                                                 _range=M_rangeMeshElements );
    }

    M_fieldVelocity.reset( new element_velocity_type(M_XhVelocity,"velocity") );
    M_fieldPressure.reset( new element_pressure_type(M_XhPressure,"pressure") );

    double tElapsed = this->timerTool("Constructor").stop("createSpaces");
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

        M_meshALE = meshale( _mesh=M_mesh,_prefix=this->prefix(),_directory=this->repository() );
        this->log("FluidMechanics","createALE", "create meshale object done" );
        // mesh displacement only on moving
        M_meshDisplacementOnInterface.reset( new element_mesh_disp_type(M_meshALE->displacement()->functionSpace(),"mesh_disp_on_interface") );
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
                    BodyBoundaryCondition bpbc;
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
    if (this->hasMarkerDirichletBClm())
    {
        //std::cout << "createTraceMesh\n"<<std::endl;
        bool useSubMeshRelation = boption(_name="dirichletbc.lm.use-submesh-relation",_prefix=this->prefix());
        size_type useSubMeshRelationKey = (useSubMeshRelation)? EXTRACTION_KEEP_MESH_RELATION : 0;
        M_meshDirichletLM = createSubmesh(_mesh=this->mesh(),_range=markedfaces(this->mesh(),this->markerDirichletBClm()), _context=useSubMeshRelationKey,_view=true );
        if ( boption(_name="dirichletbc.lm.savemesh",_prefix=this->prefix()) )
        {
            std::string nameMeshDirichletLM = "nameMeshDirichletLM.msh";
            saveGMSHMesh(_mesh=M_meshDirichletLM,_filename=nameMeshDirichletLM);
        }

        M_XhDirichletLM = space_trace_velocity_type::New( _mesh=M_meshDirichletLM, _worldscomm=this->localNonCompositeWorldsComm() );
        //std::cout << "M_XhDirichletLM->nDof()"<< M_XhDirichletLM->nDof() <<std::endl;
    }

    // lagrange multiplier for pressure bc
    if ( this->hasMarkerPressureBC() )
    {
        M_meshLagrangeMultiplierPressureBC = createSubmesh(_mesh=this->mesh(),_range=markedfaces(this->mesh(),this->markerPressureBC()),_view=true );
        M_spaceLagrangeMultiplierPressureBC = space_trace_velocity_component_type::New( _mesh=M_meshLagrangeMultiplierPressureBC, _worldscomm=this->localNonCompositeWorldsComm() );
        M_fieldLagrangeMultiplierPressureBC1.reset( new element_trace_velocity_component_type( M_spaceLagrangeMultiplierPressureBC ) );
        if ( nDim == 3 )
            M_fieldLagrangeMultiplierPressureBC2.reset( new element_trace_velocity_component_type( M_spaceLagrangeMultiplierPressureBC ) );
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
    std::string geoExportType="static";//change_coords_only, change, static

    if ( nOrderGeo == 1 /*&& doExport*/ )
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
                                       _backend=M_backend,
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
                                       _backend=M_backend,
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
                                        _backend=M_backend,
                                        _type=InterpolationNonConforme(false,true,false,15) );

        this->log("FluidMechanics","createPostProcessExporters", "step1 done" );

        M_opIpressure = opInterpolation(_domainSpace=M_XhPressure,
                                        _imageSpace=M_XhScalarVisuHO,
                                        _range=elements(M_XhScalarVisuHO->mesh()),
                                        _backend=M_backend,
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
                                          _backend=M_backend,
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
                                            _backend=M_backend,
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::createFunctionSpacesNormalStress()
{
    if ( M_XhNormalBoundaryStress ) return;
    M_XhNormalBoundaryStress = space_normalstress_type::New( _mesh=M_meshTrace );
    M_fieldNormalStress.reset(new element_normalstress_type(M_XhNormalBoundaryStress));
    M_fieldWallShearStress.reset(new element_normalstress_type(M_XhNormalBoundaryStress));
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::createFunctionSpacesVorticity()
{
    if ( M_materialProperties->isDefinedOnWholeMesh() )
        M_XhVorticity = space_vorticity_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm());
    else
        M_XhVorticity = space_vorticity_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm(),
                                                   _range=M_rangeMeshElements );
    M_fieldVorticity.reset( new element_vorticity_type(M_XhVorticity));
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::createFunctionSpacesSourceAdded()
{
    if ( M_materialProperties->isDefinedOnWholeMesh() )
        M_XhSourceAdded=space_vectorial_PN_type::New( _mesh=M_mesh,_worldscomm=this->localNonCompositeWorldsComm() );
    else
        M_XhSourceAdded=space_vectorial_PN_type::New( _mesh=M_mesh,_worldscomm=this->localNonCompositeWorldsComm(),
                                                      _range=M_rangeMeshElements );
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

    if ( !M_mesh )
        this->initMesh();

    // backend
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldCommPtr() );
#if 0
    if ( M_modelName.empty() )
    {
        std::string theFluidModel = this->modelProperties().models().model( this->keyword() ).equations();
        this->setModelName( theFluidModel );
    }
#endif
    if ( M_solverName.empty() )
    {
        bool isLinear = true;
        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
            if ( physicFluidData->equation() == "Navier-Stokes" )
            {
                isLinear = false;
                break;
            }
        }
        if ( isLinear )
            M_solverName="LinearSystem";
        else
            M_solverName="Newton";
    }

#if 0
    // material properties
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().materials().setParameterValues( paramValues );
    if ( !M_materialsProperties )
    {
        M_materialsProperties.reset( new materialsproperties_type( this->prefix(), this->repository().expr() ) );
        M_materialsProperties->updateForUse( M_mesh, this->modelProperties().materials(), *this );
    }
#endif
    // functionSpaces and elements
    this->createFunctionSpaces();

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
        if ( Environment::vm().count( prefixvm(this->prefix(),"stabilization-gls.convection-diffusion.location.expressions" ) ) )
        {
            std::string locationExpression = soption(_prefix=this->prefix(),_name="stabilization-gls.convection-diffusion.location.expressions");
            auto rangeStab = elements(this->mesh(),expr(locationExpression));
            for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
                M_stabilizationGLSEltRangeConvectionDiffusion[rangeData.first] = intersect( rangeStab, rangeData.second );
        }
        else
        {
            for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
                M_stabilizationGLSEltRangeConvectionDiffusion[rangeData.first] = rangeData.second;
        }
        for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
            M_stabilizationGLSEltRangePressure[rangeData.first] = rangeData.second;
    }
    //-------------------------------------------------//
    // init function defined in json
    this->initUserFunctions();
    // init post-processinig (exporter, measure at point, ...)
    this->initPostProcess();
    //-------------------------------------------------//
    // init ALE mesh
    if (this->isMoveDomain())
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        auto itAleBC = this->markerALEMeshBC().begin();
        auto const enAleBC = this->markerALEMeshBC().end();
        for ( ; itAleBC!=enAleBC ; ++itAleBC )
        {
            std::string bcName = itAleBC->first;
            auto itAleMark = itAleBC->second.begin();
            auto const enAleMark = itAleBC->second.end();
            for ( ; itAleMark!=enAleMark ; ++itAleMark )
                M_meshALE->addBoundaryFlags( bcName, *itAleMark );
        }

        M_meshALE->init();

        this->log("FluidMechanics","init", "meshALE done" );
#endif
    }

    //-------------------------------------------------//
    // bc body (call after meshALE->init() in case of restart)
    M_bodySetBC.updateForUse( *this );

    // update constant parameters
    this->updateParameterValues();

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
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    M_blockVectorSolution.buildVector( this->backend() );

    M_algebraicFactory.reset( new model_algebraic_factory_type(this->shared_from_this(),this->backend()) );

    if ( boption(_name="use-velocity-near-null-space",_prefix=this->prefix() ) )
    {
        std::string nearNullSpacePrefix = this->prefix();
        if ( Environment::vm().count(prefixvm(this->prefix(),"use-velocity-near-null-space.prefix").c_str()) )
            nearNullSpacePrefix = soption( _name="use-velocity-near-null-space.prefix", _prefix=this->prefix() );

        NullSpace<double> userNullSpace = detail::getNullSpace(this->functionSpaceVelocity(), mpl::int_<nDim>() ) ;
        M_algebraicFactory->attachNearNullSpace( 0,userNullSpace, nearNullSpacePrefix ); // for block velocity in fieldsplit
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
            this->algebraicFactory()->preconditionerTool()->attachAuxiliarySparseMatrix( "mass-matrix", massbf.matrixPtr() );
    }

    if ( this->hasOperatorPCD() )
        this->algebraicFactory()->preconditionerTool()->attachOperatorPCD("pcd", this->operatorPCD());

    if ( M_timeStepping == "Theta" )
    {
        M_timeStepThetaSchemePreviousContrib = this->backend()->newVector(M_blockVectorSolution.vectorMonolithic()->mapPtr() );
        M_algebraicFactory->addVectorResidualAssembly( M_timeStepThetaSchemePreviousContrib, 1.0, "Theta-Time-Stepping-Previous-Contrib", true );
        M_algebraicFactory->addVectorLinearRhsAssembly( M_timeStepThetaSchemePreviousContrib, -1.0, "Theta-Time-Stepping-Previous-Contrib", false );
        if ( M_stabilizationGLS )
            M_algebraicFactory->dataInfos().addVectorInfo( prefixvm( this->prefix(),"time-stepping.previous-solution"), this->backend()->newVector( M_blockVectorSolution.vectorMonolithic()->mapPtr() ) );
    }


    if ( !M_bodySetBC.empty() )
    {
        int nBlock = this->nBlockMatrixGraph();
        BlocksBaseSparseMatrix<double> myblockMat(nBlock,nBlock);
        for (int i=0;i<nBlock;++i)
            myblockMat(i,i) = this->backend()->newIdentityMatrix( M_blockVectorSolution(i)->mapPtr(),M_blockVectorSolution(i)->mapPtr() );

        size_type startBlockIndexVelocity = this->startSubBlockSpaceIndex("velocity");
        std::set<size_type> dofsAllBodies;
        for ( auto const& [bpname,bpbc] : M_bodySetBC )
        {
            //CHECK( this->hasStartSubBlockSpaceIndex("body-bc.translational-velocity") ) << " start dof index for body-bc.translational-velocity is not present\n";
            //CHECK( this->hasStartSubBlockSpaceIndex("body-bc.angular-velocity") ) << " start dof index for body-bc.angular-velocity is not present\n";
            size_type startBlockIndexTranslationalVelocity = this->startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
            size_type startBlockIndexAngularVelocity = this->startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");

            myblockMat(startBlockIndexVelocity,startBlockIndexTranslationalVelocity) = bpbc.matrixPTilde_translational();
            myblockMat(startBlockIndexVelocity,startBlockIndexAngularVelocity) = bpbc.matrixPTilde_angular();

            auto rangeBody = bpbc.rangeMarkedFacesOnFluid();
            auto dofsBody = this->functionSpaceVelocity()->dofs( rangeBody );
            auto matFI_Id = myblockMat(startBlockIndexVelocity,startBlockIndexVelocity);
            for ( auto dofid : dofsBody )
            {
                matFI_Id->set( dofid,dofid, 0.);
                dofsAllBodies.insert( dofid );
            }
            matFI_Id->close();
        }

        auto matP = backend()->newBlockMatrix(_block=myblockMat, _copy_values=true);
        M_algebraicFactory->initSolverPtAP( matP );

        this->functionSpaceVelocity()->dof()->updateIndexSetWithParallelMissingDof( dofsAllBodies );
        std::set<size_type> dofEliminationIdsPtAP;
        matP->mapRow().dofIdToContainerId(startBlockIndexVelocity, dofsAllBodies, dofEliminationIdsPtAP );
        M_algebraicFactory->solverPtAP_setDofEliminationIds( dofEliminationIdsPtAP );

        if ( M_bodySetBC.hasElasticVelocity() )
            M_algebraicFactory->initExplictPartOfSolution();
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
            if ( nDim != 2 ) continue;
            CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
            M_fieldsUserVectorial[funcName]->on(_range=M_rangeMeshElements,_expr=funcData.expressionVectorial2() );
        }
        else if ( funcData.isVectorial3() )
        {
            if ( nDim != 3 ) continue;
            CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
            M_fieldsUserVectorial[funcName]->on(_range=M_rangeMeshElements,_expr=funcData.expressionVectorial3() );
        }
    }
}


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


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->setPostProcessExportsAllFieldsAvailable( {"velocity","pressure","vorticity","displacement"} );
    this->setPostProcessExportsPidName( "pid" );
    this->setPostProcessExportsAllFieldsAvailable( "trace_mesh", {"trace.normal-stress","trace.wall-shear-stress" /*, "trace.body.translational-velocity", "trace.body.angular-velocity"*/ } );
    this->setPostProcessExportsPidName( "trace_mesh", "trace.pid" );
    this->setPostProcessSaveAllFieldsAvailable( {"velocity","pressure","vorticity","displacement"} );
    super_type::initPostProcess();

    // init exporters
    if ( boption(_name="exporter.export") )
    {
        //M_postProcessFieldExported = this->postProcessFieldExported( this->modelProperties().postProcess().exports( this->keyword() ).fields() );
        //if ( !M_postProcessFieldExported.empty() )
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

        //M_postProcessFieldOnTraceExported = this->postProcessFieldOnTraceExported( this->modelProperties().postProcess().exports( this->keyword() ).fields() );
        //if ( !M_postProcessFieldOnTraceExported.empty() && nOrderGeo == 1 )
        if ( !this->postProcessExportsFields( "trace_mesh" ).empty() && nOrderGeo == 1  )
        {
            if ( !M_materialProperties->isDefinedOnWholeMesh() )
                this->functionSpaceVelocity()->dof()->meshSupport()->updateBoundaryInternalFaces();
#if 1
            auto rangeTrace = ( M_materialProperties->isDefinedOnWholeMesh() )? boundaryfaces(this->mesh()) : this->functionSpaceVelocity()->dof()->meshSupport()->rangeBoundaryFaces(); // not very nice, need to store the meshsupport
            M_meshTrace = createSubmesh( _mesh=this->mesh(), _range=rangeTrace, _context=size_type(EXTRACTION_KEEP_MESH_RELATION|EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT),_view=true );
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

            if ( this->hasPostProcessExportsField( "trace_mesh", "trace.normal-stress" ) ||
                 this->hasPostProcessExportsField( "trace_mesh", "trace.wall-shear-stress" ) )
                this->createFunctionSpacesNormalStress();
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

    for ( auto const& [bpname,bpbc] : M_bodySetBC )
    {
        this->setStartSubBlockSpaceIndex( "body-bc."+bpbc.name()+".translational-velocity", currentStartIndex++ );
        this->setStartSubBlockSpaceIndex( "body-bc."+bpbc.name()+".angular-velocity", currentStartIndex++ );
    }


    return currentStartIndex;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::buildBlockVector()
{
    this->initBlockVector();
    //M_blockVectorSolution.buildVector( this->backend() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
int
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initBlockVector()
{
    int nBlock = this->nBlockMatrixGraph();
    M_blockVectorSolution.resize( nBlock );
    int cptBlock = 0;
    M_blockVectorSolution(cptBlock++) = this->fieldVelocityPtr();
    M_blockVectorSolution(cptBlock++) = this->fieldPressurePtr();
    // impose mean pressure by lagrange multiplier
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        for ( int k=0;k<M_XhMeanPressureLM.size();++k )
            M_blockVectorSolution(cptBlock++) = this->backend()->newVector( M_XhMeanPressureLM[k] );
    }
    // lagrange multiplier for Dirichlet BC
    if (this->hasMarkerDirichletBClm())
    {
        M_blockVectorSolution(cptBlock) = this->backend()->newVector( this->XhDirichletLM() );
        ++cptBlock;
    }
    if ( this->hasMarkerPressureBC() )
    {
        M_blockVectorSolution(cptBlock++) = M_fieldLagrangeMultiplierPressureBC1;
        if ( nDim == 3 )
            M_blockVectorSolution(cptBlock++) = M_fieldLagrangeMultiplierPressureBC2;
    }
    // windkessel outel with implicit scheme
    if ( this->hasFluidOutletWindkesselImplicit() )
    {
        for (int k=0;k<this->nFluidOutletWindkesselImplicit();++k)
        {
            M_blockVectorSolution(cptBlock) = this->backend()->newVector( M_fluidOutletWindkesselSpace );
            ++cptBlock;
        }
    }

    if ( !M_bodySetBC.empty() )
    {
        for ( auto const& [bpname,bpbc] : M_bodySetBC )
        {
            M_blockVectorSolution(cptBlock++) = bpbc.fieldTranslationalVelocityPtr();
            M_blockVectorSolution(cptBlock++) = bpbc.fieldAngularVelocityPtr();
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
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::Body::setup( pt::ptree const& p, ModelMaterials const& mats, mesh_ptrtype mesh, std::string const& exprRepository )
{
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

    M_materialsProperties.reset( new materialsproperties_type( "",exprRepository/*this->prefix(), this->repository().expr()*/ ) );
    M_materialsProperties->updateForUse( mesh, mats, *this, matNames, onlyMarkers );

    this->updateForUse();

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::Body::updateForUse()
{
    CHECK( M_materialsProperties ) << "no materialsProperties defined";

    M_mass = 0;
    M_massCenter = eigen_vector_type<nRealDim>::Zero();
    for ( auto const& rangeData : M_materialsProperties->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& density = M_materialsProperties->density( matName );
        auto const& densityExpr = density.exprScalar();

        M_mass += integrate(_range=range,_expr=densityExpr).evaluate()(0,0);
        M_massCenter += integrate(_range=range,_expr=densityExpr*P()).evaluate();
    }
    M_massCenter /= M_mass;

    M_momentOfInertia = moment_of_inertia_type::Zero();
    for ( auto const& rangeData : M_materialsProperties->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& density = M_materialsProperties->density( matName );
        auto const& densityExpr = density.exprScalar();

        if constexpr ( nDim == 2 )
                     {
                         M_momentOfInertia(0,0) += integrate(_range=range,_expr=densityExpr*( inner(P()-this->massCenterExpr()) ) ).evaluate()(0,0);
                     }
        else
        {
            auto rvec = P()-this->massCenterExpr();
            M_momentOfInertia += integrate(_range=range,_expr=densityExpr*( inner(rvec)*eye<nDim,nDim>() - rvec*trans(rvec) ) ).evaluate();
        }
    }

}

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
        if ( M_body.physics().empty() )
            M_body.initPhysics( "body", ModelModels{}/*fluidToolbox.modelProperties().models()*/ );
        M_body.setup( *ptmaterials, fluidToolbox.modelProperties().materials(), fluidToolbox.mesh(), fluidToolbox.repository().expr() );
    }
    else
    {
        ModelExpression massExpr, momentOfInertiaExpr, initialMassCenterExpr;
        massExpr.setExpr( "mass", pt, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
        if ( massExpr.template hasExpr<1,1>() )
            M_body.setMass( massExpr.template expr<1,1>().evaluate()(0,0) );
        momentOfInertiaExpr.setExpr( "moment-of-inertia", pt, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
        if constexpr ( nDim == 2 )
        {
            if ( momentOfInertiaExpr.template hasExpr<1,1>() )
                M_body.setMomentOfInertia( momentOfInertiaExpr.template expr<1,1>().evaluate()(0,0) );
        }
        else
        {
            if ( momentOfInertiaExpr.template hasExpr<nDim,nDim>() )
                M_body.setMomentOfInertia( momentOfInertiaExpr.template expr<nDim,nDim>().evaluate() );
        }
        initialMassCenterExpr.setExpr( "mass-center", pt, fluidToolbox.worldComm(), fluidToolbox.repository().expr() /*,indexes*/ );
        if ( initialMassCenterExpr.template hasExpr<nRealDim,1>() )
        {
            auto initMassCenter = initialMassCenterExpr.template expr<nRealDim,1>();
            M_massCenterRef = initMassCenter.evaluate();
            M_body.setMassCenter( M_massCenterRef );
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
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::init( self_type const& fluidToolbox )
{
    auto rangeBodyBoundary = markedfaces(fluidToolbox.mesh(), std::set<std::string>(M_markers) );
    M_rangeMarkedFacesOnFluid = rangeBodyBoundary;
    M_mesh = createSubmesh(_mesh=fluidToolbox.mesh(),_range=rangeBodyBoundary,_view=true );
    M_XhTranslationalVelocity = space_trace_p0c_vectorial_type::New( _mesh=M_mesh );
    if constexpr ( nDim == 2 )
                     M_XhAngularVelocity = space_trace_angular_velocity_type::New( _mesh=M_mesh );
    else
        M_XhAngularVelocity = M_XhTranslationalVelocity;
    M_fieldTranslationalVelocity = M_XhTranslationalVelocity->elementPtr();
    M_fieldAngularVelocity = M_XhAngularVelocity->elementPtr();

    if ( this->hasElasticVelocityFromExpr() )
    {
        M_XhElasticVelocity = space_trace_velocity_type::New( _mesh=M_mesh );
        M_fieldElasticVelocity = M_XhElasticVelocity->elementPtr();
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::updateForUse( self_type const& fluidToolbox )
{
    if ( !M_mesh )
        this->init( fluidToolbox );

    auto const& w = *M_fieldAngularVelocity;

    if (fluidToolbox.isMoveDomain())
    {
        if ( M_body.hasMaterialsProperties() )
        {
            M_body.updateForUse();
        }
        else
        {
            auto disp = mean(_range=M_rangeMarkedFacesOnFluid,_expr=idv(fluidToolbox.meshALE()->displacement()) );
            //M_massCenter = M_massCenterRef + disp;
            M_body.setMassCenter( M_massCenterRef + disp );
        }

        //if ( fluidToolbox.worldComm().isMasterRank() )
        //    std::cout << "M_massCenter=\n " << M_body.massCenter() << std::endl;
    }

    auto XhV = fluidToolbox.functionSpaceVelocity();

    XhV->rebuildDofPoints(); // TODO REMOVE this line, interpolation operator must not use dofPoint!!!

    // matrix interpolation of translational velocity
    if ( !M_matrixPTilde_translational )
    {
        auto opI_partTranslationalVelocity = opInterpolation( _domainSpace=M_XhTranslationalVelocity ,_imageSpace=XhV,_range=M_rangeMarkedFacesOnFluid );
        M_matrixPTilde_translational = opI_partTranslationalVelocity->matPtr();
    }

    // matrix interpolation with angular velocity expr (depends on mesh position and mass center -> rebuild at each call of updateForUse)
    auto massCenter = this->massCenterExpr();
    if constexpr (nDim == 2 )
    {
        auto opI_AngularVelocity = opInterpolation( _domainSpace=M_XhAngularVelocity ,_imageSpace=XhV,_range=M_rangeMarkedFacesOnFluid,
                                                    _type= makeExprInterpolation( id(w)*vec(-Py()+massCenter(1,0),Px()-massCenter(0,0) ), nonconforming_t() ) );
        M_matrixPTilde_angular = opI_AngularVelocity->matPtr();
    }
    else
    {
        auto r = vec(Px()-massCenter(0,0),Py()-massCenter(1,0),Pz()-massCenter(2,0) );
        auto opI_AngularVelocity = opInterpolation( _domainSpace=M_XhAngularVelocity,
                                                    _imageSpace=XhV,
                                                    _range=M_rangeMarkedFacesOnFluid,
                                                    _type= makeExprInterpolation( cross(id(w),r), nonconforming_t() ) );
        M_matrixPTilde_angular = opI_AngularVelocity->matPtr();
    }


#if 0
    if ( this->hasElasticVelocityFromExpr() )
    {
        bool meshIsOnRefAtBegin = fluidToolbox.meshALE()->isOnReferenceMesh();
        if ( !meshIsOnRefAtBegin )
            fluidToolbox.meshALE()->revertReferenceMesh( false );
        for ( auto const& [bcName,eve] : M_elasticVelocityExprBC )
        {
            auto eveRange = std::get<1>( eve ).empty()? elements(this->mesh())/*bpbc.rangeMarkedFacesOnFluid()*/ : markedelements(this->mesh(),std::get<1>( eve ) );
            auto eveExpr =  std::get<0>( eve ).template expr<nDim,1>();
            M_fieldElasticVelocity->on(_range=eveRange,_expr=eveExpr,_close=true ); // TODO crash if use here markedfaces of fluid with partial mesh support
        }
        if ( !meshIsOnRefAtBegin )
            fluidToolbox.meshALE()->revertMovingMesh( false );
    }
#endif

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodyBoundaryCondition::updateElasticVelocityFromExpr(  self_type const& fluidToolbox )
{
    if ( !this->hasElasticVelocityFromExpr() )
        return;

    bool meshIsOnRefAtBegin = fluidToolbox.meshALE()->isOnReferenceMesh();
    if ( !meshIsOnRefAtBegin )
        fluidToolbox.meshALE()->revertReferenceMesh( false );
    for ( auto const& [bcName,eve] : M_elasticVelocityExprBC )
    {
        auto eveRange = std::get<1>( eve ).empty()? elements(this->mesh())/*bpbc.rangeMarkedFacesOnFluid()*/ : markedelements(this->mesh(),std::get<1>( eve ) );
        auto eveExpr =  std::get<0>( eve ).template expr<nDim,1>();
        M_fieldElasticVelocity->on(_range=eveRange,_expr=eveExpr,_close=true ); // TODO crash if use here markedfaces of fluid with partial mesh support
    }
    if ( !meshIsOnRefAtBegin )
        fluidToolbox.meshALE()->revertMovingMesh( false );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodySetBoundaryCondition::updateForUse( self_type const& fluidToolbox )
{
    for ( auto & [name,bpbc] : *this )
        bpbc.updateForUse( fluidToolbox );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::BodySetBoundaryCondition::updateAlgebraicFactoryForUse( self_type const& fluidToolbox, model_algebraic_factory_ptrtype algebraicFactory )
{
    if ( this->empty() )
        return;

    // not very nice, we need to update direclty P, not rebuild

    int nBlock = fluidToolbox.nBlockMatrixGraph();
    BlocksBaseSparseMatrix<double> myblockMat(nBlock,nBlock);
    for (int i=0;i<nBlock;++i)
        myblockMat(i,i) = fluidToolbox.backend()->newIdentityMatrix( fluidToolbox.blockVectorSolution()(i)->mapPtr(),fluidToolbox.blockVectorSolution()(i)->mapPtr() );

    size_type startBlockIndexVelocity = fluidToolbox.startSubBlockSpaceIndex("velocity");
    for ( auto & [bpname,bpbc] : *this )
    {
        //CHECK( this->hasStartSubBlockSpaceIndex("body-bc.translational-velocity") ) << " start dof index for body-bc.translational-velocity is not present\n";
        //CHECK( this->hasStartSubBlockSpaceIndex("body-bc.angular-velocity") ) << " start dof index for body-bc.angular-velocity is not present\n";
        size_type startBlockIndexTranslationalVelocity = fluidToolbox.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
        size_type startBlockIndexAngularVelocity = fluidToolbox.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");

        myblockMat(startBlockIndexVelocity,startBlockIndexTranslationalVelocity) = bpbc.matrixPTilde_translational();
        myblockMat(startBlockIndexVelocity,startBlockIndexAngularVelocity) = bpbc.matrixPTilde_angular();

        auto dofsBody = fluidToolbox.functionSpaceVelocity()->dofs( bpbc.rangeMarkedFacesOnFluid() );
        auto matFI_Id = myblockMat(startBlockIndexVelocity,startBlockIndexVelocity);
        for ( auto dofid : dofsBody )
            matFI_Id->set( dofid,dofid, 0.);
        matFI_Id->close();
    }

    auto matP = fluidToolbox.backend()->newBlockMatrix(_block=myblockMat, _copy_values=true);
    algebraicFactory->initSolverPtAP( matP );


    if ( this->hasElasticVelocity() )
    {
        auto uExplictiPart = fluidToolbox.functionSpaceVelocity()->element( algebraicFactory->explictPartOfSolution(), fluidToolbox.rowStartInVector()+0);
        for ( auto const& [bpname,bpbc] : *this )
            uExplictiPart.on(_range=bpbc.rangeMarkedFacesOnFluid(),_expr=idv(bpbc.fieldElasticVelocityPtr()),_close=true ); // TODO sync all body in one call
    }

}


} // namespace FeelModels
} // namespace Feel
