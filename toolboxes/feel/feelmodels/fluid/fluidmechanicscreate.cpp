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
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::FluidMechanics( std::string const& prefix,
                                                    bool buildMesh,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    ModelBaseRepository const& modelRep )
    :
    super_type( prefix,worldComm,subPrefix, modelRep ),
    M_isUpdatedForUse(false ),
    M_materialProperties( new material_properties_type( prefix ) )
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

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"FluidMechanics.info") );
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // build  mesh, space,exporter,...
    if ( buildMesh )
        this->initMesh();
    //-----------------------------------------------------------------------------//

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());

}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::self_ptrtype
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::New( std::string const& prefix, bool buildMesh,
                                         WorldComm const& worldComm, std::string const& subPrefix,
                                         ModelBaseRepository const& modelRep )
{
    return boost::make_shared<self_type>( prefix, buildMesh, worldComm, subPrefix, modelRep );

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



namespace detailbc
{
std::list<std::string>
generateMarkerBCList(BoundaryConditions const bc,std::string const& field,std::string const& bcname, std::string const& marker, std::string const& param="number" )
{
    std::list<std::string> markerList;

    std::pair<bool,int> numberOfMarkerRead = bc.iparam( field/*"velocity"*/, bcname/*"Dirichlet"*/, marker/*(d)*/, param/*"number"*/ );
    int numberOfMarker = ( numberOfMarkerRead.first )? numberOfMarkerRead.second : 1;
    for (int k=0 ; k<numberOfMarker ; ++k )
    {
        std::string currentMarker = ( numberOfMarker == 1 )? marker : (boost::format("%1%%2%")%marker %k).str();
        markerList.push_back( currentMarker );
    }

    return markerList;
}
}

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

    if ( Environment::vm().count(prefixvm(this->prefix(),"model").c_str()) )
        this->setModelName( soption(_name="model",_prefix=this->prefix()) );

    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        this->setSolverName( soption(_name="solver",_prefix=this->prefix()) );

    //--------------------------------------------------------------//
    // fsi options
    M_useFSISemiImplicitScheme = false;
    M_couplingFSIcondition = "dirichlet-neumann";
    M_couplingFSI_RNG_useInterfaceOperator = false;

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
    std::vector<bool> extendedDT( space_fluid_type::nSpaces,false );
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
        M_Xh = space_fluid_type::New( _mesh=M_mesh,
                                      _extended_doftable=extendedDT );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), M_materialProperties->markers());
        M_Xh = space_fluid_type::New( _mesh=M_mesh,
                                      _extended_doftable=extendedDT, _range=M_rangeMeshElements );
    }

    M_Solution.reset( new element_fluid_type(M_Xh,"U"));

    double tElapsed = this->timerTool("Constructor").stop("createSpaces");
    this->log("FluidMechanics","createFunctionSpaces", (boost::format("finish in %1% s") %tElapsed).str() );
}

//---------------------------------------------------------------------------------------------------------//


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::createALE()
{
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() )
    {
        this->log("FluidMechanics","createALE", "start" );
        this->timerTool("Constructor").start();

        M_XhMeshVelocityInterface = space_meshvelocityonboundary_type::New(_mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm());
        M_XhMeshALEmapDisc = space_alemapdisc_type::New(_mesh=this->mesh(), _worldscomm=this->localNonCompositeWorldsComm() );

#if 0
        // init
        bool moveGhostEltFromExtendedStencil=false;
        for ( bool hasExt : M_Xh->extendedDofTableComposite() )
            moveGhostEltFromExtendedStencil = moveGhostEltFromExtendedStencil || hasExt;
#else
        bool moveGhostEltFromExtendedStencil = this->useExtendedDofTable();
#endif
        if ( moveGhostEltFromExtendedStencil )
            this->log("FluidMechanics","createALE", "use moveGhostEltFromExtendedStencil" );

        M_meshALE.reset(new mesh_ale_type( M_mesh,
                                           this->prefix(),
                                           this->localNonCompositeWorldsComm()[0],
                                           moveGhostEltFromExtendedStencil,
                                           this->repository() ));
        this->log("FluidMechanics","createALE", "--1--" );
        // mesh displacement only on moving
        M_meshDisplacementOnInterface.reset( new element_mesh_disp_type(M_meshALE->displacement()->functionSpace(),"mesh_disp_on_interface") );
        this->log("FluidMechanics","createALE", "--2--" );
        // mesh velocity only on moving interface
        M_meshVelocityInterface.reset(new element_meshvelocityonboundary_type( M_XhMeshVelocityInterface, "mesh_velocity_interface" ) );

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
        std::pair<bool,std::string> dirichletbcTypeRead = this->modelProperties().boundaryConditions().sparam( "velocity", "Dirichlet", marker(d), "method" );
        std::string dirichletbcType = ( dirichletbcTypeRead.first )? dirichletbcTypeRead.second : soption(_name="dirichletbc.type",_prefix=this->prefix());
        CHECK( dirichletbcType=="elimination" || dirichletbcType=="nitsche" || dirichletbcType=="lm" ) << "invalid dirichletbc.type " << dirichletbcType;

        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "velocity", "Dirichlet", marker(d) );
        this->setMarkerDirichletBCByNameId( dirichletbcType, marker(d), markerList,ComponentType::NO_COMPONENT );

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Dirichlet", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
    }
    for ( ComponentType comp : std::vector<ComponentType>( { ComponentType::X, ComponentType::Y, ComponentType::Z } ) )
    {
        std::string compTag = ( comp ==ComponentType::X )? "x" : (comp == ComponentType::Y )? "y" : "z";
        std::string bcDirichletCompField = (boost::format("velocity_%1%")%compTag).str();
        std::string bcDirichletCompKeyword = "Dirichlet";
        this->M_bcDirichletComponents[comp] = this->modelProperties().boundaryConditions().getScalarFields( { { bcDirichletCompField, bcDirichletCompKeyword } } );
        for( auto const& d : this->M_bcDirichletComponents.find(comp)->second )
        {
            std::pair<bool,std::string> dirichletbcTypeRead = this->modelProperties().boundaryConditions().sparam( bcDirichletCompField, bcDirichletCompKeyword, marker(d), "method" );
            std::string dirichletbcType = ( dirichletbcTypeRead.first )? dirichletbcTypeRead.second : soption(_name="dirichletbc.type",_prefix=this->prefix());
            CHECK( dirichletbcType=="elimination" || dirichletbcType=="nitsche" || dirichletbcType=="lm" ) << "invalid dirichletbc.type " << dirichletbcType;

            std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), bcDirichletCompField, bcDirichletCompKeyword, marker(d) );
            this->setMarkerDirichletBCByNameId( dirichletbcType, marker(d), markerList, comp );

            std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( bcDirichletCompField, bcDirichletCompKeyword, marker(d), "alemesh_bc" );
            std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
            for (std::string const& currentMarker : markerList )
                this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
        }
    }

    for( std::string const& bcMarker : this->modelProperties().boundaryConditions().markers( { { "velocity", "interface_fsi" }, { "velocity","moving_boundary"} } ) )
    {
        this->addMarkerALEMeshBC("moving",bcMarker);
        this->M_isMoveDomain=true;
    }

    this->M_bcNeumannScalar = this->modelProperties().boundaryConditions().getScalarFields( "velocity", "Neumann_scalar" );
    for( auto const& d : this->M_bcNeumannScalar )
    {
        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "velocity", "Neumann_scalar", marker(d) );
        this->setMarkerNeumannBC(NeumannBCShape::SCALAR,marker(d),markerList);

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Neumann_scalar", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
    }
    this->M_bcNeumannVectorial = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "velocity", "Neumann_vectorial" );
    for( auto const& d : this->M_bcNeumannVectorial )
    {
        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "velocity", "Neumann_vectorial", marker(d) );
        this->setMarkerNeumannBC(NeumannBCShape::VECTORIAL,marker(d),markerList);

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Neumann_vectorial", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
    }
    this->M_bcNeumannTensor2 = this->modelProperties().boundaryConditions().template getMatrixFields<nDim>( "velocity", "Neumann_tensor2" );
    for( auto const& d : this->M_bcNeumannTensor2 )
    {
        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "velocity", "Neumann_tensor2", marker(d) );
        this->setMarkerNeumannBC(NeumannBCShape::TENSOR2,marker(d),markerList);

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "velocity", "Neumann_tensor2", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
    }

    this->M_bcPressure = this->modelProperties().boundaryConditions().getScalarFields( "pressure", "Dirichlet" );
    for( auto const& d : this->M_bcPressure )
    {
        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "pressure", "Dirichlet", marker(d) );
        this->setMarkerPressureBC(marker(d),markerList);

        std::pair<bool,std::string> bcTypeMeshALERead = this->modelProperties().boundaryConditions().sparam( "pressure", "Dirichlet", marker(d), "alemesh_bc" );
        std::string bcTypeMeshALE = ( bcTypeMeshALERead.first )? bcTypeMeshALERead.second : std::string("fixed");
        for (std::string const& currentMarker : markerList )
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
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

        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "fluid", "outlet", bcMarker );
        for (std::string const& currentMarker : markerList )
        {
            this->M_fluidOutletsBCType.push_back(std::make_tuple(currentMarker,typeOutlet, windkesselParam ));
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
        }
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

        std::list<std::string> markerList = detailbc::generateMarkerBCList( this->modelProperties().boundaryConditions(), "fluid", "inlet", bcMarker );
        for (std::string const& currentMarker : markerList )
        {
            this->M_fluidInletDesc.push_back(std::make_tuple(currentMarker,fullTypeInlet, expr<2>( exprFluidInlet,"",this->worldComm(),this->repository().expr() )) );
            this->addMarkerALEMeshBC(bcTypeMeshALE,currentMarker);
        }
    }

    this->M_volumicForcesProperties = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "fluid", "VolumicForces" );


    // Dirichlet bc using a lagrange multiplier
    if (this->hasMarkerDirichletBClm())
    {
        //std::cout << "createTraceMesh\n"<<std::endl;
        bool useSubMeshRelation = boption(_name="dirichletbc.lm.use-submesh-relation",_prefix=this->prefix());
        size_type useSubMeshRelationKey = (useSubMeshRelation)? EXTRACTION_KEEP_MESH_RELATION : 0;
        M_meshDirichletLM = createSubmesh(this->mesh(),markedfaces(this->mesh(),this->markerDirichletBClm()), useSubMeshRelationKey );
        if ( boption(_name="dirichletbc.lm.savemesh",_prefix=this->prefix()) )
        {
            std::string nameMeshDirichletLM = "nameMeshDirichletLM.msh";
            saveGMSHMesh(_mesh=M_meshDirichletLM,_filename=nameMeshDirichletLM);
        }

        M_XhDirichletLM = space_dirichletlm_velocity_type::New( _mesh=M_meshDirichletLM, _worldscomm=this->localNonCompositeWorldsComm() );
        //std::cout << "M_XhDirichletLM->nDof()"<< M_XhDirichletLM->nDof() <<std::endl;
    }

    // lagrange multiplier for pressure bc
    if ( this->hasMarkerPressureBC() )
    {
        M_meshLagrangeMultiplierPressureBC = createSubmesh(this->mesh(),markedfaces(this->mesh(),this->markerPressureBC()) );
        M_spaceLagrangeMultiplierPressureBC = space_trace_velocity_component_type::New( _mesh=M_meshLagrangeMultiplierPressureBC, _worldscomm=this->localNonCompositeWorldsComm() );
        M_fieldLagrangeMultiplierPressureBC1.reset( new element_trace_velocity_component_type( M_spaceLagrangeMultiplierPressureBC ) );
        if ( nDim == 3 )
            M_fieldLagrangeMultiplierPressureBC2.reset( new element_trace_velocity_component_type( M_spaceLagrangeMultiplierPressureBC ) );
    }

    // init fluid outlet
    this->initFluidOutlet();
    // init fluid inlet
    this->initFluidInlet();
}
//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::createPostProcessExporters()
{
    this->log("FluidMechanics","createPostProcessExporters", "start" );

    bool doExport = boption(_name="exporter.export");
    //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;//(this->isMoveDomain())?ExporterGeometry::EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY:ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
    std::string geoExportType="static";//change_coords_only, change, static

    if ( nOrderGeo == 1 && doExport )
    {
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export",
                               //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export")),
                               _geo=geoExportType,
                               _worldcomm=M_Xh->worldComm(),
                               _path=this->exporterPath() );
    }

    if ( M_isHOVisu && doExport )
    {
#if 1 //defined(FEELPP_HAS_VTK)
        //M_exporter_ho = export_ho_type::New( this->application()->vm(), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Export_HO"))/*.c_str()*/, M_Xh->worldComm() );

// #if defined( FEELPP_MODELS_HAS_MESHALE )
//         if (M_isMoveDomain) this->meshALE()->revertReferenceMesh();
// #endif
        //auto Xh_create_ho = space_create_ho_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );

        boost::shared_ptr<mesh_visu_ho_type> meshVisuHO;
        std::string hovisuSpaceUsed = soption(_name="hovisu.space-used",_prefix=this->prefix());
        bool doLagP1parallel=false;
        if ( hovisuSpaceUsed == "velocity" )
        {
            // with velocity field
            auto Xh_create_ho = M_Xh->template functionSpace<0>()->compSpace();
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
            auto Xh_create_ho = M_Xh->template functionSpace<1>();
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
                                  _worldcomm=M_Xh->worldComm(),
                                  _path=this->exporterPath() );


        M_XhVectorialVisuHO = space_vectorial_visu_ho_type::New(_mesh=meshVisuHO/*opLagP1->mesh()*/, _worldscomm=this->localNonCompositeWorldsComm());
        //M_XhScalarVisuHO = space_scalar_visu_ho_type::New(_mesh=opLagP1->mesh(),_worldscomm=this->localNonCompositeWorldsComm());
        M_XhScalarVisuHO = M_XhVectorialVisuHO->compSpace();

        M_velocityVisuHO.reset( new element_vectorial_visu_ho_type(M_XhVectorialVisuHO,"u_visuHO"));
        M_pressureVisuHO.reset( new element_scalar_visu_ho_type(M_XhScalarVisuHO,"p_visuHO"));
        if (M_isMoveDomain) M_meshdispVisuHO.reset( new element_vectorial_visu_ho_type(M_XhVectorialVisuHO,"meshdisp_visuHO"));

        this->log("FluidMechanics","createPostProcessExporters", "start opInterpolation" );
        boost::mpi::timer timerOpI;

        M_opIvelocity = opInterpolation(_domainSpace=M_Xh->template functionSpace<0>(),
                                        _imageSpace=M_XhVectorialVisuHO,
                                        _range=elements(M_XhVectorialVisuHO->mesh()),
                                        _backend=M_backend,
                                        _type=InterpolationNonConforme(false,true,false,15) );

        this->log("FluidMechanics","createPostProcessExporters", "step1 done" );

        M_opIpressure = opInterpolation(_domainSpace=M_Xh->template functionSpace<1>(),
                                        _imageSpace=M_XhScalarVisuHO,
                                        _range=elements(M_XhScalarVisuHO->mesh()),
                                        _backend=M_backend,
                                        _type=InterpolationNonConforme(false,true,false,15) );

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

        this->log("FluidMechanics","createPostProcessExporters", "step2 done" );

        if (M_isMoveDomain )
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            std::vector<int> saveActivities_meshALE;
            if (!M_Xh->hasEntriesForAllSpaces())
            {
                saveActivities_meshALE = M_meshALE->functionSpace()->worldComm().activityOnWorld();
                M_meshALE->functionSpace()->worldComm().applyActivityOnlyOn(0/*VelocityWorld*/);
            }
            M_opImeshdisp = opInterpolation(_domainSpace=M_meshALE->functionSpace(),
                                            _imageSpace=M_XhVectorialVisuHO,
                                            _range=elements(M_XhVectorialVisuHO->mesh()),
                                            _backend=M_backend,
                                            _type=InterpolationNonConforme(false,true,false,15) );
            if (!M_Xh->hasEntriesForAllSpaces())
                M_meshALE->functionSpace()->worldComm().setIsActive(saveActivities_meshALE);
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

    if ( M_materialProperties->isDefinedOnWholeMesh() )
        M_XhNormalBoundaryStress = space_stress_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    else
        M_XhNormalBoundaryStress = space_stress_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm(), _range=M_rangeMeshElements );
    M_fieldNormalStress.reset(new element_stress_type(M_XhNormalBoundaryStress));
    M_fieldNormalStressRefMesh.reset(new element_stress_type(M_XhNormalBoundaryStress));
#if defined( FEELPP_MODELS_HAS_MESHALE )
    M_normalStressFromStruct.reset(new element_stress_type(M_XhNormalBoundaryStress));
#endif
    M_fieldWallShearStress.reset(new element_stress_type(M_XhNormalBoundaryStress));
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
        auto meshinlet = createSubmesh( this->mesh(),markedfaces(this->mesh(),marker) );
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

            auto backendinlet = backend_type::build( soption( _name="backend" ), prefixvm(this->prefix(),"fluidinlet"), M_Xh->worldComm() );
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
        double evalExprFluidInlet = exprFluidInlet.evaluate();

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
    if ( M_isUpdatedForUse ) return;

    this->log("FluidMechanics","init", "start" );
    this->timerTool("Constructor").start();

    boost::timer thetimer;

    if ( !M_mesh )
        this->initMesh();

    // backend
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), this->worldComm() );

    if ( M_modelName.empty() )
    {
        std::string theFluidModel = this->modelProperties().models().model("fluid").equations();
        this->setModelName( theFluidModel );
    }
    if ( M_solverName.empty() )
    {
        if ( M_modelName == "Stokes" || M_modelName == "StokesTransient" )
            M_solverName="LinearSystem";
        else
            M_solverName="Newton";
    }

    // functionSpaces and elements
    this->createFunctionSpaces();

    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();

    this->initBoundaryConditions();

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
            M_stabilizationGLSEltRangeConvectionDiffusion = elements(this->mesh(),expr(locationExpression));
        }
        else
        {
            M_stabilizationGLSEltRangeConvectionDiffusion = elements(this->mesh());
        }
        M_stabilizationGLSEltRangePressure = elements(this->mesh());
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

        // if restart else move submesh define from fluid mesh
        if ( this->hasFluidOutletWindkesselImplicit() )
        {
            M_fluidOutletWindkesselSpaceMeshDisp = space_fluidoutlet_windkessel_mesh_disp_type::New( _mesh=M_fluidOutletWindkesselMesh,
                                                                                                     _worldscomm=this->localNonCompositeWorldsComm() );
            M_fluidOutletWindkesselMeshDisp = M_fluidOutletWindkesselSpaceMeshDisp->elementPtr();
            M_fluidOutletWindkesselOpMeshDisp = opInterpolation(_domainSpace=M_meshALE->functionSpace(),
                                                                _imageSpace=M_fluidOutletWindkesselSpaceMeshDisp,
                                                                _range=elements(M_fluidOutletWindkesselMesh),
                                                                _backend=M_backend );
            if ( this->doRestart() )
            {
                // interpolate disp
                M_fluidOutletWindkesselOpMeshDisp->apply( *M_meshALE->displacement(), *M_fluidOutletWindkesselMeshDisp );
                // apply disp
                M_fluidOutletWindkesselMeshMover.apply( M_fluidOutletWindkesselMesh, *M_fluidOutletWindkesselMeshDisp );
            }
        }

        // space usefull to tranfert sigma*N()
        this->createFunctionSpacesNormalStress();

        this->log("FluidMechanics","init", "meshALE done" );
#endif
    }

    // call here because need meshale markers
    this->updateBoundaryConditionsForUse();

    //-------------------------------------------------//
    // define start dof index ( lm , windkessel )
    this->initStartBlockIndexFieldsInMatrix();
    //-------------------------------------------------//
    // build solution block vector
    this->buildBlockVector();

    //-------------------------------------------------//
    if ( buildModelAlgebraicFactory )
    {
        this->initAlgebraicFactory();
    }

    //-------------------------------------------------//
    //-------------------------------------------------//
    M_isUpdatedForUse = true;

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("FluidMechanics","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    M_algebraicFactory.reset( new model_algebraic_factory_type(this->shared_from_this(),this->backend()) );

    if ( boption(_name="use-velocity-near-null-space",_prefix=this->prefix() ) )
    {
        NullSpace<double> userNullSpace = detail::getNullSpace(this->functionSpaceVelocity(), mpl::int_<nDim>() ) ;
        M_algebraicFactory->attachNearNullSpace( 0,userNullSpace ); // for block velocity in fieldsplit
    }
    this->initInHousePreconditioner();
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

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
         suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_bdf_fluid = bdf(  _space=M_Xh,
                       _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity-pressure"+suffixName)),
                       _prefix=this->prefix(),
                       // don't use the fluid.bdf {initial,final,step}time but the general bdf info, the order will be from fluid.bdf
                       _initial_time=this->timeInitial(),
                       _final_time=this->timeFinal(),
                       _time_step=this->timeStep(),
                       _restart=this->doRestart(),
                       _restart_path=this->restartPath(),
                       _restart_at_last_save=this->restartAtLastSave(),
                       _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
    M_bdf_fluid->setfileFormat( myFileFormat );
    M_bdf_fluid->setPathSave( (fs::path(this->rootRepository()) /
                               fs::path( prefixvm(this->prefix(), (boost::format("bdf_o_%1%_dt_%2%")%M_bdf_fluid->bdfOrder() %this->timeStep()).str() ) ) ).string() );

    // start or restart time step scheme
    if ( !this->doRestart() )
    {
        // start time step
        M_bdf_fluid->start(*M_Solution);
        // up current time
        this->updateTime( M_bdf_fluid->time() );
    }
    else
    {
        // start time step
        M_bdf_fluid->restart();
        // load a previous solution as current solution
        *M_Solution = M_bdf_fluid->unknown(0);
        // up initial time
        this->setTimeInitial( M_bdf_fluid->timeInitial() );
        // up current time
        this->updateTime( M_bdf_fluid->time() );

        this->log("FluidMechanics","initTimeStep", "restart bdf/exporter done" );
    }

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

        M_fluidOutletWindkesselMesh = createSubmesh( this->mesh(), markedfaces(this->mesh(),markerNameBFOutletForSubmesh) );
        M_fluidOutletWindkesselSpace = space_fluidoutlet_windkessel_type::New( _mesh=M_fluidOutletWindkesselMesh,
                                                                               _worldscomm=std::vector<WorldComm>(2,this->worldComm()) );
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
                int askedIter = M_bdf_fluid->iteration() - 1;
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
    if ( this->modelProperties().functions().empty() )
        return;

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
            M_fieldsUserScalar[funcName]->on(_range=elements(this->mesh()),_expr=funcData.expressionScalar() );
        }
        else if ( funcData.isVectorial2() )
        {
            if ( nDim != 2 ) continue;
            CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
            M_fieldsUserVectorial[funcName]->on(_range=elements(this->mesh()),_expr=funcData.expressionVectorial2() );
        }
        else if ( funcData.isVectorial3() )
        {
            if ( nDim != 3 ) continue;
            CHECK( this->hasFieldUserVectorial( funcName ) ) << "user function " << funcName << "not registered";
            M_fieldsUserVectorial[funcName]->on(_range=elements(this->mesh()),_expr=funcData.expressionVectorial3() );
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
        if ( o == prefixvm(prefix,"normal-stress") || o == prefixvm(prefix,"all") )
            res.insert( "normal-stress" );
        if ( o == prefixvm(prefix,"wall-shear-stress") || o == prefixvm(prefix,"all") )
            res.insert( "wall-shear-stress" );
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
    std::string modelName = "fluid";

    // update post-process expression
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    bool hasMeasure = false;

    M_postProcessFieldExported = this->postProcessFieldExported( this->modelProperties().postProcess().exports( modelName ).fields() );
    // init exporter
    if ( !M_postProcessFieldExported.empty() )
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

    // forces (lift, drag) and flow rate measures
    pt::ptree ptree = this->modelProperties().postProcess().pTree( modelName );
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
                    this->postProcessMeasuresIO().setMeasure("drag_"+name,0.);
                    this->postProcessMeasuresIO().setMeasure("lift_"+name,0.);
                    hasMeasure = true;
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
                    this->postProcessMeasuresIO().setMeasure("flowrate_"+name,0.);
                    hasMeasure = true;
                }
            }
            else if ( ptreeLevel1Name == "Pressure" )
            {
                // this->modelProperties().postProcess().operator[](ppTypeMeasures).push_back( "Pressure" );
                M_postProcessMeasuresFields["pressure"] = "";
                this->postProcessMeasuresIO().setMeasure("pressure_sum",0.);
                this->postProcessMeasuresIO().setMeasure("pressure_mean",0.);
                hasMeasure = true;
            }
            else if ( ptreeLevel1Name == "VelocityDivergence" )
            {
                M_postProcessMeasuresFields["velocity-divergence"] = "";
                // this->modelProperties().postProcess().operator[](ppTypeMeasures).push_back( "VelocityDivergence" );
                this->postProcessMeasuresIO().setMeasure("velocity_divergence_sum",0.);
                this->postProcessMeasuresIO().setMeasure("velocity_divergence_mean",0.);
                this->postProcessMeasuresIO().setMeasure("velocity_divergence_normL2",0.);
                hasMeasure = true;
            }
        }
    }

    // point measures
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( modelName ) )
    {
        auto const& ptPos = evalPoints.pointPosition();
        node_type ptCoord(3);
        for ( int c=0;c<3;++c )
            ptCoord[c]=ptPos.value()(c);

        auto const& fields = evalPoints.fields();
        for ( std::string const& field : fields )
        {
            if ( field == "velocity" )
            {
                if ( !M_postProcessMeasuresContextVelocity )
                    M_postProcessMeasuresContextVelocity.reset( new context_velocity_type( functionSpaceVelocity()->context() ) );
                int ctxId = M_postProcessMeasuresContextVelocity->nPoints();
                M_postProcessMeasuresContextVelocity->add( ptCoord );
                std::string ptNameExport = (boost::format("velocity_%1%")%ptPos.name()).str();
                this->postProcessMeasuresEvaluatorContext().add("velocity", ctxId, ptNameExport );

                std::vector<double> vecValues = { 0. };
                if ( nDim > 1 ) vecValues.push_back( 0. );
                if ( nDim > 2 ) vecValues.push_back( 0. );
                this->postProcessMeasuresIO().setMeasureComp( ptNameExport, vecValues );
                hasMeasure = true;
            }
            else if ( field == "pressure" )
            {
                if ( !M_postProcessMeasuresContextPressure )
                    M_postProcessMeasuresContextPressure.reset( new context_pressure_type( functionSpacePressure()->context() ) );
                int ctxId = M_postProcessMeasuresContextPressure->nPoints();
                M_postProcessMeasuresContextPressure->add( ptCoord );
                std::string ptNameExport = (boost::format("pressure_%1%")%ptPos.name()).str();
                this->postProcessMeasuresEvaluatorContext().add("pressure", ctxId, ptNameExport );

                this->postProcessMeasuresIO().setMeasure(ptNameExport,0.);
                hasMeasure = true;
            }
        }
    }

    if ( hasMeasure )
    {
        this->postProcessMeasuresIO().setParameter( "time", this->timeInitial() );
        // start or restart measure file
        if (!this->doRestart())
            this->postProcessMeasuresIO().start();
        else if ( !this->isStationary() )
            this->postProcessMeasuresIO().restart( "time", this->timeInitial() );
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
size_type
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initStartBlockIndexFieldsInMatrix()
{
    size_type currentStartIndex = 2;// velocity and pressure before
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        M_startBlockIndexFieldsInMatrix["define-pressure-cst-lm"] = currentStartIndex;
        currentStartIndex += M_XhMeanPressureLM.size();
    }
    if (this->hasMarkerDirichletBClm())
    {
        M_startBlockIndexFieldsInMatrix["dirichletlm"] = currentStartIndex++;
    }
    if ( this->hasMarkerPressureBC() )
    {
        M_startBlockIndexFieldsInMatrix["pressurelm1"] = currentStartIndex++;
        if ( nDim == 3 )
            M_startBlockIndexFieldsInMatrix["pressurelm2"] = currentStartIndex++;
    }
    if ( this->hasFluidOutletWindkesselImplicit() )
    {
        M_startBlockIndexFieldsInMatrix["windkessel"] = currentStartIndex++;
    }

    return currentStartIndex;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::buildBlockVector()
{
    this->initBlockVector();
    M_blockVectorSolution.buildVector( this->backend() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
int
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initBlockVector()
{
    int nBlock = this->nBlockMatrixGraph();
    M_blockVectorSolution.resize( nBlock );
    M_blockVectorSolution(0) = this->fieldVelocityPressurePtr();
    int cptBlock=1;
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

    return cptBlock;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::initInHousePreconditioner()
{

    bool attachMassMatrix = boption(_prefix=this->prefix(),_name="preconditioner.attach-mass-matrix");
    if ( attachMassMatrix )
    {
        auto massbf = form2( _trial=this->functionSpaceVelocity(), _test=this->functionSpaceVelocity());
        auto const& u = this->fieldVelocity();
        if ( this->isStationaryModel() )
            massbf += integrate( _range=elements( this->mesh() ), _expr=inner( idt(u),id(u) ) );
        else
        {
            //double coeff = this->materialProperties()->cstRho()*this->timeStepBDF()->polyDerivCoefficient(0);
            auto coeff = idv(this->materialProperties()->fieldDensity())*this->timeStepBDF()->polyDerivCoefficient(0);
            massbf += integrate( _range=elements( this->mesh() ), _expr=coeff*inner( idt(u),id(u) ) );
        }
        massbf.matrixPtr()->close();
        this->algebraicFactory()->preconditionerTool()->attachAuxiliarySparseMatrix( "mass-matrix", massbf.matrixPtr() );
    }

    if ( M_preconditionerAttachPCD )
    {
        BoundaryConditions bcPrecPCD;
        bcPrecPCD.clear();

        auto itFindFieldVelocity = this->modelProperties().boundaryConditions().find("velocity");
        bool hasFindFieldVelocity = itFindFieldVelocity != this->modelProperties().boundaryConditions().end();
        if ( hasFindFieldVelocity )
        {
            auto itFindDirichletType = itFindFieldVelocity->second.find("Dirichlet");
            if ( itFindDirichletType != itFindFieldVelocity->second.end() )
            {
                for ( auto const& myBcDesc : itFindDirichletType->second )
                {
                    auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",myBcDesc.marker() ) );
                    auto const& listMarkerFaces = std::get<0>( ret );
                    ExpressionStringAtMarker myBcDesc2( myBcDesc );
                    myBcDesc2.setMeshMarkers( listMarkerFaces );
                    bcPrecPCD["velocity"]["Dirichlet"].push_back( myBcDesc2 );
                }
            }
            // For weak Dirichlet (Nitche,Magrange Multiplier ) ???
            // TODO Dirchlet component

            auto itFindNeumannScalType = itFindFieldVelocity->second.find("Neumann_scalar");
            if ( itFindNeumannScalType != itFindFieldVelocity->second.end() )
            {
                for ( auto const& myBcDesc : itFindNeumannScalType->second )
                {
                    auto markList = this->markerNeumannBC( NeumannBCShape::SCALAR,myBcDesc.marker() );
                    if ( markList.empty() ) continue;
                    ExpressionStringAtMarker myBcDesc2( myBcDesc );
                    myBcDesc2.setMeshMarkers( markList );
                    bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc2 );
                }
            }
            auto itFindNeumannVecType = itFindFieldVelocity->second.find("Neumann_vectorial");
            if ( itFindNeumannVecType != itFindFieldVelocity->second.end() )
            {
                for ( auto const& myBcDesc : itFindNeumannVecType->second )
                {
                    auto markList = this->markerNeumannBC( NeumannBCShape::VECTORIAL,myBcDesc.marker() );
                    if ( markList.empty() ) continue;
                    ExpressionStringAtMarker myBcDesc2( myBcDesc );
                    myBcDesc2.setMeshMarkers( markList );
                    bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc2 );
                }
            }
            auto itFindNeumannTensor2Type = itFindFieldVelocity->second.find("Neumann_tensor2");
            if ( itFindNeumannTensor2Type != itFindFieldVelocity->second.end() )
            {
                for ( auto const& myBcDesc : itFindNeumannTensor2Type->second )
                {
                    auto markList = this->markerNeumannBC( NeumannBCShape::TENSOR2,myBcDesc.marker() );
                    if ( markList.empty() ) continue;
                    ExpressionStringAtMarker myBcDesc2( myBcDesc );
                    myBcDesc2.setMeshMarkers( markList );
                    bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc2 );
                }
            }
        }
#if 0
        auto itFindFieldFluid = this->modelProperties().boundaryConditions().find("fluid");
        if ( itFindFieldFluid != this->modelProperties().boundaryConditions().end() )
        {
            auto itFindOutletType = itFindFieldFluid->second.find("outlet");
            if ( itFindOutletType != itFindFieldFluid->second.end() )
            {
                for ( auto const& myBcDesc : itFindOutletType->second )
                    bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc );
            }
        }
#else
        if ( !this->M_fluidOutletsBCType.empty() )
        {
            std::list<std::string> markList;
            for ( auto const& bcOutlet : this->M_fluidOutletsBCType )
                markList.push_back( std::get<0>(bcOutlet) );
            ExpressionStringAtMarker myBcDesc2( std::make_tuple( "expression","wind","0","","" ) );
            myBcDesc2.setMeshMarkers( markList );
            bcPrecPCD["velocity"]["Neumann"].push_back( myBcDesc2 );
        }
#endif

        // TODO other bc (fsi,...)
#if 1
        if ( this->worldComm().isMasterRank() && this->verbose() )
        {
            for( auto const& s : bcPrecPCD )
            {
                std::cout << "field " << s.first << "\n";
                for( auto const& t : s.second )
                {
                    std::cout << " - type " << t.first << "\n";
                    for( auto const& c : t.second )
                    {
                        std::ostringstream ostrMarkers;
                        ostrMarkers << "(";
                        for ( std::string const& mark : c.meshMarkers() )
                            ostrMarkers << mark << " ";
                        ostrMarkers << ")";
                        if ( c.hasExpression2() )
                            std::cout << "  . boundary  " << c.marker() << " " << ostrMarkers.str() << " expr : " << c.expression1() << " expr2:" << c.expression2() << "\n";
                        else
                            std::cout << "  . boundary  " << c.marker() << " " << ostrMarkers.str() << " expr : " << c.expression() << "\n";
                    }
                }
            }
        }
#endif
        //CHECK( this->algebraicFactory()->preconditionerTool()->matrix() ) << "no matrix define in preconditionerTool";

        // build pcd operator
        boost::shared_ptr<OperatorPCD<space_fluid_type>> opPCD;
        opPCD = boost::make_shared<OperatorPCD<space_fluid_type>>( this->functionSpace(),this->backend(),bcPrecPCD,"velocity",false,true);
        this->algebraicFactory()->preconditionerTool()->attachOperatorPCD("pcd",opPCD);
    }

}

} // namespace FeelModels
} // namespace Feel
