/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmecbase.hpp>

#include <feel/feelfilters/savegmshmesh.hpp>
//#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
//#include <feel/feelvf/inv.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>
#include <feel/feelmodels/modelmesh/markedmeshtool.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>

namespace Feel {
namespace FeelModels {

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::FluidMechanicsBase( std::string const& prefix,
                                                            bool buildMesh,
                                                            WorldComm const& worldComm,
                                                            std::string const& subPrefix,
                                                            std::string const& rootRepository )
    :
    super_type( prefix,worldComm,subPrefix, self_type::expandStringFromSpec( rootRepository ) ),
    M_hasBuildFromMesh( false ), M_isUpdatedForUse(false ),
    M_densityViscosityModel( new densityviscosity_model_type( prefix ) )
{
    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".FluidMechanicsConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".FluidMechanicsSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".FluidMechanicsPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".FluidMechanicsTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
std::string
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::expandStringFromSpec( std::string const& expr )
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
FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
const uint16_type FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::nOrderVelocity;
FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
const uint16_type FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::nOrderPressure;
FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
const uint16_type FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::nOrderGeo;


//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::build()
{
    this->log("FluidMechanics","build", "start");
    //-----------------------------------------------------------------------------//
    // create or reload mesh
    this->createMesh();
    //-----------------------------------------------------------------------------//
    // functionSpaces and elements
    this->createFunctionSpaces();
    //-----------------------------------------------------------------------------//
    // bdf time schema
    this->createTimeDiscretisation();
    //-----------------------------------------------------------------------------//
    // ALE mode (maybe)
    this->createALE();
    //-----------------------------------------------------------------------------//
    // physical parameters
    this->createOthers();
    //-----------------------------------------------------------------------------//
    //export
    this->createPostProcess();
    //-----------------------------------------------------------------------------//
    M_hasBuildFromMesh = true;
    this->log("FluidMechanics","build", "finish");
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::loadMesh( mesh_ptrtype __mesh )
{
    this->log("FluidMechanics","loadMesh", "start");
    //-----------------------------------------------------------------------------//
    // create or reload mesh
    if (this->doRestart() && !__mesh) this->createMesh();
    else M_mesh = __mesh;
    //-----------------------------------------------------------------------------//
    // functionSpaces and elements
    this->createFunctionSpaces();
    //-----------------------------------------------------------------------------//
    // bdf time schema
    this->createTimeDiscretisation();
    //-----------------------------------------------------------------------------//
    // ALE mode (maybe)
    this->createALE();
    //-----------------------------------------------------------------------------//
    // physical parameters
    this->createOthers();
    //-----------------------------------------------------------------------------//
    //export
    this->createPostProcess();
    //-----------------------------------------------------------------------------//
    M_hasBuildFromMesh = true;
    this->log("FluidMechanics","loadMesh", "finish");
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    this->log("FluidMechanics","loadParameterFromOptionsVm", "start");

    //--------------------------------------------------------------//
    // exporters options
    M_isHOVisu = nOrderGeo > 1;
    if ( Environment::vm().count(prefixvm(this->prefix(),"hovisu").c_str()) )
        M_isHOVisu = boption(_name="hovisu",_prefix=this->prefix());

    // overwrite export field options in json if given in cfg
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_velocity").c_str()) )
        if ( boption(_name="do_export_velocity",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Velocity );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_pressure").c_str()) )
        if ( boption(_name="do_export_pressure",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Pressure );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_displacement").c_str()) )
        if ( boption(_name="do_export_displacement",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Displacement );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_vorticity").c_str()) )
        if ( boption(_name="do_export_vorticity",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Vorticity );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_normalstress").c_str()) )
        if ( boption(_name="do_export_normalstress",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::NormalStress );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_wallshearstress").c_str()) )
        if ( boption(_name="do_export_wallshearstress",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::WallShearStress );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_density").c_str()) )
        if ( boption(_name="do_export_density",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Density );
    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_viscosity").c_str()) )
        if ( boption(_name="do_export_viscosity",_prefix=this->prefix()) )
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Viscosity );

    if ( boption(_name="do_export_meshale",_prefix=this->prefix()) )
        this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::ALEMesh );

    if ( Environment::vm().count(prefixvm(this->prefix(),"do_export_all").c_str()) )
        if ( boption(_name="do_export_all",_prefix=this->prefix()) )
        {
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Velocity );
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Pressure );
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Displacement );
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Vorticity );
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::NormalStress );
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::WallShearStress );
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Density );
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Viscosity );
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::ALEMesh );
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::Pid );
            this->M_postProcessFieldExported.insert( FluidMechanicsPostProcessFieldExported::LagrangeMultiplierPressureBC );
        }

    //--------------------------------------------------------------//
    M_haveSourceAdded=false;//true when update
    M_velocityDivIsEqualToZero=true;

    M_dirichletBCnitscheGamma = doption(_name="dirichletbc.nitsche.gamma",_prefix=this->prefix());

    std::string theFluidModel = this->modelProperties().model();
    if ( Environment::vm().count(prefixvm(this->prefix(),"model").c_str()) )
        theFluidModel = soption(_name="model",_prefix=this->prefix());
    this->setModelName( theFluidModel );

    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        this->setSolverName( soption(_name="solver",_prefix=this->prefix()) );

    //--------------------------------------------------------------//
    // fsi options
    M_useFSISemiImplicitScheme = false;
    M_couplingFSIcondition = "dirichlet-neumann";
    M_couplingFSI_Nitsche_gamma = 2500;
    M_couplingFSI_Nitsche_gamma0 = 1;
    M_couplingFSI_Nitsche_alpha = 1;
    M_couplingFSI_RNG_useInterfaceOperator = false;
    M_couplingFSI_solidIs1dReduced=false;

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

    //--------------------------------------------------------------//
    // gravity
    std::string gravityStr;
    if ( Environment::vm().count(prefixvm(this->prefix(),"gravity-force").c_str()) )
        gravityStr = soption(_name="gravity-force",_prefix=this->prefix());
    else if (nDim == 2 )
        gravityStr = "{0,-9.80665}";
    else if (nDim == 3 )
        gravityStr = "{0,0,-9.80665}";
    M_gravityForce = expr<nDim,1,2>( gravityStr );
    M_useGravityForce = boption(_name="use-gravity-force",_prefix=this->prefix());

    // thermodynamics coupling
    M_useThermodynModel = boption(_name="use-thermodyn",_prefix=this->prefix());
    M_BoussinesqRefTemperature = doption(_name="Boussinesq.ref-temperature",_prefix=this->prefix());

    this->log("FluidMechanics","loadParameterFromOptionsVm", "finish");
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createWorldsComm()
{
    this->log("FluidMechanics","createWorldsComm", "start");

    if (this->worldComm().localSize()==this->worldComm().globalSize())
    {
        std::vector<WorldComm> vecWorldComm(space_fluid_type::nSpaces,this->worldComm());
        std::vector<WorldComm> vecLocalWorldComm(1,this->worldComm());
        this->setWorldsComm(vecWorldComm);
        this->setLocalNonCompositeWorldsComm(vecLocalWorldComm);
    }
    else
    {
        // manage world comm : WARNING only true without the lagrange multiplier
        const int VelocityWorld=0;
        const int PressureWorld=1;
        const int LagrangeWorld=2;
        int CurrentWorld=0;
        if (this->worldComm().globalRank() < this->worldComm().globalSize()/2 )
            CurrentWorld=VelocityWorld;
        else
            CurrentWorld=PressureWorld;

        std::vector<WorldComm> vecWorldComm(space_fluid_type::nSpaces,this->worldComm());
        std::vector<WorldComm> vecLocalWorldComm(1,this->worldComm());
        if (this->worldComm().globalSize()>1)
        {
            vecWorldComm[0]=this->worldComm().subWorldComm(VelocityWorld);
            vecWorldComm[1]=this->worldComm().subWorldComm(PressureWorld);

            vecLocalWorldComm[0]=this->worldComm().subWorldComm(CurrentWorld);
        }
        else
        {
            vecWorldComm[0]=WorldComm();
            vecWorldComm[1]=WorldComm();
            //vecWorldCommVel[0]=WorldComm();
            vecLocalWorldComm[0]=WorldComm();
        }
        this->setWorldsComm(vecWorldComm);
        this->setLocalNonCompositeWorldsComm(vecLocalWorldComm);
    }
    this->log("FluidMechanics","createWorldsComm", "finish");
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createMesh()
{
    this->log("FluidMechanics","createMesh", "start");
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    double tElapsed = this->timerTool("Constructor").stop("createMesh");
    this->log("FluidMechanics","createMesh", (boost::format("finish in %1% s") %tElapsed).str() );
}

//---------------------------------------------------------------------------------------------------------//
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
FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createFunctionSpaces()
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
    M_densityViscosityModel->updateForUse( this->mesh(), this->modelProperties().materials(),  this->localNonCompositeWorldsComm(), hasExtendedDofTable );

    // fluid mix space : velocity and pressure
    if ( M_densityViscosityModel->isDefinedOnWholeMesh() )
    {
        M_rangeMeshElements = elements(this->mesh());
        M_Xh = space_fluid_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm(),
                                      _extended_doftable=extendedDT );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), M_densityViscosityModel->markers());
        M_Xh = space_fluid_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm(),
                                      _extended_doftable=extendedDT, _range=M_rangeMeshElements );
    }

    
    
#if 0
#if 1
    M_Xh = detail::createFluidFunctionSpaces(*this,extendedDT,mpl::bool_<UsePeriodicity>());
#else
#warning TODOhere!!!
#if !FLUIDMECHANICS_USE_PERIODICITY
    M_Xh = space_fluid_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm(),
                                  _extended_doftable=extendedDT );
#else
    node_type translat( nDim );
    translat[0] = doption(_name="periodicity.translate-x",_prefix=this->prefix());
    if ( nDim >=2 )
        translat[1] = doption(_name="periodicity.translate-y",_prefix=this->prefix());
    if ( nDim == 3 )
        translat[2]= doption(_name="periodicity.translate-z",_prefix=this->prefix());
    std::string marker1 = soption(_name="periodicity.marker1",_prefix=this->prefix());
    std::string marker2 = soption(_name="periodicity.marker2",_prefix=this->prefix());
    auto theperiodicity = periodicity( Periodic<>( this->mesh()->markerName(marker1),this->mesh()->markerName(marker2), translat), NoPeriodicity() );
    M_Xh = space_fluid_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm(),
                                  _extended_doftable=extendedDT,
                                  _periodicity=theperiodicity );
#endif
#endif
#endif
    M_Solution.reset( new element_fluid_type(M_Xh,"U"));

    // space usefull to tranfert sigma*N()
    //this->createFunctionSpacesNormalStress();

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh->worldComm() );

    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        if ( M_densityViscosityModel->isDefinedOnWholeMesh() )
            M_XhMeanPressureLM = space_meanpressurelm_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
        else
            M_XhMeanPressureLM = space_meanpressurelm_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm(),
                                                                 _range=M_rangeMeshElements );
    }


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

    double tElapsed = this->timerTool("Constructor").stop("createSpaces");
    this->log("FluidMechanics","createFunctionSpaces", (boost::format("finish in %1% s") %tElapsed).str() );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createTimeDiscretisation()
{
    this->log("FluidMechanics","createTimeDiscretisation", "start" );
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
                               fs::path( prefixvm(this->prefix(), (boost::format("bdf_o_%1%_dt_%2%")%this->timeStep() %M_bdf_fluid->bdfOrder()).str() ) ) ).string() );

    double tElapsed = this->timerTool("Constructor").stop("createTimeDiscr");
    this->log("FluidMechanics","createTimeDiscretisation", (boost::format("finish in %1% s") %tElapsed).str() );
}

//---------------------------------------------------------------------------------------------------------//


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createALE()
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
                                           this->rootRepositoryWithoutNumProc() ));
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

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createPostProcess()
{
    this->log("FluidMechanics","createPostProcess", "start" );
    this->timerTool("Constructor").start();
    this->createPostProcessExporters();
    //this->createPostProcessMeasures();
    double tElapsed = this->timerTool("Constructor").stop("createPostProcess");
    this->log("FluidMechanics","createPostProcess", (boost::format("finish in %1% s") %tElapsed).str() );
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createPostProcessExporters()
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


#if 0
        M_exporter_gmsh = gmsh_export_type::New( this->application()->vm(),
                                                 prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Export_gmshHo")), M_Xh->worldComm() );
#endif
    }

    if ( M_isHOVisu && doExport )
    {
#if 1 //defined(FEELPP_HAS_VTK)
        //M_exporter_ho = export_ho_type::New( this->application()->vm(), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Export_HO"))/*.c_str()*/, M_Xh->worldComm() );

#if defined( FEELPP_MODELS_HAS_MESHALE )
        if (M_isMoveDomain) this->meshALE()->revertReferenceMesh();
#endif
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
        if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::NormalStress ) ||
             this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::WallShearStress ) )
            M_XhVectorialDiscVisuHO = space_vectorialdisc_visu_ho_type::New(_mesh=meshVisuHO/*opLagP1->mesh()*/,_worldscomm=this->localNonCompositeWorldsComm());

        M_velocityVisuHO.reset( new element_vectorial_visu_ho_type(M_XhVectorialVisuHO,"u_visuHO"));
        M_pressureVisuHO.reset( new element_scalar_visu_ho_type(M_XhScalarVisuHO,"p_visuHO"));
        if (M_isMoveDomain) M_meshdispVisuHO.reset( new element_vectorial_visu_ho_type(M_XhVectorialVisuHO,"meshdisp_visuHO"));
        if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::NormalStress ) )
            M_normalStressVisuHO.reset( new element_vectorialdisc_visu_ho_type(M_XhVectorialDiscVisuHO,"normalstress_visuHO") );
        if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::WallShearStress ) )
            M_fieldWallShearStressVisuHO.reset( new element_vectorialdisc_visu_ho_type(M_XhVectorialDiscVisuHO,"wallshearstress_visuHO") );

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

        if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::NormalStress ) ||
             this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::WallShearStress ) )
        {
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

#if defined( FEELPP_MODELS_HAS_MESHALE )
        if (M_isMoveDomain) this->meshALE()->revertMovingMesh();
#endif

#endif
    }

    this->log("FluidMechanics","createPostProcessExporters", "finish" );

} // createPostProcessExporters

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createOthers()
{
    this->log("FluidMechanics","createOthers", "start" );
    this->timerTool("Constructor").start();

    //----------------------------------------------------------------------------//
    // space usefull to tranfert sigma*N()
    if (this->isMoveDomain()) this->createFunctionSpacesNormalStress();
    //----------------------------------------------------------------------------//
    // fluid inlet
    this->createBCFluidInlet();
    //----------------------------------------------------------------------------//
    // fluid outlet
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
        if ( M_isMoveDomain )
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            M_fluidOutletWindkesselSpaceMeshDisp = space_fluidoutlet_windkessel_mesh_disp_type::New( _mesh=M_fluidOutletWindkesselMesh,
                                                                                                     _worldscomm=this->localNonCompositeWorldsComm() );
            M_fluidOutletWindkesselMeshDisp = M_fluidOutletWindkesselSpaceMeshDisp->elementPtr();
            M_fluidOutletWindkesselOpMeshDisp = opInterpolation(_domainSpace=M_meshALE->functionSpace(),
                                                                _imageSpace=M_fluidOutletWindkesselSpaceMeshDisp,
                                                                _range=elements(M_fluidOutletWindkesselMesh),
                                                                _backend=M_backend );
#endif
        }
    }
    this->timerTool("Constructor").stop("createOthers");
    this->log("FluidMechanics","createOthers", "finish" );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createFunctionSpacesNormalStress()
{
    if ( M_XhNormalBoundaryStress ) return;

    if ( M_densityViscosityModel->isDefinedOnWholeMesh() )
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

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createFunctionSpacesVorticity()
{
    if ( M_densityViscosityModel->isDefinedOnWholeMesh() )
        M_XhVorticity = space_vorticity_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm());
    else
        M_XhVorticity = space_vorticity_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm(),
                                                   _range=M_rangeMeshElements );
    M_fieldVorticity.reset( new element_vorticity_type(M_XhVorticity));
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createFunctionSpacesSourceAdded()
{
    if ( M_densityViscosityModel->isDefinedOnWholeMesh() )
        M_XhSourceAdded=space_vectorial_PN_type::New( _mesh=M_mesh,_worldscomm=this->localNonCompositeWorldsComm() );
    else
        M_XhSourceAdded=space_vectorial_PN_type::New( _mesh=M_mesh,_worldscomm=this->localNonCompositeWorldsComm(),
                                                      _range=M_rangeMeshElements );
    M_SourceAdded.reset( new element_vectorial_PN_type(M_XhSourceAdded,"SourceAdded"));
}

//---------------------------------------------------------------------------------------------------------//


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createBCFluidInlet()
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

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateFluidInletVelocity()
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


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::init( bool buildMethodNum,
                                              typename model_algebraic_factory_type::appli_ptrtype const& app )
{
    if ( M_isUpdatedForUse ) return;

    this->log("FluidMechanics","init", "start" );
    this->timerTool("Constructor").start();

    boost::timer thetimer;

    if ( !M_hasBuildFromMesh )
        this->build();

    // update definePressureCst respect to the method choosen
    if ( this->definePressureCst() )
        this->updateDefinePressureCst();

    // lagrange multiplier for pressure bc
    if ( this->hasMarkerPressureBC() )
    {
        M_meshLagrangeMultiplierPressureBC = createSubmesh(this->mesh(),markedfaces(this->mesh(),this->markerPressureBC()) );
        M_spaceLagrangeMultiplierPressureBC = space_trace_velocity_component_type::New( _mesh=M_meshLagrangeMultiplierPressureBC, _worldscomm=this->localNonCompositeWorldsComm() );
        M_fieldLagrangeMultiplierPressureBC1.reset( new element_trace_velocity_component_type( M_spaceLagrangeMultiplierPressureBC ) );
        if ( nDim == 3 )
            M_fieldLagrangeMultiplierPressureBC2.reset( new element_trace_velocity_component_type( M_spaceLagrangeMultiplierPressureBC ) );
    }

    // update marker in mesh (mainly used with CIP stab)
    if ( (this->doCIPStabConvection() || this->doCIPStabDivergence() || this->doCIPStabPressure() ) && !this->applyCIPStabOnlyOnBoundaryFaces() )
        this->updateMarkedZonesInMesh();

    if ( M_useThermodynModel )
    {
        M_thermodynModel.reset( new thermodyn_model_type(prefixvm(this->prefix(),"thermo"), false, this->worldComm(),
                                                         this->subPrefix(), this->rootRepositoryWithoutNumProc() ) );
        M_thermodynModel->setFieldVelocityConvectionIsUsed( !M_useGravityForce/*false*/ );
        M_thermodynModel->loadMesh( this->mesh() );
        // disable thermo exporter if we use fluid exporter
        M_thermodynModel->setDoExportResults( false );
        M_thermodynModel->init( !M_useGravityForce/*false*/ );

        M_rangeMeshElementsAeroThermal = intersect( M_rangeMeshElements, M_thermodynModel->rangeMeshElements() );
    }
    //-------------------------------------------------//
    // add ALE markers
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
        if (this->doRestart())
        {
            if ( this->hasFluidOutletWindkesselImplicit() )
            {
                // interpolate disp
                M_fluidOutletWindkesselOpMeshDisp->apply( *M_meshALE->displacement(), *M_fluidOutletWindkesselMeshDisp );
                // apply disp
                M_fluidOutletWindkesselMeshMover.apply( M_fluidOutletWindkesselMesh, *M_fluidOutletWindkesselMeshDisp );
            }
        }

        this->log("FluidMechanics","init", "meshALE done" );
#endif
    }

    //-------------------------------------------------//
    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();
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
    this->initFluidOutlet();
    // init function defined in json
    this->initUserFunctions();
    // init post-processinig (exporter, measure at point, ...)
    this->initPostProcess();
    //-------------------------------------------------//
    // define start dof index ( lm , windkessel )
    size_type currentStartIndex = 2;// velocity and pressure before
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        M_startBlockIndexFieldsInMatrix["define-pressure-cst-lm"] = currentStartIndex++;
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
    if ( M_useThermodynModel && M_useGravityForce )
    {
        M_thermodynModel->setRowStartInMatrix( currentStartIndex );
        M_thermodynModel->setColStartInMatrix( currentStartIndex );
        M_thermodynModel->setRowStartInVector( currentStartIndex );
        ++currentStartIndex;
    }
    //-------------------------------------------------//
    // prepare block vector
    int nBlock = this->nBlockMatrixGraph();
    M_blockVectorSolution.resize( nBlock );
    M_blockVectorSolution(0) = this->fieldVelocityPressurePtr();
    int cptBlock=1;
    // impose mean pressure by lagrange multiplier
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        M_blockVectorSolution(cptBlock) = this->backend()->newVector( M_XhMeanPressureLM );
        ++cptBlock;
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
    // thermodynamics model
    if ( M_useThermodynModel && M_useGravityForce )
    {
        M_blockVectorSolution(cptBlock++) = M_thermodynModel->fieldTemperaturePtr();
    }

    // init vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );
    //-------------------------------------------------//
    if (buildMethodNum)
    {
        M_algebraicFactory.reset( new model_algebraic_factory_type(app,this->backend()) );
#if 1
        bool attachMassMatrix = boption(_prefix=this->prefix(),_name="preconditioner.attach-mass-matrix");
        if ( attachMassMatrix )
        {
            auto massbf = form2( _trial=this->functionSpaceVelocity(), _test=this->functionSpaceVelocity());
            auto const& u = this->fieldVelocity();
            double coeff = this->densityViscosityModel()->cstRho()*this->timeStepBDF()->polyDerivCoefficient(0);
            if ( this->isStationary() ) coeff=1.;
            massbf += integrate( _range=elements( this->mesh() ), _expr=coeff*inner( idt(u),id(u) ) );
            massbf.matrixPtr()->close();
            M_algebraicFactory->preconditionerTool()->attachAuxiliarySparseMatrix( "mass-matrix", massbf.matrixPtr() );
        }
#endif
    }
    //-------------------------------------------------//
    this->updateBoundaryConditionsForUse();
    //-------------------------------------------------//
    M_isUpdatedForUse = true;

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("FluidMechanics","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateMarkedZonesInMesh()
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


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    // start or restart time step scheme
    if (!this->doRestart())
    {
        if ( ( !this->startBySolveStokesStationary() ) ||
             ( this->startBySolveStokesStationary() && this->hasSolveStokesStationaryAtKickOff() ) )
        {
            // start time step
            M_bdf_fluid->start(*M_Solution);
            // up current time
            this->updateTime( M_bdf_fluid->time() );
        }
    }
    else
    {
        // start time step
        M_bdf_fluid->restart();
        // load a previous solution as current solution
        *M_Solution = M_bdf_fluid->unknown(0);
        // up initial time
        this->setTimeInitial( M_bdf_fluid->timeInitial() );
        // restart exporter
        //this->restartPostProcess();
        // up current time
        this->updateTime( M_bdf_fluid->time() );

        this->log("FluidMechanics","initTimeStep", "restart bdf/exporter done" );
    }
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::initFluidOutlet()
{
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

}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::initUserFunctions()
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

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateUserFunctions( bool onlyExprWithTimeSymbol )
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


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    // update post-process expression
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    bool hasMeasure = false;

    // clean doExport with fields not available
    if ( !this->isMoveDomain() )
    {
        M_postProcessFieldExported.erase( FluidMechanicsPostProcessFieldExported::Displacement );
        M_postProcessFieldExported.erase( FluidMechanicsPostProcessFieldExported::ALEMesh );
    }

    // restart exporters if restart is activated
    if ( this->doRestart() && this->restartPath().empty() )
    {
        // if restart and same directory, update the exporter for new value, else nothing (create a new exporter)
        if ( M_exporter && M_exporter->doExport() )
            M_exporter->restart( this->timeInitial() );
        if ( M_exporter_ho && M_exporter_ho->doExport() )
            M_exporter_ho->restart( this->timeInitial() );
    }

    // add user functions
    if ( this->modelProperties().postProcess().find("Fields") != this->modelProperties().postProcess().end() )
    {
        for ( auto const& o : this->modelProperties().postProcess().find("Fields")->second )
        {
            if ( this->hasFieldUserScalar( o ) || this->hasFieldUserVectorial( o ) )
                M_postProcessUserFieldExported.insert( o );
        }
    }

    // forces (lift, drag) and flow rate measures
    auto const& ptree = this->modelProperties().postProcess().pTree();
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
                this->modelProperties().postProcess().operator[](ppTypeMeasures).push_back( "Pressure" );
                this->postProcessMeasuresIO().setMeasure("pressure_sum",0.);
                this->postProcessMeasuresIO().setMeasure("pressure_mean",0.);
                hasMeasure = true;
            }
            else if ( ptreeLevel1Name == "VelocityDivergence" )
            {
                this->modelProperties().postProcess().operator[](ppTypeMeasures).push_back( "VelocityDivergence" );
                this->postProcessMeasuresIO().setMeasure("velocity_divergence_sum",0.);
                this->postProcessMeasuresIO().setMeasure("velocity_divergence_mean",0.);
                this->postProcessMeasuresIO().setMeasure("velocity_divergence_normL2",0.);
                hasMeasure = true;
            }
        }
    }

    // point measures
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint() )
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




} // namespace FeelModels
} // namespace Feel
