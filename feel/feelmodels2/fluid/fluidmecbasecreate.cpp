/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels2/fluid/fluidmecbase.hpp>

#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
//#include <feel/feelvf/inv.hpp>

#include <feel/feelmodels2/modelmesh/reloadmesh.hpp>
#include <feel/feelmodels2/modelmesh/markedmeshtool.hpp>

namespace Feel {
namespace FeelModels {

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::FluidMechanicsBase( //bool __isStationary,
                                                            std::string __prefix,
                                                            bool __buildMesh,
                                                            WorldComm const& __worldComm,
                                                            std::string __subPrefix,
                                                            std::string __appliShortRepository )
    :
    super_type( __prefix,__worldComm,__subPrefix,__appliShortRepository),
    M_hasBuildFromMesh( false ), M_isUpdatedForUse(false ),
    M_densityViscosityModel( new densityviscosity_model_type(  __prefix ) )
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
    this->createExporters();
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
    if (this->doRestart()) this->createMesh();
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
    this->createExporters();
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

    M_meshSize = doption(_name="hsize",_prefix=this->prefix());
#if defined(FEELPP_HAS_VTK)
    M_isHOVisu =nOrderGeo > 1;
    if ( Environment::vm().count(prefixvm(this->prefix(),"hovisu").c_str()) )
        M_isHOVisu = boption(_name="hovisu",_prefix=this->prefix());
#else
    M_isHOVisu=false;
    if ( Environment::vm().count(prefixvm(this->prefix(),"hovisu").c_str()) )
        FeelModels::Log(this->prefix()+".FluidMechanics","constructor", "WARNING : hovisu disable because VTK not find",
                        this->worldComm(),this->verboseAllProc());
#endif
    M_doExportMeshALE = boption(_name="do_export_meshale",_prefix=this->prefix());
    M_doExportVorticity = boption(_name="do_export_vorticity",_prefix=this->prefix());
    M_doExportNormalStress = boption(_name="do_export_normalstress",_prefix=this->prefix());
    M_doExportWallShearStress = boption(_name="do_export_wallshearstress",_prefix=this->prefix());
    M_doExportMeshDisplacementOnInterface = boption(_name="do_export_meshdisplacementoninterface",_prefix=this->prefix());
    M_doExportViscosity = boption(_name="do_export_viscosity",_prefix=this->prefix());
    M_doExportAll = boption(_name="do_export_all",_prefix=this->prefix());

    M_haveSourceAdded=false;//true when update
    M_velocityDivIsEqualToZero=true;

    M_dirichletBCnitscheGamma = doption(_name="dirichletbc.nitsche.gamma",_prefix=this->prefix());

    this->pdeType( soption(_name="model",_prefix=this->prefix()) ); // up M_pdeType
    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        M_pdeSolver = soption(_name="solver",_prefix=this->prefix());
    //M_stressTensorLaw = soption(_name="stress_tensor_law",_prefix=this->prefix());
    M_useFSISemiImplicitScheme = false;
    M_couplingFSIcondition = "dirichlet";
    M_gammaNitschFSI = 2500;
    M_gamma0NitschFSI = 1;

    M_startBySolveNewtonian = boption(_prefix=this->prefix(),_name="start-by-solve-newtonian");
    M_hasSolveNewtonianAtKickOff = false;
    M_startBySolveStokesStationary = boption(_prefix=this->prefix(),_name="start-by-solve-stokes-stationary");
    M_hasSolveStokesStationaryAtKickOff = false;

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
    M_definePressureCstPenalisationBeta = doption(_name="define-pressure-cst.penalisation-beta",_prefix=this->prefix());

    // fluid outlet
    M_nFluidOutlet = this->hasFluidOutlet()? ioption(_name="fluid-outlet.number", _prefix=this->prefix()) : 0;
    if ( this->hasFluidOutlet() )
    {
        M_fluidOutletMarkerName.resize(M_nFluidOutlet);

        auto itMarkerFluidOutlet = M_fluidOutletsBCType[this->fluidOutletType()].begin();
        CHECK( itMarkerFluidOutlet != M_fluidOutletsBCType[this->fluidOutletType()].end() ) << "no fluid outlet bc found with " << this->fluidOutletType() << "\n";

        std::string markerNameBFOutletBase = *itMarkerFluidOutlet;
        bool useGenericMarkerName = (M_fluidOutletsBCType[this->fluidOutletType()].size() == 1) && (M_nFluidOutlet > 1);

        if ( !useGenericMarkerName )
            CHECK( M_fluidOutletsBCType[this->fluidOutletType()].size() == M_nFluidOutlet )  << "invalid nFluidOutlet : "
                                                                                             << M_fluidOutletsBCType[this->fluidOutletType()].size() << " and "
                                                                                             << M_nFluidOutlet << "\n";

        for (int k=0;k<M_nFluidOutlet;++k)
        {
            if ( useGenericMarkerName )
            {
                M_fluidOutletMarkerName[k] = markerNameBFOutletBase+(boost::format("%1%") %k).str();
            }
            else
            {
                M_fluidOutletMarkerName[k] = *itMarkerFluidOutlet;
                ++itMarkerFluidOutlet;
                //M_fluidOutletMarkerName[k] = markerNameBFOutletBase;
            }
            //M_meshAleBCType["free"].push_back( M_fluidOutletMarkerName[k] );
            this->addMarkerALEMeshBC("free",M_fluidOutletMarkerName[k]);
        }

        if ( this->fluidOutletType()=="windkessel" )
        {
            M_fluidOutletWindkesselCoupling = soption(_name="fluid-outlet.windkessel.coupling", _prefix=this->prefix());

            M_fluidOutletWindkesselPressureDistal.resize(M_nFluidOutlet,0);
            M_fluidOutletWindkesselPressureDistal_old.resize(M_nFluidOutlet,std::vector<double>(Feel::BDF_MAX_ORDER));
            M_fluidOutletWindkesselPressureProximal.resize(M_nFluidOutlet,0);
            //M_fluidOutletMarkerName.resize(M_nFluidOutlet);
            M_fluidOutletWindkesselRd.resize(M_nFluidOutlet);
            M_fluidOutletWindkesselRp.resize(M_nFluidOutlet);
            M_fluidOutletWindkesselCd.resize(M_nFluidOutlet);

            for (int k=0;k<M_nFluidOutlet;++k)
            {
                // windkessel parameter
                M_fluidOutletWindkesselRd[k] = doption(_name=(boost::format("fluid-outlet.windkessel.Rd%1%")%k).str(), _prefix=this->prefix());
                M_fluidOutletWindkesselRp[k] = doption(_name=(boost::format("fluid-outlet.windkessel.Rp%1%")%k).str(), _prefix=this->prefix());
                M_fluidOutletWindkesselCd[k] = doption(_name=(boost::format("fluid-outlet.windkessel.Cd%1%")%k).str(), _prefix=this->prefix());
            }
        } // windkessel
    } // hasFluidOutlet



    // physical parameters
    //M_CstRho = doption(_name="rho",_prefix=this->prefix());
    //M_CstMu = doption(_name="mu",_prefix=this->prefix());
    //M_CstNu = this->viscosityModel()->cstMu()/M_CstRho;

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

    // save path of file mesh
    //auto fmpath = this->fileNameMeshPath();//prefixvm(this->prefix(),"FluidMechanicsMesh.path");
    std::string fmpath = (fs::path( this->appliRepository() ) / fs::path(this->fileNameMeshPath())).string();
    if (this->doRestart())
    {
        this->log("FluidMechanics","createMesh", "restart with : "+fmpath);

        if ( !this->restartPath().empty() )
        {
            fmpath = (fs::path( this->restartPath() ) / fs::path(this->fileNameMeshPath())).string();
        }
        M_mesh = reloadMesh<mesh_type>(fmpath,this->worldComm());
    }
    else
    {
        if (this->hasMshfileStr())
        {
            std::string path = this->appliRepository();
            std::string mshfileRebuildPartitions = path + "/" + this->prefix() + ".msh";

            this->log("FluidMechanics","createMesh", "load msh file : " + this->mshfileStr());

            M_mesh = loadGMSHMesh(_mesh=new mesh_type,
                                  _filename=this->mshfileStr(),
                                  _worldcomm=this->worldComm(),
                                  _rebuild_partitions=this->rebuildMeshPartitions(),
                                  _rebuild_partitions_filename=mshfileRebuildPartitions,
                                  _partitions=this->worldComm().localSize(),
                                  _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK);

            if (this->rebuildMeshPartitions()) this->setMshfileStr(mshfileRebuildPartitions);
        }
        else if (this->hasGeofileStr())
        {
            std::string path = this->appliRepository();
            std::string mshfile = path + "/" + this->prefix() + ".msh";
            this->setMshfileStr(mshfile);

            fs::path curPath=fs::current_path();
            bool hasChangedRep=false;
            if ( curPath != fs::path(this->appliRepository()) )
            {
                this->log("FluidMechanics","createMesh", "change repository (temporary) for build mesh from geo : "+ this->appliRepository() );
                bool hasChangedRep=true;
                Environment::changeRepository( _directory=boost::format(this->appliRepository()), _subdir=false );
            }

            M_mesh = GeoTool::createMeshFromGeoFile<mesh_type>(this->geofileStr(),this->prefix(),M_meshSize,1,
                                                               this->worldComm().localSize(),this->worldComm());
            // go back to previous repository
            if ( hasChangedRep )
                Environment::changeRepository( _directory=boost::format(curPath.string()), _subdir=false );
        }
        else
        {
            std::string geotoolSavePath;
            if ( this->geotoolSaveDirectory()!=this->appliShortRepository() )
            {
                this->log("FluidMechanics","createMesh", "change rep -> "+ this->geotoolSaveDirectory() );
                Environment::changeRepository( _directory=boost::format(this->geotoolSaveDirectory()), _subdir=false );
                geotoolSavePath = Environment::rootRepository()+"/"+ this->geotoolSaveDirectory();
            }
            else
            {
                geotoolSavePath = this->appliRepository();
            }

            std::string geotoolSaveName = this->geotoolSaveName();
            std::string mshfile = geotoolSavePath + "/" + geotoolSaveName + ".msh";
            std::string geofilename = geotoolSavePath + "/" + geotoolSaveName;// without .geo
            this->setMshfileStr(mshfile);

            this->loadConfigMeshFile(geofilename);

            if ( this->geotoolSaveDirectory()!=this->appliShortRepository() )
            {
                this->log("FluidMechanics","createMesh", "change rep -> " + this->appliRepository() );
                Environment::changeRepository( _directory=boost::format(this->appliShortRepository()), _subdir=true );
            }

        }
        this->saveMSHfilePath(fmpath);
    }


    double timeElapsedCreateMesh = this->timerTool("Constructor").stop("createMesh");
    this->log("FluidMechanics","createMesh", (boost::format("finish in %1% s") % timeElapsedCreateMesh).str() );
} // createMesh()

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
    if ( (this->doCIPStabConvection() || this->doCIPStabDivergence()) && !this->applyCIPStabOnlyOnBoundaryFaces() )
    {
        this->log("FluidMechanics","createFunctionSpaces", "use buildDofTableMPIExtended on velocity" );
        extendedDT[0] = true;
    }
    if ( this->doCIPStabPressure() )
    {
        this->log("FluidMechanics","createFunctionSpaces", "use buildDofTableMPIExtended on pressure" );
        extendedDT[1] = true;
    }
    // fluid mix space : velocity and pressure
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
    M_Solution.reset( new element_fluid_type(M_Xh,"U"));

    // space usefull to tranfert sigma*N()
    //this->createFunctionSpacesNormalStress();

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh->worldComm() );

    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
        M_XhMeanPressureLM = space_meanpressurelm_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );


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

    this->timerTool("Constructor").stop("createSpaces");
    this->log("FluidMechanics","createFunctionSpaces", "finish" );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createTimeDiscretisation()
{
    this->log("FluidMechanics","createTimeDiscretisation", "start" );
    this->timerTool("Constructor").start();

    // bdf time schema
    std::string suffixName = "";
    if ( soption(_name="ts.file-format",_prefix=this->prefix()) == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_bdf_fluid = bdf( _vm=Environment::vm(), _space=M_Xh,
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

    M_bdf_fluid->setPathSave( (fs::path(this->appliRepository()) /
                               fs::path( prefixvm(this->prefix(), (boost::format("bdf_o_%1%_dt_%2%")%this->timeStep() %M_bdf_fluid->bdfOrder()).str() ) ) ).string() );

    this->timerTool("Constructor").stop("createTimeDiscr");
    this->log("FluidMechanics","createTimeDiscretisation", "finish" );
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
                                           this->appliShortRepository() ));
        // mesh displacement only on moving
        M_meshDisplacementOnInterface.reset( new element_mesh_disp_type(M_meshALE->displacement()->functionSpace(),"mesh_disp_on_interface") );
        // mesh velocity only on moving interface
        M_meshVelocityInterface.reset(new element_meshvelocityonboundary_type( M_XhMeshVelocityInterface, "mesh_velocity_interface" ) );

        this->timerTool("Constructor").stop("createALE");
        this->log("FluidMechanics","createALE", "finish");
    }
#endif

}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createExporters()
{
    this->log("FluidMechanics","createExporters", "start" );
    this->timerTool("Constructor").start();

    //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;//(this->isMoveDomain())?ExporterGeometry::EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY:ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
    std::string geoExportType="static";//change_coords_only, change, static

#if 0
    if ( !fs::exists(this->exporterPath()) )
    {
        // master rank create directories
        if ( this->worldComm().isMasterRank() )
            fs::create_directories( this->exporterPath() );
        // wait for all process
        this->worldComm().globalComm().barrier();
    }
#endif

    if (!M_isHOVisu)
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
    else
    {
#if defined(FEELPP_HAS_VTK)
        //M_exporter_ho = export_ho_type::New( this->application()->vm(), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"Export_HO"))/*.c_str()*/, M_Xh->worldComm() );

#if defined( FEELPP_MODELS_HAS_MESHALE )
        if (M_isMoveDomain) this->getMeshALE()->revertReferenceMesh();
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
                                       _path=this->appliRepository(),
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
                                       _path=this->appliRepository(),
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
        if (M_doExportNormalStress || M_doExportWallShearStress )
            M_XhVectorialDiscVisuHO = space_vectorialdisc_visu_ho_type::New(_mesh=meshVisuHO/*opLagP1->mesh()*/,_worldscomm=this->localNonCompositeWorldsComm());

        M_velocityVisuHO.reset( new element_vectorial_visu_ho_type(M_XhVectorialVisuHO,"u_visuHO"));
        M_pressureVisuHO.reset( new element_scalar_visu_ho_type(M_XhScalarVisuHO,"p_visuHO"));
        if (M_isMoveDomain) M_meshdispVisuHO.reset( new element_vectorial_visu_ho_type(M_XhVectorialVisuHO,"meshdisp_visuHO"));
        if (M_doExportNormalStress) M_normalStressVisuHO.reset( new element_vectorialdisc_visu_ho_type(M_XhVectorialDiscVisuHO,"normalstress_visuHO") );
        if (M_doExportWallShearStress) M_wallShearStressVisuHO.reset( new element_vectorialdisc_visu_ho_type(M_XhVectorialDiscVisuHO,"wallshearstress_visuHO") );

        this->log("FluidMechanics","createExporters", "start opInterpolation" );
        boost::mpi::timer timerOpI;

        M_opIvelocity = opInterpolation(_domainSpace=M_Xh->template functionSpace<0>(),
                                        _imageSpace=M_XhVectorialVisuHO,
                                        _range=elements(M_XhVectorialVisuHO->mesh()),
                                        _backend=M_backend,
                                        _type=InterpolationNonConforme(false,true,false,15) );

        this->log("FluidMechanics","createExporters", "step1 done" );

        M_opIpressure = opInterpolation(_domainSpace=M_Xh->template functionSpace<1>(),
                                        _imageSpace=M_XhScalarVisuHO,
                                        _range=elements(M_XhScalarVisuHO->mesh()),
                                        _backend=M_backend,
                                        _type=InterpolationNonConforme(false,true,false,15) );

        if (M_doExportNormalStress || M_doExportWallShearStress )
        {
            M_opIstress = opInterpolation(_domainSpace=M_XhNormalBoundaryStress,
                                          _imageSpace=M_XhVectorialDiscVisuHO,
                                          _range=elements(M_XhVectorialDiscVisuHO->mesh()),
                                          _backend=M_backend,
                                          _type=InterpolationNonConforme(false,true,false,15) );
        }

        this->log("FluidMechanics","createExporters", "step2 done" );

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
        this->log("FluidMechanics","createExporters", "finish all opInterpolation in " + (boost::format("%1% s") % timeElapsedOpI).str() );

#if defined( FEELPP_MODELS_HAS_MESHALE )
        if (M_isMoveDomain) this->getMeshALE()->revertMovingMesh();
#endif

#endif
    }

    this->timerTool("Constructor").stop("createExporters");
    this->log("FluidMechanics","createExporters", "finish" );

} // createExporters

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createOthers()
{
    this->log("FluidMechanics","createOthers", "start" );
    this->timerTool("Constructor").start();
    //----------------------------------------------------------------------------//
    // rho, mu, nu with scalar P0 space
    //std::vector<bool> extendedDT(1,this->useExtendedDofTable() );
    //M_XhScalarP0 = space_densityviscosity_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm(),
    //                                                 _extended_doftable=extendedDT );
    //M_P0Rho.reset( new element_densityviscosity_type(M_XhScalarP0,"rho"));
    //M_P0Mu.reset( new element_densityviscosity_type(M_XhScalarP0,"mu"));
    //M_P0Nu.reset( new element_densityviscosity_type(M_XhScalarP0,"nu"));
    //*M_P0Rho= vf::project(_space=M_XhScalarP0, _range=elements( M_mesh),
    //                        _expr=vf::cst(M_CstRho),_geomap=this->geomap());
    //*M_P0Mu= vf::project(_space=M_XhScalarP0, _range=elements( M_mesh),
    //                     _expr=vf::cst(M_CstMu),_geomap=this->geomap());
    //*M_P0Nu= vf::project(_space=M_XhScalarP0, _range=elements( M_mesh),
    //                       _expr=vf::cst(M_CstNu),_geomap=this->geomap());
    // viscosity model
    //M_viscosityModelDesc.reset( new viscosity_model_type( this->stressTensorLawType(), *M_P0Mu, this->prefix() ) );
    //M_densityViscosityModel->initFromSpace(M_XhScalarP0);//,this->fieldVelocity(),this->fieldPressure());
    M_densityViscosityModel->initFromMesh( this->mesh(), this->useExtendedDofTable() );
    //----------------------------------------------------------------------------//
    // space usefull to tranfert sigma*N()
    if (this->isMoveDomain()) this->createFunctionSpacesNormalStress();
    //----------------------------------------------------------------------------//
    // fluid outlet
    if ( this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" && this->fluidOutletWindkesselCoupling() == "implicit" )
    {
        // list usefull to create the outlets submesh
        std::list<std::string> markerNameBFOutletForSubmesh;
        for (int k=0;k<M_nFluidOutlet;++k)
            markerNameBFOutletForSubmesh.push_back(M_fluidOutletMarkerName[k]);

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
    M_XhNormalBoundaryStress = space_stress_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    M_normalBoundaryStress.reset(new element_stress_type(M_XhNormalBoundaryStress));
#if defined( FEELPP_MODELS_HAS_MESHALE )
    M_normalStressFromStruct.reset(new element_stress_type(M_XhNormalBoundaryStress));
#endif
    M_wallShearStress.reset(new element_stress_type(M_XhNormalBoundaryStress));
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createFunctionSpacesVorticity()
{
    M_Xh_vorticity = space_vorticity_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm());
    M_vorticity.reset( new element_vorticity_type(M_Xh_vorticity));
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createFunctionSpacesSourceAdded()
{
    M_XhSourceAdded=space_vectorial_PN_type::New( _mesh=M_mesh );
    M_SourceAdded.reset( new element_vectorial_PN_type(M_XhSourceAdded,"SourceAdded"));
}

//---------------------------------------------------------------------------------------------------------//

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


    // build definePressureCst space if not done yet
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" && !M_XhMeanPressureLM )
        M_XhMeanPressureLM = space_meanpressurelm_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );

    // update marker in mesh (mainly used with CIP stab)
    if ( (this->doCIPStabConvection() || this->doCIPStabDivergence() || this->doCIPStabPressure() ) && !this->applyCIPStabOnlyOnBoundaryFaces() )
        this->updateMarkedZonesInMesh();

    //-------------------------------------------------//
    // add ALE markers
    if (M_isMoveDomain)
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
            if ( this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" && this->fluidOutletWindkesselCoupling() == "implicit" )
            {
                // interpolate disp
                M_fluidOutletWindkesselOpMeshDisp->apply( *M_meshALE->displacement(), *M_fluidOutletWindkesselMeshDisp );
                // apply disp
                M_fluidOutletWindkesselMeshMover.apply( /*this->*/M_fluidOutletWindkesselMesh, *M_fluidOutletWindkesselMeshDisp );
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
    // windkessel outlet
    if ( this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" )
    {
        for (int k=0;k<M_nFluidOutlet;++k)
        {
            //std::string nameFile = this->appliRepository() + "/" + prefixvm(this->prefix(),"bloodFlowOutlet.windkessel.data");
            //std::string nameFile =boost::format(this->appliRepository() + "/" + prefixvm(this->prefix(),"bloodFlowOutlet.windkessel%1%.data") %k ).str();
            std::string nameFile =this->appliRepository() + "/" + prefixvm(this->prefix(),(boost::format("bloodFlowOutlet.windkessel%1%.data") %k ).str());

            if (!this->doRestart())
            {
                M_fluidOutletWindkesselPressureDistal[k] = 0;
                M_fluidOutletWindkesselPressureProximal[k] = 0;
                for (int l=0;l<M_fluidOutletWindkesselPressureDistal_old[k].size();++l)
                    M_fluidOutletWindkesselPressureDistal_old[k][l] = 0;
                if (this->worldComm().isMasterRank())
                {
                    std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::trunc);
                    file << 0 << " " << this->timeInitial() << " " << M_fluidOutletWindkesselPressureDistal[k] << " " << M_fluidOutletWindkesselPressureProximal[k] << "\n";
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
                    while ( !fileI.eof() && !find )
                    {
                        fileI >> cptIter >> timeIter >> valPresDistal >> valPresProximal;
                        buffer << cptIter << " " << timeIter << " " << valPresDistal << " " << valPresProximal << "\n";

                        for (int l=0 ; l< M_fluidOutletWindkesselPressureDistal_old[k].size() ; ++l)
                            if (cptIter == askedIter - l) M_fluidOutletWindkesselPressureDistal_old[k][l] = valPresDistal;

                        if (cptIter == askedIter) find=true;

                        //if (cptIter == askedIter) { M_fluidOutletWindkesselPressureDistal_old[k][0]; find=true;}
                        //if (cptIter == askedIter-1) { M_fluidOutletWindkesselPressureDistal_old[k][1]; find=true;}
                    }
                    fileI.close();
                    std::ofstream fileW(nameFile.c_str(), std::ios::out | std::ios::trunc);
                    fileW << buffer.str();
                    fileW.close();
                    //std::cout << cptIter <<" " << timeIter << " " << valPresDistal << std::endl;
                    //M_fluidOutletWindkesselPressureDistal_old[k][0] = valPresDistal;
                }
            }
            mpi::broadcast( this->worldComm().globalComm(), M_fluidOutletWindkesselPressureDistal_old, this->worldComm().masterRank() );
        } // for (int k=0;k<M_nFluidOutlet;++k)

        this->log("FluidMechanics","init", "restart windkessel done" );

    } // if (M_hasFluidOutlet)

    //-------------------------------------------------//
    // define start dof index ( lm , windkessel )
    size_type currentStartIndex = 0;
    currentStartIndex += M_Xh->nLocalDofWithGhost();
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        M_startDofIndexFieldsInMatrix["define-pressure-cst-lm"] = currentStartIndex;
        currentStartIndex += M_XhMeanPressureLM->nLocalDofWithGhost() ;
    }
    if (this->hasMarkerDirichletBClm())
    {
        M_startDofIndexFieldsInMatrix["dirichletlm"] = currentStartIndex;
        currentStartIndex += this->XhDirichletLM()->nLocalDofWithGhost() ;
    }
    if (this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" && this->fluidOutletWindkesselCoupling() == "implicit" )
    {
        M_startDofIndexFieldsInMatrix["windkessel"] = currentStartIndex;
        currentStartIndex += 2*this->nFluidOutlet();
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
    // windkessel outel with implicit scheme
    if ( this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" &&
         this->fluidOutletWindkesselCoupling() == "implicit" )
    {
        for (int k=0;k<this->nFluidOutlet();++k)
        {
            M_blockVectorSolution(cptBlock) = this->backend()->newVector( M_fluidOutletWindkesselSpace );
            ++cptBlock;
        }
    }
    // init vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );
    //-------------------------------------------------//
    if (buildMethodNum)
    {
        M_algebraicFactory.reset( new model_algebraic_factory_type(app,this->backend()) );

    }
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
        this->restartExporters();
        // up current time
        this->updateTime( M_bdf_fluid->time() );

        this->log("FluidMechanics","initTimeStep", "restart bdf/exporter done" );
    }
}

//---------------------------------------------------------------------------------------------------------//

} // namespace FeelModels
} // namespace Feel
