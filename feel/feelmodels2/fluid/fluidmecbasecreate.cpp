/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels2/fluid/fluidmecbase.hpp>

#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
//#include <feel/feelvf/inv.hpp>

#include <feel/feelmodels2/modelmesh/reloadmesh.hpp>


namespace Feel {
namespace FeelModels {

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::FluidMechanicsBase( bool __isStationary,
                                                                                                                     std::string __prefix,
                                                                                                                     WorldComm const& __worldComm,
                                                                                                                     bool __buildMesh,
                                                                                                                     std::string __subPrefix,
                                                                                                                     std::string __appliShortRepository )
:
super_type( __isStationary,__prefix,__worldComm,__subPrefix,__appliShortRepository)
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

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::build()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","build", "start",
                                        this->worldComm(),this->verboseAllProc());

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

    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","build", "finish",
                                        this->worldComm(),this->verboseAllProc());

}

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::loadMesh( mesh_ptrtype __mesh )
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","loadMesh", "start",
                                        this->worldComm(),this->verboseAllProc());

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

    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","loadMesh", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::loadParameterFromOptionsVm()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","loadParameterFromOptionsVm", "start",
                                        this->worldComm(),this->verboseAllProc());

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
    M_stressTensorLaw = soption(_name="stress_tensor_law",_prefix=this->prefix());
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
    M_CstRho = doption(_name="rho",_prefix=this->prefix());
    M_CstMu = doption(_name="mu",_prefix=this->prefix());
    M_CstNu = M_CstMu/M_CstRho;

    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","loadParameterFromOptionsVm", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::createWorldsComm()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createWorldsComm", "start",
                                        this->worldComm(),this->verboseAllProc());

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
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createWorldsComm", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::createMesh()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createMesh", "start",
                                        this->worldComm(),this->verboseAllProc());
    this->timerTool("Constructor").start();

    // save path of file mesh
    //auto fmpath = this->fileNameMeshPath();//prefixvm(this->prefix(),"FluidMechanicsMesh.path");
    std::string fmpath = (fs::path( this->appliRepository() ) / fs::path(this->fileNameMeshPath())).string();
    if (this->doRestart())
    {
        if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createMesh",
                                            "restart with : "+fmpath,
                                            this->worldComm(),this->verboseAllProc());

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

            if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createMesh", "load msh file : " + this->mshfileStr(),
                                                this->worldComm(),this->verboseAllProc());

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
            M_mesh = GeoTool::createMeshFromGeoFile<mesh_type>(this->geofileStr(),this->prefix(),M_meshSize,1,
                                                               this->worldComm().localSize(),this->worldComm());
        }
        else
        {
            std::string geotoolSavePath;
            if ( this->geotoolSaveDirectory()!=this->appliShortRepository() )
            {
                if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createMesh", "change rep -> "+ this->geotoolSaveDirectory(),
                                                    this->worldComm(),this->verboseAllProc());
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
                if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createMesh", "change rep -> " + this->appliRepository() ,
                                                    this->worldComm(),this->verboseAllProc());
                Environment::changeRepository( _directory=boost::format(this->appliShortRepository()), _subdir=true );
            }

        }
        this->saveMSHfilePath(fmpath);
    }


    double timeElapsedCreateMesh = this->timerTool("Constructor").stop("createMesh");
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createMesh",
                                        (boost::format("finish in %1% s") % timeElapsedCreateMesh).str(),
                                        this->worldComm(),this->verboseAllProc());

} // createMesh()

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::createFunctionSpaces()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createFunctionSpaces", "start",
                                        this->worldComm(),this->verboseAllProc());
    this->timerTool("Constructor").start();

    // maybe build extended dof table
    std::vector<bool> extendedDT( space_fluid_type::nSpaces,false );
    if ( (this->doCIPStabConvection() || this->doCIPStabDivergence()) && !this->applyCIPStabOnlyOnBoundaryFaces() )
    {
        FeelModels::Log(this->prefix()+".FluidMechanics","createFunctionSpaces", "use buildDofTableMPIExtended on velocity",
                       this->worldComm(),this->verboseAllProc());
        extendedDT[0] = true;
    }
    if ( this->doCIPStabPressure() )
    {
        FeelModels::Log(this->prefix()+".FluidMechanics","createFunctionSpaces", "use buildDofTableMPIExtended on pressure",
                       this->worldComm(),this->verboseAllProc());
        extendedDT[1] = true;
    }
    // fluid mix space : velocity and pressure
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
        M_meshDirichletLM = createSubmesh(this->mesh(),markedfaces(this->mesh(),this->markerDirichletBClm()/*M_dirichletBCType["lm"]*/), useSubMeshRelationKey );
        if ( boption(_name="dirichletbc.lm.savemesh",_prefix=this->prefix()) )
        {
            std::string nameMeshDirichletLM = "nameMeshDirichletLM.msh";
            saveGMSHMesh(_mesh=M_meshDirichletLM,_filename=nameMeshDirichletLM);
        }

        M_XhDirichletLM = space_dirichletlm_velocity_type::New( _mesh=M_meshDirichletLM, _worldscomm=this->localNonCompositeWorldsComm() );
        //std::cout << "M_XhDirichletLM->nDof()"<< M_XhDirichletLM->nDof() <<std::endl;

    }

    this->timerTool("Constructor").stop("createSpaces");
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createFunctionSpaces", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::createTimeDiscretisation()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createTimeDiscretisation", "start",
                                        this->worldComm(),this->verboseAllProc());
    this->timerTool("Constructor").start();

    // bdf time schema
    std::string suffixName = "";
    if ( soption(_name="bdf.file-format",_prefix=this->prefix()) == "binary" )
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
                       _save=this->bdfSaveInFile(), _freq=this->bdfSaveFreq() );

    this->timerTool("Constructor").stop("createTimeDiscr");
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createTimeDiscretisation", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//


template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::createALE()
{
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() )
    {
        if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createALE", "start",
                                            this->worldComm(),this->verboseAllProc());
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
        if ( this->verbose() && moveGhostEltFromExtendedStencil )
            FeelModels::Log(this->prefix()+".FluidMechanics","createALE", "use moveGhostEltFromExtendedStencil",
                           this->worldComm(),this->verboseAllProc());

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
        if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createALE", "finish",
                                            this->worldComm(),this->verboseAllProc());
    }
#endif

}

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::createExporters()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createExporters", "start",
                                        this->worldComm(),this->verboseAllProc());
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

        if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createExporters", "start opInterpolation",
                                            this->worldComm(),this->verboseAllProc());
        boost::mpi::timer timerOpI;

        M_opIvelocity = opInterpolation(_domainSpace=M_Xh->template functionSpace<0>(),
                                        _imageSpace=M_XhVectorialVisuHO,
                                        _range=elements(M_XhVectorialVisuHO->mesh()),
                                        _backend=M_backend,
                                        _type=InterpolationNonConforme(false,true,false,15) );

        if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createExporters", "step1 done",
                                            this->worldComm(),this->verboseAllProc());

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

        if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createExporters", "step2 done",
                                            this->worldComm(),this->verboseAllProc());

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
        if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createExporters", "finish all opInterpolation in " + (boost::format("%1% s") % timeElapsedOpI).str(),
                                            this->worldComm(),this->verboseAllProc());

#if defined( FEELPP_MODELS_HAS_MESHALE )
        if (M_isMoveDomain) this->getMeshALE()->revertMovingMesh();
#endif

#endif
    }

    this->timerTool("Constructor").stop("createExporters");
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createExporters", "finish",
                                        this->worldComm(),this->verboseAllProc());

} // createExporters

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::createOthers()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createOthers", "start",
                                        this->worldComm(),this->verboseAllProc());
    this->timerTool("Constructor").start();
    //----------------------------------------------------------------------------//
    // rho, mu, nu with scalar P0 space
    std::vector<bool> extendedDT(1,this->useExtendedDofTable() );
    M_XhScalarP0 = space_scalar_P0_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm(),
                                              _extended_doftable=extendedDT );
    M_P0Rho.reset( new element_scalar_P0_type(M_XhScalarP0,"rho"));
    M_P0Mu.reset( new element_scalar_P0_type(M_XhScalarP0,"mu"));
    M_P0Nu.reset( new element_scalar_P0_type(M_XhScalarP0,"nu"));
    *M_P0Rho= vf::project(_space=M_XhScalarP0, _range=elements( M_mesh),
                          _expr=vf::cst(M_CstRho),_geomap=this->geomap());
    *M_P0Mu= vf::project(_space=M_XhScalarP0, _range=elements( M_mesh),
                         _expr=vf::cst(M_CstMu),_geomap=this->geomap());
    *M_P0Nu= vf::project(_space=M_XhScalarP0, _range=elements( M_mesh),
                         _expr=vf::cst(M_CstNu),_geomap=this->geomap());
    // viscosity model
    M_viscosityModelDesc.reset( new viscosity_model_type( this->stressTensorLawType(), *M_P0Mu, this->prefix() ) );
    //----------------------------------------------------------------------------//
    // space usefull to tranfert sigma*N()
    if (this->isMoveDomain()) this->createFunctionSpacesNormalStress();
    //----------------------------------------------------------------------------//
    // fluid outlet
#if 0
    if ( this->hasFluidOutlet() )
    {

        if ( this->fluidOutletType()=="windkessel")
        {
            M_fluidOutletWindkesselPressureDistal.resize(M_nFluidOutlet,0);
            M_fluidOutletWindkesselPressureDistal_old.resize(M_nFluidOutlet,std::vector<double>(Feel::BDF_MAX_ORDER));
            M_fluidOutletWindkesselPressureProximal.resize(M_nFluidOutlet,0);
            M_fluidOutletMarkerName.resize(M_nFluidOutlet);
            M_fluidOutletWindkesselRd.resize(M_nFluidOutlet);
            M_fluidOutletWindkesselRp.resize(M_nFluidOutlet);
            M_fluidOutletWindkesselCd.resize(M_nFluidOutlet);

            //std::list<std::string> markerNameBFOutlet;
            //ForEachBC( bcDef, cl::fluid_outlet,markerNameBFOutlet.push_back(PhysicalName));
            //M_fluidOutletsBCType["windkessel"]
            //std::string markerNameBFOutletBase = markerNameBFOutlet.front();

            auto itMarkerFluidOutlet = M_fluidOutletsBCType["windkessel"].begin();
            CHECK( itMarkerFluidOutlet != M_fluidOutletsBCType["windkessel"].end() ) << "no fluid outlet bc found\n";

            std::string markerNameBFOutletBase = *itMarkerFluidOutlet;

            // list usefull to create the outlets submesh
            std::list<std::string> markerNameBFOutletForSubmesh;

            bool useGenericMarkerName = (M_fluidOutletsBCType["windkessel"].size() == 1) && (M_nFluidOutlet > 1);

            if ( !useGenericMarkerName )
                CHECK( M_fluidOutletsBCType["windkessel"].size() == M_nFluidOutlet )  << "invalid nFluidOutlet : "
                                                                                      << M_fluidOutletsBCType["windkessel"].size() << " and "
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

                markerNameBFOutletForSubmesh.push_back(M_fluidOutletMarkerName[k]);

                M_fluidOutletWindkesselRd[k] = option(_name=(boost::format("fluid-outlet.windkessel.Rd%1%")%k).str(), _prefix=this->prefix()).as<double>();
                M_fluidOutletWindkesselRp[k] = option(_name=(boost::format("fluid-outlet.windkessel.Rp%1%")%k).str(), _prefix=this->prefix()).as<double>();
                M_fluidOutletWindkesselCd[k] = option(_name=(boost::format("fluid-outlet.windkessel.Cd%1%")%k).str(), _prefix=this->prefix()).as<double>();

                //std::cout << "M_fluidOutletMarkerName["<<k<<"] "<<M_fluidOutletMarkerName[k]<<std::endl;
                //std::cout << "M_fluidOutletWindkesselRd["<<k<<"] "<<M_fluidOutletWindkesselRd[k]<<std::endl;
                //std::cout << "M_fluidOutletWindkesselRp["<<k<<"] "<<M_fluidOutletWindkesselRp[k]<<std::endl;
                //std::cout << "M_fluidOutletWindkesselCd["<<k<<"] "<<M_fluidOutletWindkesselCd[k]<<std::endl;
            }
            if ( this->fluidOutletWindkesselCoupling() == "implicit" )
            {
                M_fluidOutletWindkesselMesh = createSubmesh( this->mesh(), markedfaces(this->mesh(),markerNameBFOutletForSubmesh) );
                M_fluidOutletWindkesselSpace = space_fluidoutlet_windkessel_type::New( _mesh=M_fluidOutletWindkesselMesh,
                                                                                       _worldscomm=std::vector<WorldComm>(2,this->worldComm()) );
                if ( M_isMoveDomain )
                {
                    M_fluidOutletWindkesselSpaceMeshDisp = space_fluidoutlet_windkessel_mesh_disp_type::New( _mesh=M_fluidOutletWindkesselMesh,
                                                                                                             _worldscomm=this->localNonCompositeWorldsComm() );

                    M_fluidOutletWindkesselMeshDisp = M_fluidOutletWindkesselSpaceMeshDisp->elementPtr();
                    M_fluidOutletWindkesselOpMeshDisp = opInterpolation(_domainSpace=M_meshALE->functionSpace(),
                                                                        _imageSpace=M_fluidOutletWindkesselSpaceMeshDisp,
                                                                        _range=elements(M_fluidOutletWindkesselMesh),
                                                                        _backend=M_backend );
                }

            }

        }



    } // if ( this->hasFluidOutlet() )
    else { M_nFluidOutlet=0; }
#else
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
#endif
    this->timerTool("Constructor").stop("createOthers");
    if (this->verbose()) FeelModels::Log(this->prefix()+".FluidMechanics","createOthers", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::createFunctionSpacesNormalStress()
{
    M_XhNormalBoundaryStress = space_stress_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    M_normalBoundaryStress.reset(new element_stress_type(M_XhNormalBoundaryStress));
#if defined( FEELPP_MODELS_HAS_MESHALE )
    M_normalStressFromStruct.reset(new element_stress_type(M_XhNormalBoundaryStress));
#endif
    M_wallShearStress.reset(new element_stress_type(M_XhNormalBoundaryStress));
}

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::createFunctionSpacesVorticity()
{
    M_Xh_vorticity = space_vorticity_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm());
    M_vorticity.reset( new element_vorticity_type(M_Xh_vorticity));
}

//---------------------------------------------------------------------------------------------------------//

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType, bool UsePeriodicity>
void
FluidMechanicsBase< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType,UsePeriodicity >::createFunctionSpacesSourceAdded()
{
    M_XhSourceAdded=space_vectorial_PN_type::New( _mesh=M_mesh );
    M_SourceAdded.reset( new element_vectorial_PN_type(M_XhSourceAdded,"SourceAdded"));
}

//---------------------------------------------------------------------------------------------------------//


} // namespace FeelModels
} // namespace Feel
