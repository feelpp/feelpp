/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels2/solid/solidmecbase.hpp>

#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>

#include <feel/feelmodels2/modelmesh/reloadmesh.hpp>


namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::SolidMechanicsBase( bool __isStationary,
                                                                             std::string __prefix,
                                                                             WorldComm const& __worldComm,
                                                                             bool __buildMesh,
                                                                             std::string __subPrefix,
                                                                             std::string __appliShortRepository )
    :
    super_type(__isStationary,__prefix,__worldComm,__subPrefix, __appliShortRepository),
    M_mechanicalProperties( new mechanicalproperties_type( __prefix ) )
{
    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanicsConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanicsSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanicsPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanicsTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::build()
{
    if ( M_pdeType == "Elasticity" )                      { buildStandardModel(); }
    else if ( M_pdeType == "Elasticity-Large-Deformation" ) { buildStandardModel(); }
    else if ( M_pdeType == "Hyper-Elasticity" )           { buildStandardModel(); }
    else if ( M_pdeType == "Generalised-String" )         { build1dReducedModel(); }
    else
        CHECK(false) << "invalid pdeType" << M_pdeType;
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::build( mesh_ptrtype mesh )
{
    if ( M_pdeType == "Elasticity" || M_pdeType == "Elasticity-Large-Deformation" || M_pdeType == "Hyper-Elasticity" )
        this->buildStandardModel(mesh);
    else
        CHECK(false) << "invalid pdeType" << M_pdeType << "need standard model";
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::build( mesh_1d_reduced_ptrtype mesh )
{
    if ( M_pdeType == "Generalised-String" )
        this->build1dReducedModel(mesh);
    else
        CHECK(false) << "invalid pdeType" << M_pdeType << " (need 1d_reduced model)";
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::buildStandardModel( mesh_ptrtype mesh )
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","initStandardModel", "start",
                                        this->worldComm(),this->verboseAllProc());
    //-----------------------------------------------------------------------------//
    M_isStandardModel=true;
    M_is1dReducedModel=false;
    //-----------------------------------------------------------------------------//
    // create or reload mesh
    if ( mesh )
        M_mesh = mesh;
    else
        this->createMesh();
    //-----------------------------------------------------------------------------//
    // functionSpaces and elements
    this->createFunctionSpaces();
    //-----------------------------------------------------------------------------//
    // time schema
    this->createTimeDiscretisation();
    //-----------------------------------------------------------------------------//
    // physical parameters
    this->createOthers();
    //-----------------------------------------------------------------------------//
    // exporters
    this->createExporters();
    //-----------------------------------------------------------------------------//
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","initStandardModel", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::loadMesh( mesh_ptrtype mesh )
{
    if (this->doRestart())
        this->build();
    else
        this->build(mesh);
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::build1dReducedModel( mesh_1d_reduced_ptrtype mesh )
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","init1dReducedModel", "start",
                                        this->worldComm(),this->verboseAllProc());
    //-----------------------------------------------------------------------------//
    M_isStandardModel=false;
    M_is1dReducedModel=true;
    //-----------------------------------------------------------------------------//
    if ( mesh )
        M_mesh_1d_reduced = mesh;
    else
        this->createMesh1dReduced();
    //-----------------------------------------------------------------------------//
    this->createFunctionSpaces1dReduced();
    //-----------------------------------------------------------------------------//
    this->createTimeDiscretisation1dReduced();
    //-----------------------------------------------------------------------------//
    this->createExporters1dReduced();
    //-----------------------------------------------------------------------------//
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","init1dReducedModel", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::loadParameterFromOptionsVm()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","loadParameterFromOptionsVm", "start",
                                        this->worldComm(),this->verboseAllProc());

    M_meshSize = doption(_name="hsize",_prefix=this->prefix());
    M_useDisplacementPressureFormulation = boption(_name="use-incompressibility-constraint",_prefix=this->prefix());
    M_mechanicalProperties->setUseDisplacementPressureFormulation(M_useDisplacementPressureFormulation);
    this->pdeType( soption(_name="model",_prefix=this->prefix()) );
    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        M_pdeSolver = soption(_name="solver",_prefix=this->prefix());
    M_useFSISemiImplicitScheme = false;
    M_couplingFSIcondition = "neumann";
    M_gammaNitschFSI = 2500;
    //M_weakCL = boption(_name="useweakbc",_prefix=this->prefix());
    M_penalbc = doption(_name="weakbccoeff",_prefix=this->prefix());
    M_isHOVisu = nOrderGeo > 1;
    if ( Environment::vm().count(prefixvm(this->prefix(),"hovisu").c_str()) )
        M_isHOVisu = boption(_name="hovisu",_prefix=this->prefix());
    M_doExportVelocity = boption(_name="do_export_velocity",_prefix=this->prefix());
    M_doExportAcceleration = boption(_name="do_export_acceleration",_prefix=this->prefix());
    M_doExportNormalStress = boption(_name="do_export_normalstress",_prefix=this->prefix());
    M_doExportVelocityInterfaceFromFluid = boption(_name="do_export_velocityinterfacefromfluid",_prefix=this->prefix());

    //-----------------------------------------------//
#if 0 // ASUP
    M_rho = doption(_name="rho",_prefix=this->prefix());// rho
    M_materialLaw = soption(_name="material_law",_prefix=this->prefix());
    M_youngmodulus = doption(_name="youngmodulus",_prefix=this->prefix());// E
    M_coeffpoisson = doption(_name="coeffpoisson",_prefix=this->prefix());// sigma
    this->updateLameCoeffFromYoungPoisson();
#endif
#if 0
    if (std::abs(0.5-M_coeffpoisson) > 1e-6 )
    {
        M_coefflame2 = M_youngmodulus/(2*(1+M_coeffpoisson));// mu
        M_coefflame1 = M_youngmodulus*M_coeffpoisson/((1+M_coeffpoisson)*(1-2*M_coeffpoisson));// lambda
    }
    else
    {
        M_coefflame2 = -M_coeffpoisson/M_youngmodulus; // mu
        M_coefflame1 = 0.5*(1+M_coeffpoisson)/M_youngmodulus; // lambda
    }
#endif
    //-----------------------------------------------//

    //time schema parameters
    std::string timeSchema = soption(_name="time-schema",_prefix=this->prefix());
    if (timeSchema == "Newmark")
    {
        M_genAlpha_alpha_m=1.0;
        M_genAlpha_alpha_f=1.0;
    }
    else if (timeSchema == "Generalized-Alpha")
    {
        M_genAlpha_rho=doption(_name="time-rho",_prefix=this->prefix());
        M_genAlpha_alpha_m=(2.- M_genAlpha_rho)/(1.+M_genAlpha_rho);
        M_genAlpha_alpha_f=1./(1.+M_genAlpha_rho);
    }
    else CHECK( false ) << "time scheme not supported : " << timeSchema << "\n";

    M_genAlpha_gamma=0.5+M_genAlpha_alpha_m-M_genAlpha_alpha_f;
    M_genAlpha_beta=0.25*(1+M_genAlpha_alpha_m-M_genAlpha_alpha_f)*(1+M_genAlpha_alpha_m-M_genAlpha_alpha_f);

    //-----------------------------------------------//

    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","loadParameterFromOptionsVm", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//
template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createWorldsComm()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createWorldsComm", "start",
                                        this->worldComm(),this->verboseAllProc());

    if (this->worldComm().localSize()==this->worldComm().globalSize())
    {
        std::vector<WorldComm> vecWorldComm(space_displacement_type::nSpaces,this->worldComm());
        std::vector<WorldComm> vecLocalWorldComm(1,this->worldComm());
        this->setWorldsComm(vecWorldComm);
        this->setLocalNonCompositeWorldsComm(vecLocalWorldComm);
    }
    else
    {
#if 0
        // manage world comm : WARNING only true without the lagrange multiplier
        const int DisplacementWorld=0;
        const int PressureWorld=1;
        int CurrentWorld=0;
        if (this->worldComm().globalRank() < this->worldComm().globalSize()/2 )
            CurrentWorld=DisplacementWorld;
        else
            CurrentWorld=PressureWorld;

        std::vector<WorldComm> vecWorldComm(2);
        std::vector<WorldComm> vecLocalWorldComm(1);
        if (this->worldComm().globalSize()>1)
        {
            vecWorldComm[0]=this->worldComm().subWorldComm(DisplacementWorld);
            vecWorldComm[1]=this->worldComm().subWorldComm(PressureWorld);
            vecLocalWorldComm[0]=this->worldComm().subWorldComm(CurrentWorld);
        }
        else
        {
            vecWorldComm[0]=WorldComm();
            vecWorldComm[1]=WorldComm();
            vecLocalWorldComm[0]=WorldComm();
        }
        this->setWorldsComm(vecWorldComm);
        this->setLocalNonCompositeWorldsComm(vecLocalWorldComm);
#else
        // non composite case
        std::vector<WorldComm> vecWorldComm(1,this->worldComm());
        std::vector<WorldComm> vecLocalWorldComm(1,this->worldComm());
        this->setWorldsComm(vecWorldComm);
        this->setLocalNonCompositeWorldsComm(vecLocalWorldComm);
#endif

    }

    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createWorldsComm", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createMesh()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createMesh", "start",
                                        this->worldComm(),this->verboseAllProc());
    this->timerTool("Constructor").start();

    //this->changeRepository();

    //auto smpath = this->fileNameMeshPath();//prefixvm(this->prefix(),"SolidMechanicsMesh.path");
    std::string smpath = (fs::path( this->appliRepository() ) / fs::path(this->fileNameMeshPath())).string();
    if (this->doRestart())
    {
        if ( !this->restartPath().empty() )
        {
            smpath = (fs::path( this->restartPath() ) / fs::path(this->fileNameMeshPath())).string();
        }
        M_mesh = reloadMesh<mesh_type>(smpath,this->worldComm());
    }
    else
    {
        if (this->hasMshfileStr())
        {
            std::string path = this->appliRepository();
            std::string mshfileRebuildPartitions = path + "/" + this->prefix() + ".msh";

            M_mesh = loadGMSHMesh(_mesh=new mesh_type,
                                  _filename=this->mshfileStr(),
                                  _worldcomm=this->worldComm(),
                                  _rebuild_partitions=this->rebuildMeshPartitions(),
                                  _rebuild_partitions_filename=mshfileRebuildPartitions,
                                  _partitions=this->worldComm().localSize(),
                                  _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );

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
                if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createMesh",
                                                     "change repository (temporary) for build mesh from geo : "+ this->appliRepository(),
                                                     this->worldComm(),this->verboseAllProc());
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
                if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createMesh", "change rep -> "+ this->geotoolSaveDirectory(),
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
            //std::string path = this->appliRepository();
            //std::string mshfile = path + "/" + this->prefix() + ".msh";
            this->setMshfileStr(mshfile);


            this->loadConfigMeshFile( geofilename );


            if ( this->geotoolSaveDirectory()!=this->appliShortRepository() )
            {
                if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createMesh", "change rep -> " + this->appliRepository() ,
                                                    this->worldComm(),this->verboseAllProc());
                Environment::changeRepository( _directory=boost::format(this->appliShortRepository()), _subdir=true );
            }

        }
        this->saveMSHfilePath(smpath);
    }

    this->timerTool("Constructor").stop("createMesh");
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createMesh", "finish",
                                        this->worldComm(),this->verboseAllProc());

} // createMesh()

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createMesh1dReduced()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createMesh1dReduced", "start",
                                        this->worldComm(),this->verboseAllProc());
    auto name = prefixvm(this->prefix(),"1dreduced");

    std::string smpath = prefixvm(this->prefix(),"SolidMechanics1dreducedMesh.path");
    if (this->doRestart())
    {
        if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createMesh1dReduced", "reload mesh (because restart)",
                                            this->worldComm(),this->verboseAllProc());

        if ( !this->restartPath().empty() ) smpath = this->restartPath()+"/"+smpath;

#if defined(SOLIDMECHANICS_1D_REDUCED_CREATESUBMESH)
        auto meshSM2dClass = reloadMesh<mesh_type>(smpath,this->worldComm());
        SOLIDMECHANICS_1D_REDUCED_CREATESUBMESH(meshSM2dClass);
        M_mesh_1d_reduced=mesh;
#else
        M_mesh_1d_reduced = reloadMesh<mesh_1d_reduced_type>(smpath,this->worldComm());
#endif
    }
    else
    {
#if 0
        // create a worldcomm only on master rank proc (no partition in 1d struct)
        std::vector<int> MapWorld(this->worldComm().globalSize(),0);
        MapWorld[this->worldComm().masterRank()] = 1;
        auto worldComm1dReduced = WorldComm(MapWorld,
                                            this->worldComm().localRank(),
                                            this->worldComm().globalComm(),
                                            this->worldComm().godComm()).subWorldComm(1);
#endif

        if (Environment::vm().count(prefixvm(this->prefix(),"1dreduced-geofile")))
        {
            if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createMesh1dReduced", "use 1dreduced-geofile",
                                                this->worldComm(),this->verboseAllProc());

            std::string geofile=soption(_name="1dreduced-geofile",_prefix=this->prefix() );
            auto mesh = GeoTool::createMeshFromGeoFile<mesh_1d_reduced_type>(geofile,name,M_meshSize);
            M_mesh_1d_reduced = mesh;
        }
        else
        {
            this->loadConfigMeshFile1dReduced( name );
        }
        // write msh file path
        std::ofstream file(smpath.c_str(), std::ios::out);
        std::string path = this->appliRepository();
        std::string mshfile = path + "/" + name + ".msh";
        file << mshfile;
        file.close();

    }

    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createMesh1dReduced", "finish",
                                        this->worldComm(),this->verboseAllProc());

} // createMesh1dReduced


//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createFunctionSpaces()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createFunctionSpaces", "start",
                                        this->worldComm(),this->verboseAllProc());
    this->timerTool("Constructor").start();

    //--------------------------------------------------------//
    // function space for displacement (and maybe pressure)
    M_Xh = space_displacement_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    //--------------------------------------------------------//
    // displacement
    U_displ_struct.reset( new element_displacement_type( M_Xh, "structure displacement" ));
    //--------------------------------------------------------//
    if ( M_useDisplacementPressureFormulation )
    {
        M_XhPressure = space_pressure_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
        M_fieldPressure.reset( new element_pressure_type( M_XhPressure, "pressure" ) );
    }
    //subfunctionspace vectorial
    M_XhVectorial = M_Xh;

    //--------------------------------------------------------//
    // pre-stress ( not functional )
    if (false)
        U_displ_struct_prestress.reset(new element_vectorial_type( M_XhVectorial, "structure displacement prestress" ));
    //--------------------------------------------------------//

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh->worldComm() );

    this->timerTool("Constructor").stop("createSpaces");
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createFunctionSpaces", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createFunctionSpaces1dReduced()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createFunctionSpaces1dReduced", "start",
                                        this->worldComm(),this->verboseAllProc());

    // function space and elements
    M_Xh_vect_1d_reduced = space_vect_1d_reduced_type::New(_mesh=M_mesh_1d_reduced,
                                                           _worldscomm=std::vector<WorldComm>(1,M_mesh_1d_reduced->worldComm()));
#if 0
    M_Xh_1d_reduced = space_1d_reduced_type::New(_mesh=M_mesh_1d_reduced,
                                                 _worldscomm=std::vector<WorldComm>(1,M_mesh_1d_reduced->worldComm()));
#else
    M_Xh_1d_reduced = M_Xh_vect_1d_reduced->compSpace();
#endif
    // scalar field
    M_disp_1d_reduced.reset( new element_1d_reduced_type( M_Xh_1d_reduced, "structure displacement" ));
    M_velocity_1d_reduced.reset( new element_1d_reduced_type( M_Xh_1d_reduced, "structure velocity" ));
    M_acceleration_1d_reduced.reset( new element_1d_reduced_type( M_Xh_1d_reduced, "structure acceleration" ));
    // vectorial field
    M_disp_vect_1d_reduced.reset(new element_vect_1d_reduced_type( M_Xh_vect_1d_reduced, "structure vect 1d displacement" ));
    M_velocity_vect_1d_reduced.reset(new element_vect_1d_reduced_type( M_Xh_vect_1d_reduced, "velocity vect 1d displacement" ));

#if 0
#if 0 // with intersection if 1d mesh is on middle surface
    M_map_stress_1d_reduced.reset(new std::vector<mesh_type::node_type>(M_disp_1d_reduced->nDof()));
    M_map_disp_1d_reduced.reset(new std::vector<mesh_type::node_type>(M_Xh->nDof()));
    auto bcDef = SOLIDMECHANICS_BC(this);
    ForEachBC( bcDef,cl::paroi_mobile,
               precomputeTransfertStress2dTo1d(PhysicalName);
               precomputeTransfertDisp1dTo2d(PhysicalName); );
#else // with interpolation if 1d mesh is on boundary
    //precomputeDisp1dTo2dWithInterpolation();
    //precomputeNormalStress2dTo1dWithInterpolation();
#endif
#endif

    // backend : use worldComm of Xh_1d_reduced
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh_1d_reduced->worldComm() );

    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createFunctionSpaces1dReduced", "finish",
                                        this->worldComm(),this->verboseAllProc());
}


//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createAdditionalFunctionSpacesFSI()
{
    if ( this->isStandardModel() )
        this->createAdditionalFunctionSpacesFSIStandard();
    else
        this->createAdditionalFunctionSpacesFSI1dReduced();
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createAdditionalFunctionSpacesFSIStandard()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createAdditionalFunctionSpacesFSIStandard", "start",
                                        this->worldComm(),this->verboseAllProc());

    //--------------------------------------------------------//
    // function space for normal stress
    if ( !M_XhStress )
        M_XhStress = space_stress_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    if ( !M_normalStressFromFluid )
        M_normalStressFromFluid.reset(new element_stress_type( M_XhStress, "normalStressBoundaryFromFluid" ));
    if ( !M_normalStressFromStruct )
        M_normalStressFromStruct.reset(new element_stress_type( M_XhStress, "normalStressBoundaryFromStruct" ));

    //--------------------------------------------------------//
    if ( this->couplingFSIcondition() == "robin" && !M_velocityInterfaceFromFluid )
        M_velocityInterfaceFromFluid.reset( new element_vectorial_type( M_XhVectorial, "velocityInterfaceFromFluid" ));


    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createAdditionalFunctionSpacesFSIStandard", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createAdditionalFunctionSpacesFSI1dReduced()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createAdditionalFunctionSpacesFSI1dReduced", "start",
                                        this->worldComm(),this->verboseAllProc());

    // normal stress as source term
    M_XhStressVect_1d_reduced = space_stress_vect_1d_reduced_type::New(_mesh=M_mesh_1d_reduced,
                                                                       _worldscomm=std::vector<WorldComm>(1,M_mesh_1d_reduced->worldComm()));
    M_stress_1d_reduced.reset( new element_stress_scal_1d_reduced_type( M_XhStressVect_1d_reduced->compSpace(), "structure stress" ));
    M_stress_vect_1d_reduced.reset(new element_stress_vect_1d_reduced_type( M_XhStressVect_1d_reduced, "stress 1d vect displacement" ));

    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createAdditionalFunctionSpacesFSI1dReduced", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createTimeDiscretisation()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createTimeDiscretisation", "start",
                                        this->worldComm(),this->verboseAllProc());
    this->timerTool("Constructor").start();

    auto ti = this->timeInitial();
    auto tf = this->timeFinal();
    auto dt = this->timeStep();

    std::string suffixName = "";
    if ( soption(_name="ts.file-format",_prefix=this->prefix()) == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_newmark_displ_struct = newmark( _vm=Environment::vm(), _space=M_Xh,
                                      _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"newmark"+suffixName)),
                                      _prefix=this->prefix(),
                                      _initial_time=ti, _final_time=tf, _time_step=dt,
                                      _restart=this->doRestart(), _restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                      _save=this->bdfSaveInFile(), _freq=this->bdfSaveFreq() );
    M_newmark_displ_struct->setPathSave( (fs::path(this->appliRepository()) /
                                          fs::path( prefixvm(this->prefix(), (boost::format("newmark_dt_%1%")%dt).str() ) ) ).string() );

    if ( M_useDisplacementPressureFormulation )
    {
        M_savetsPressure = bdf( _vm=Environment::vm(), _space=this->functionSpacePressure(),
                                _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"pressure"+suffixName)),
                                _prefix=this->prefix(),
                                _initial_time=ti, _final_time=tf, _time_step=dt,
                                _restart=this->doRestart(), _restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                _save=this->bdfSaveInFile(), _freq=this->bdfSaveFreq() );
        M_savetsPressure->setPathSave( (fs::path(this->appliRepository()) /
                                        fs::path( prefixvm(this->prefix(),"save-pressure" ) ) ).string() );
    }


    this->timerTool("Constructor").stop("createTimeDiscr");
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createTimeDiscretisation", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createTimeDiscretisation1dReduced()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createTimeDiscretisation1dReduced", "start",
                                        this->worldComm(),this->verboseAllProc());

    auto ti = this->timeInitial();
    auto tf = this->timeFinal();
    auto dt = this->timeStep();

    std::string suffixName = "";
    if ( soption(_name="ts.file-format",_prefix=this->prefix()) == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_newmark_displ_1d_reduced = newmark( _vm=Environment::vm(),
                                          _space=M_Xh_1d_reduced,
                                          _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"structure-1dreduced."+suffixName)),
                                          _prefix=this->prefix(),
                                          _initial_time=ti, _final_time=tf, _time_step=dt,
                                          _restart=this->doRestart(),_restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                          _save=this->bdfSaveInFile(), _freq=this->bdfSaveFreq() );

    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createTimeDiscretisation1dReduced", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createExporters()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createExporters", "start",
                                        this->worldComm(),this->verboseAllProc());
    this->timerTool("Constructor").start();

    // maybe need to build additional spaces
    if ( M_doExportNormalStress && !M_XhStress )
    {
        M_XhStress = space_stress_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
        M_normalStressFromFluid.reset(new element_stress_type( M_XhStress, "normalStressBoundaryFromFluid" ));
    }

    //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
    std::string geoExportType="static";//change_coords_only, change, static
    if (!M_isHOVisu)
    {
        M_exporter = exporter( _mesh=this->mesh(),
                               //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export")),
                               _name="Export",
                               _geo=geoExportType,
                               _worldcomm=M_Xh->worldComm(),
                               _path=this->exporterPath() );
    }
    else
    {
#if defined(FEELPP_HAS_VTK)
        boost::shared_ptr<mesh_visu_ho_type> meshVisuHO;
        std::string hovisuSpaceUsed = soption(_name="hovisu.space-used",_prefix=this->prefix());
        bool doLagP1parallel=false;
        if ( hovisuSpaceUsed == "displacement" )
        {
            //auto Xh_create_ho = space_create_ho_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
            auto Xh_create_ho = this->functionSpaceDisplacement()->compSpace();

            auto opLagP1 = lagrangeP1(_space=Xh_create_ho,
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
            CHECK( false ) << "not implement\n";
        }
        else if ( hovisuSpaceUsed == "p1" )
        {
            meshVisuHO = this->mesh()->createP1mesh();
        }
        else CHECK( false ) << "invalid hovisu.space-used " << hovisuSpaceUsed;

        M_exporter_ho = exporter( _mesh=meshVisuHO,
                                  //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"ExportHO")),
                                  _name="ExportHO",
                                  _geo=geoExportType,
                                  _worldcomm=M_Xh->worldComm(),
                                  _path=this->exporterPath() );

        M_XhVectorialVisuHO = space_vectorial_visu_ho_type::New( _mesh=meshVisuHO, _worldscomm=this->localNonCompositeWorldsComm());
        M_displacementVisuHO.reset( new element_vectorial_visu_ho_type(M_XhVectorialVisuHO,"u_visuHO"));

        M_opIdisplacement = opInterpolation(_domainSpace=this->functionSpaceDisplacement(),
                                            _imageSpace=M_XhVectorialVisuHO,
                                            _range=elements(M_XhVectorialVisuHO->mesh()),
                                            _backend=M_backend,
                                            _type=InterpolationNonConforme(false) );

        if (M_doExportNormalStress)
        {
            M_opInormalstress = opInterpolation(_domainSpace=this->M_normalStressFromFluid->functionSpace(),
                                                _imageSpace=M_XhVectorialVisuHO,
                                                //_range=elements(M_XhVectorialVisuHO->mesh()),
                                                _range=boundaryfaces(M_XhVectorialVisuHO->mesh()),
                                                _backend=M_backend,
                                                _type=InterpolationNonConforme(false) );
        }

        if ( M_useDisplacementPressureFormulation )
        {
            //M_XhScalarVisuHO = space_scalar_visu_ho_type::New(_mesh=opLagP1->mesh(), _worldscomm=this->localNonCompositeWorldsComm());
            M_XhScalarVisuHO = M_XhVectorialVisuHO->compSpace();
            M_pressureVisuHO.reset( new element_scalar_visu_ho_type(M_XhScalarVisuHO,"p_visuHO"));

            M_opIpressure = opInterpolation(_domainSpace=M_XhPressure,
                                            _imageSpace=M_XhScalarVisuHO,
                                            _range=elements(M_XhScalarVisuHO->mesh()),
                                            _backend=M_backend,
                                            _type=InterpolationNonConforme(false) );
        }
#endif // FEELPP_HAS_VTK
    }
    this->timerTool("Constructor").stop("createExporters");
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createExporters", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createExporters1dReduced()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createExporters1dReduced", "start",
                                        this->worldComm(),this->verboseAllProc());

    //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
    std::string geoExportType="static";

    M_exporter_1d_reduced = exporter( _mesh=this->mesh1dReduced(),
                                      //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export-1dReduced")),
                                      _name="Export-1dReduced",
                                      _geo=geoExportType,
                                      _worldcomm=M_Xh_1d_reduced->worldComm(),
                                      _path=this->exporterPath() );

    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createExporters1dReduced", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

template< typename ConvexType, int OrderDisp,bool UseCstMechProp >
void
SolidMechanicsBase<ConvexType,OrderDisp,UseCstMechProp>::createOthers()
{
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createOthers", "start",
                                        this->worldComm(),this->verboseAllProc());
    this->timerTool("Constructor").start();

    M_XhScalarP0 = space_scalar_P0_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    //M_P0Rho.reset( new element_scalar_P0_type(M_XhScalarP0,"rho"));
    //M_P0Coefflame1.reset( new element_scalar_P0_type(M_XhScalarP0,"coefflame1"));
    //M_P0Coefflame2.reset( new element_scalar_P0_type(M_XhScalarP0,"coefflame2"));
#if 0
    this->updateRho( vf::cst(M_rho) );
    this->updateCoefflame1( vf::cst(M_coefflame1) );
    this->updateCoefflame2( vf::cst(M_coefflame2) );
#endif
#if 0
    M_mechanicalProperties.reset( new mechanicalproperties_type( M_XhScalarP0,this->useDisplacementPressureFormulation(),this->prefix() ) );
    /*this->materialLaw(),M_useDisplacementPressureFormulation,
     this->youngModulus(),this->coeffPoisson(),*/
    /**M_P0Coefflame1,*M_P0Coefflame2*/
#else
    M_mechanicalProperties->initFromSpace( M_XhScalarP0 );
#endif

    this->timerTool("Constructor").stop("createOthers");
    if (this->verbose()) FeelModels::Log(this->prefix()+".SolidMechanics","createOthers", "finish",
                                        this->worldComm(),this->verboseAllProc());
}


//---------------------------------------------------------------------------------------------------//


} //FeelModels

} // Feel




