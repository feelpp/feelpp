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

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::SolidMechanicsBase( std::string __prefix,
                                                            bool __buildMesh,
                                                            WorldComm const& __worldComm,
                                                            std::string __subPrefix,
                                                            std::string __appliShortRepository )
    :
    super_type(__prefix,__worldComm,__subPrefix, __appliShortRepository),
    M_hasBuildFromMesh( false ), M_hasBuildFromMesh1dReduced( false ), M_isUpdatedForUse( false ),
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

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::build()
{
    /**/ if ( M_pdeType == "Elasticity" )                   { buildStandardModel(); }
    else if ( M_pdeType == "Elasticity-Large-Deformation" ) { buildStandardModel(); }
    else if ( M_pdeType == "Hyper-Elasticity" )             { buildStandardModel(); }
    else if ( M_pdeType == "Generalised-String" )           { build1dReducedModel(); }
    else
        CHECK(false) << "invalid pdeType" << M_pdeType;
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::build( mesh_ptrtype mesh )
{
    if ( M_pdeType == "Elasticity" || M_pdeType == "Elasticity-Large-Deformation" || M_pdeType == "Hyper-Elasticity" )
        this->buildStandardModel(mesh);
    else
        CHECK(false) << "invalid pdeType" << M_pdeType << "need standard model";
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::build( mesh_1dreduced_ptrtype mesh )
{
    if ( M_pdeType == "Generalised-String" )
        this->build1dReducedModel(mesh);
    else
        CHECK(false) << "invalid pdeType" << M_pdeType << " (need 1d_reduced model)";
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::buildStandardModel( mesh_ptrtype mesh )
{
    this->log("SolidMechanics","buildStandardModel", "start" );
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
    M_hasBuildFromMesh = true;
    this->log("SolidMechanics","buildStandardModel", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::loadMesh( mesh_ptrtype mesh )
{
    if (this->doRestart())
        this->build();
    else
        this->build(mesh);
}
SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::loadMesh( mesh_1dreduced_ptrtype mesh )
{
    // no restart case here :
    // we consider that the mesh given is identicaly (from createsubmesh of fluid mesh)
    this->build(mesh);
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::build1dReducedModel( mesh_1dreduced_ptrtype mesh )
{
    this->log("SolidMechanics","build1dReducedModel", "start" );
    //-----------------------------------------------------------------------------//
    M_isStandardModel=false;
    M_is1dReducedModel=true;
    //-----------------------------------------------------------------------------//
    if ( mesh )
        M_mesh_1dReduced = mesh;
    else
        this->createMesh1dReduced();
    //-----------------------------------------------------------------------------//
    this->createFunctionSpaces1dReduced();
    //-----------------------------------------------------------------------------//
    this->createTimeDiscretisation1dReduced();
    //-----------------------------------------------------------------------------//
    this->createExporters1dReduced();
    //-----------------------------------------------------------------------------//
    M_hasBuildFromMesh1dReduced = true;
    this->log("SolidMechanics","build1dReducedModel", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    this->log("SolidMechanics","loadParameterFromOptionsVm", "start" );

    M_meshSize = doption(_name="hsize",_prefix=this->prefix());
    M_useDisplacementPressureFormulation = boption(_name="use-incompressibility-constraint",_prefix=this->prefix());
    M_mechanicalProperties->setUseDisplacementPressureFormulation(M_useDisplacementPressureFormulation);
    this->pdeType( soption(_name="model",_prefix=this->prefix()) );
    if ( Environment::vm().count(prefixvm(this->prefix(),"solver").c_str()) )
        M_pdeSolver = soption(_name="solver",_prefix=this->prefix());
    M_useFSISemiImplicitScheme = false;
    M_couplingFSIcondition = "neumann";
    M_gammaNitschFSI = 2500;
    //M_penalbc = doption(_name="weakbccoeff",_prefix=this->prefix());
    M_isHOVisu = nOrderGeo > 1;
    if ( Environment::vm().count(prefixvm(this->prefix(),"hovisu").c_str()) )
        M_isHOVisu = boption(_name="hovisu",_prefix=this->prefix());
    M_doExportVelocity = boption(_name="do_export_velocity",_prefix=this->prefix());
    M_doExportAcceleration = boption(_name="do_export_acceleration",_prefix=this->prefix());
    M_doExportNormalStress = boption(_name="do_export_normalstress",_prefix=this->prefix());
    M_doExportVelocityInterfaceFromFluid = boption(_name="do_export_velocityinterfacefromfluid",_prefix=this->prefix());

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

    this->log("SolidMechanics","loadParameterFromOptionsVm", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createWorldsComm()
{
    this->log("SolidMechanics","createWorldsComm", "start" );

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

    this->log("SolidMechanics","createWorldsComm", "finish" );
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createMesh()
{
    this->log("SolidMechanics","createMesh", "start" );
    this->timerTool("Constructor").start();

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
                this->log( "SolidMechanics","createMesh",
                           "change repository (temporary) for build mesh from geo : "+ this->appliRepository() );
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
                this->log("SolidMechanics","createMesh", "change rep -> "+ this->geotoolSaveDirectory() );
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
                this->log("SolidMechanics","createMesh", "change rep -> " + this->appliRepository() );
                Environment::changeRepository( _directory=boost::format(this->appliShortRepository()), _subdir=true );
            }

        }
        this->saveMSHfilePath(smpath);
    }

    this->timerTool("Constructor").stop("createMesh");
    this->log("SolidMechanics","createMesh", "finish" );
} // createMesh()

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createMesh1dReduced()
{
    this->log("SolidMechanics","createMesh1dReduced", "start" );
    auto name = prefixvm(this->prefix(),"1dreduced");

    std::string smpath = prefixvm(this->prefix(),"SolidMechanics1dreducedMesh.path");
    if (this->doRestart())
    {
        this->log("SolidMechanics","createMesh1dReduced", "reload mesh (because restart)" );

        if ( !this->restartPath().empty() ) smpath = this->restartPath()+"/"+smpath;

#if defined(SOLIDMECHANICS_1D_REDUCED_CREATESUBMESH)
        auto meshSM2dClass = reloadMesh<mesh_type>(smpath,this->worldComm());
        SOLIDMECHANICS_1D_REDUCED_CREATESUBMESH(meshSM2dClass);
        M_mesh_1dReduced=mesh;
#else
        M_mesh_1dReduced = reloadMesh<mesh_1dreduced_type>(smpath,this->worldComm());
#endif
    }
    else
    {
        if (Environment::vm().count(prefixvm(this->prefix(),"1dreduced-geofile")))
        {
            this->log("SolidMechanics","createMesh1dReduced", "use 1dreduced-geofile" );
            std::string geofile=soption(_name="1dreduced-geofile",_prefix=this->prefix() );
            auto mesh = GeoTool::createMeshFromGeoFile<mesh_1dreduced_type>(geofile,name,M_meshSize);
            M_mesh_1dReduced = mesh;
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

    this->log("SolidMechanics","createMesh1dReduced", "finish" );
} // createMesh1dReduced


//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createFunctionSpaces()
{
    this->log("SolidMechanics","createFunctionSpaces", "start" );
    this->timerTool("Constructor").start();

    //--------------------------------------------------------//
    // function space for displacement (and maybe pressure)
    M_Xh = space_displacement_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    //--------------------------------------------------------//
    // displacement
    M_fieldDisplacement.reset( new element_displacement_type( M_Xh, "structure displacement" ));
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
    this->log("SolidMechanics","createFunctionSpaces", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createFunctionSpaces1dReduced()
{
    this->log("SolidMechanics","createFunctionSpaces1dReduced", "start" );

    // function space and elements
    M_Xh_vect_1dReduced = space_vect_1dreduced_type::New(_mesh=M_mesh_1dReduced,
                                                         _worldscomm=std::vector<WorldComm>(1,M_mesh_1dReduced->worldComm()));
    M_Xh_1dReduced = M_Xh_vect_1dReduced->compSpace();
    // scalar field
    M_disp_1dReduced.reset( new element_1dreduced_type( M_Xh_1dReduced, "structure displacement" ));
    M_velocity_1dReduced.reset( new element_1dreduced_type( M_Xh_1dReduced, "structure velocity" ));
    M_acceleration_1dReduced.reset( new element_1dreduced_type( M_Xh_1dReduced, "structure acceleration" ));
    // vectorial field
    M_disp_vect_1dReduced.reset(new element_vect_1dreduced_type( M_Xh_vect_1dReduced, "structure vect 1d displacement" ));
    M_velocity_vect_1dReduced.reset(new element_vect_1dreduced_type( M_Xh_vect_1dReduced, "velocity vect 1d displacement" ));

    // backend : use worldComm of Xh_1dReduced
    M_backend_1dReduced = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh_1dReduced->worldComm() );

    this->log("SolidMechanics","createFunctionSpaces1dReduced", "finish" );
}


//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createAdditionalFunctionSpacesFSI()
{
    if ( this->isStandardModel() )
        this->createAdditionalFunctionSpacesFSIStandard();
    else
        this->createAdditionalFunctionSpacesFSI1dReduced();
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createAdditionalFunctionSpacesFSIStandard()
{
    this->log("SolidMechanics","createAdditionalFunctionSpacesFSIStandard", "start" );

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


    this->log("SolidMechanics","createAdditionalFunctionSpacesFSIStandard", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createAdditionalFunctionSpacesFSI1dReduced()
{
    if ( M_XhStressVect_1dReduced ) return;
    this->log("SolidMechanics","createAdditionalFunctionSpacesFSI1dReduced", "start" );

    // normal stress as source term
    M_XhStressVect_1dReduced = space_stress_vect_1dreduced_type::New(_mesh=M_mesh_1dReduced,
                                                                       _worldscomm=std::vector<WorldComm>(1,M_mesh_1dReduced->worldComm()));
    M_stress_1dReduced.reset( new element_stress_scal_1dreduced_type( M_XhStressVect_1dReduced->compSpace(), "structure stress" ));
    M_stress_vect_1dReduced.reset(new element_stress_vect_1dreduced_type( M_XhStressVect_1dReduced, "stress 1d vect displacement" ));

    this->log("SolidMechanics","createAdditionalFunctionSpacesFSI1dReduced", "finish" );
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createTimeDiscretisation()
{
    this->log("SolidMechanics","createTimeDiscretisation", "start" );
    this->timerTool("Constructor").start();

    double ti = this->timeInitial();
    double tf = this->timeFinal();
    double dt = this->timeStep();

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
    this->log("SolidMechanics","createTimeDiscretisation", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createTimeDiscretisation1dReduced()
{
    this->log("SolidMechanics","createTimeDiscretisation1dReduced", "start" );

    double ti = this->timeInitial();
    double tf = this->timeFinal();
    double dt = this->timeStep();

    std::string suffixName = "";
    if ( soption(_name="ts.file-format",_prefix=this->prefix()) == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_newmark_displ_1dReduced = newmark( _vm=Environment::vm(),
                                          _space=M_Xh_1dReduced,
                                          _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"structure-1dreduced."+suffixName)),
                                          _prefix=this->prefix(),
                                          _initial_time=ti, _final_time=tf, _time_step=dt,
                                          _restart=this->doRestart(),_restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                          _save=this->bdfSaveInFile(), _freq=this->bdfSaveFreq() );
    M_newmark_displ_1dReduced->setPathSave( (fs::path(this->appliRepository()) /
                                             fs::path( prefixvm(this->prefix(), (boost::format("newmark_dt_%1%")%dt).str() ) ) ).string() );

    this->log("SolidMechanics","createTimeDiscretisation1dReduced", "finish" );
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createExporters()
{
    this->log("SolidMechanics","createExporters", "start" );
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
    this->log("SolidMechanics","createExporters", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createExporters1dReduced()
{
    this->log("SolidMechanics","createExporters1dReduced", "start" );

    //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
    std::string geoExportType="static";

    M_exporter_1dReduced = exporter( _mesh=this->mesh1dReduced(),
                                      //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export-1dReduced")),
                                      _name="Export-1dReduced",
                                      _geo=geoExportType,
                                      _worldcomm=M_Xh_1dReduced->worldComm(),
                                      _path=this->exporterPath() );

    this->log("SolidMechanics","createExporters1dReduced", "finish" );
}

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::createOthers()
{
    this->log(this->prefix()+".SolidMechanics","createOthers", "start" );
    this->timerTool("Constructor").start();

    M_XhScalarP0 = space_scalar_P0_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    M_mechanicalProperties->initFromSpace( M_XhScalarP0 );

    this->timerTool("Constructor").stop("createOthers");
    this->log("SolidMechanics","createOthers", "finish" );
}


//---------------------------------------------------------------------------------------------------//

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
template< typename TheBackendType >
NullSpace<double> extendNullSpace( NullSpace<double> const& ns,
                                   boost::shared_ptr<TheBackendType> const& mybackend,
                                   boost::shared_ptr<DataMap> const& dm )
{
    std::vector< typename NullSpace<double>::vector_ptrtype > myvecbasis(ns.size());
    for ( int k=0;k< ns.size();++k )
    {
        myvecbasis[k] = mybackend->newVector(dm);
        for( int i = 0 ; i < ns.basisVector(k).map().nLocalDofWithGhost() ; ++i )
            myvecbasis[k]->set(i, ns.basisVector(k)(i) );
        myvecbasis[k]->close();
    }
    NullSpace<double> userNullSpace( myvecbasis, mybackend );
    return userNullSpace;
}
} // detail

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::init( bool buildAlgebraicFactory, typename model_algebraic_factory_type::appli_ptrtype const& app )
{
    if ( M_isUpdatedForUse ) return;

    this->log("SolidMechanics","init", "start" );
    this->timerTool("Constructor").start();

    if ( (this->isStandardModel() && !M_hasBuildFromMesh ) ||
         (this->is1dReducedModel() && !M_hasBuildFromMesh1dReduced ) )
        this->build();

    if ( this->getMarkerNameFSI().size()>0 )
        this->createAdditionalFunctionSpacesFSI();

    // update timediscr and exporters
    if (!this->doRestart())
    {
        if (this->isStandardModel())
        {
            // start time step
            M_newmark_displ_struct->start(*M_fieldDisplacement);
            if ( M_useDisplacementPressureFormulation ) M_savetsPressure->start( M_XhPressure );
            // up current time
            this->updateTime( M_newmark_displ_struct->time() );
        }
        else if (this->is1dReducedModel())
        {
            // start time step
            M_newmark_displ_1dReduced->start(*M_disp_1dReduced);
            // up current time
            this->updateTime( M_newmark_displ_1dReduced->time() );
        }
    }
    else if (!this->isStationary()) // do a restart and transient mode
    {
        if (this->isStandardModel())
        {
            // restart time step
            M_newmark_displ_struct->restart();
            if ( M_useDisplacementPressureFormulation ) M_savetsPressure->restart();
            // load a previous solution as current solution
            *M_fieldDisplacement = M_newmark_displ_struct->previousUnknown();
            if ( M_useDisplacementPressureFormulation ) *M_fieldPressure = M_savetsPressure->unknown(0);
            // up initial time
            this->setTimeInitial(M_newmark_displ_struct->timeInitial());
            // restart exporter
            this->restartExporters();
            // up current time
            this->updateTime( M_newmark_displ_struct->time() );

        }
        else  if (this->is1dReducedModel())
        {
            // restart time step
            M_newmark_displ_1dReduced->restart();
            // load a previous solution as current solution
            *M_disp_1dReduced = M_newmark_displ_1dReduced->previousUnknown();
            // up initial time
            this->setTimeInitial(M_newmark_displ_1dReduced->timeInitial());
            // restart exporter
            this->restartExporters1dReduced();
            // up current time
            this->updateTime( M_newmark_displ_1dReduced->time() );
        }
    }

    this->log("SolidMechanics","init", "start/restart timeStep scheme done" );

    // update block vector (index + data struct)
    if (this->isStandardModel())
    {
        // define start dof index ( lm , windkessel )
        size_type currentStartIndex = 0;
        currentStartIndex += M_Xh->nLocalDofWithGhost();
        if ( M_useDisplacementPressureFormulation )
        {
            M_startDofIndexFieldsInMatrix["pressure"] = currentStartIndex;
            currentStartIndex += M_XhPressure->nLocalDofWithGhost() ;
        }
        // prepare block vector
        int nBlock = this->nBlockMatrixGraph();
        M_blockVectorSolution.resize( nBlock );
        M_blockVectorSolution(0) = this->fieldDisplacementPtr();
        int cptBlock=1;
        if ( M_useDisplacementPressureFormulation )
        {
            M_blockVectorSolution(cptBlock) = M_fieldPressure;
            ++cptBlock;
        }
        // init vector associated to the block
        M_blockVectorSolution.buildVector( this->backend() );
    }

    // update algebraic model
    if (buildAlgebraicFactory)
    {
        if (this->isStandardModel())
        {
            M_algebraicFactory.reset( new model_algebraic_factory_type(app,this->backend() ) );

            if ( this->nBlockMatrixGraph() == 1 )
            {
                NullSpace<double> userNullSpace = detail::getNullSpace(this->functionSpaceDisplacement(), mpl::int_<nDim>() ) ;
                if ( boption(_name="use-null-space",_prefix=this->prefix() ) )
                    M_algebraicFactory->attachNullSpace( userNullSpace );
                if ( boption(_name="use-near-null-space",_prefix=this->prefix() ) )
                    M_algebraicFactory->attachNearNullSpace( userNullSpace );
            }
            else
            {
                NullSpace<double> userNullSpace = detail::getNullSpace(this->functionSpaceDisplacement(), mpl::int_<nDim>() ) ;
                NullSpace<double> userNullSpaceFull = detail::extendNullSpace( userNullSpace, M_algebraicFactory->backend(), M_algebraicFactory->sparsityMatrixGraph()->mapRowPtr() );
                if ( boption(_name="use-near-null-space",_prefix=this->prefix() ) )
                {
                    M_algebraicFactory->attachNearNullSpace( 0,userNullSpace ); // for block disp in fieldsplit
                    M_algebraicFactory->attachNearNullSpace( userNullSpaceFull ); // for multigrid on full system
                }
            }

        }
        else if (this->is1dReducedModel())
        {
            M_algebraicFactory_1dReduced.reset( new model_algebraic_factory_type(app,this->backend1dReduced()) );
            M_algebraicFactory_1dReduced->initFromFunctionSpace( this->functionSpace1dReduced() );
        }
    }

    M_isUpdatedForUse = true;

    this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("SolidMechanics","init", "finish" );
}

//---------------------------------------------------------------------------------------------------//


} //FeelModels

} // Feel




