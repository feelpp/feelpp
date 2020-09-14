/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/solid/solidmechanics1dreduced.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

namespace Feel
{
namespace FeelModels
{

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::SolidMechanics1dReduced( std::string const& prefix,
                                                                       std::string const& keyword,
                                                                       worldcomm_ptr_t const& worldComm,
                                                                       std::string const& subPrefix,
                                                                       ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelPhysics<nRealDim>( "solid" ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep )
{
    this->log("SolidMechanics1dReduced","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanics1dReducedConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanics1dReducedSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanics1dReducedPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".SolidMechanics1dReducedTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    // option in cfg files
    this->loadParameterFromOptionsVm();

    this->log("SolidMechanics1dReduced","constructor", "finish" );
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::init()
{

    if ( !M_mesh )
        this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    this->initFunctionSpaces();

    this->initBoundaryConditions();

    if ( !this->isStationary() )
        this->initTimeStep();

    this->initPostProcess();

    M_blockVectorSolution.resize( 1 );
    M_blockVectorSolution(0) = M_fieldDisp;

}


SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_timeStepping = "Newmark";
    // axi-sym
    M_thickness_1dReduced = doption(_name="1dreduced-thickness",_prefix=this->prefix());
    M_radius_1dReduced = doption(_name="1dreduced-radius",_prefix=this->prefix());

}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initMaterialProperties()
{}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("SolidMechanics1dReduced","initMesh", "start");
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh, prefixvm(this->prefix(),"Solid1dReducedMesh.path") /*this->fileNameMeshPath()*/);
    CHECK( M_mesh ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("initMesh");
    this->log("SolidMechanics1dReduced","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

#if 0
        this->log("SolidMechanics","initMesh1dReduced", "start" );
    std::string prefix1dreduced = prefixvm(this->prefix(),"1dreduced");

    std::string modelMeshRestartFile = prefixvm(this->prefix(),"SolidMechanics1dreducedMesh.path");
    std::string smpath = (fs::path( this->rootRepository() ) / fs::path( modelMeshRestartFile)).string();

    if (this->doRestart())
    {
        this->log("SolidMechanics","createMesh1dReduced", "reload mesh (because restart)" );

        if ( !this->restartPath().empty() )
            smpath = (fs::path( this->restartPath() ) / fs::path( modelMeshRestartFile)).string();

#if defined(SOLIDMECHANICS_1D_REDUCED_CREATESUBMESH)
        auto meshSM2dClass = reloadMesh<mesh_type>(smpath,this->worldCommPtr());
        SOLIDMECHANICS_1D_REDUCED_CREATESUBMESH(meshSM2dClass);
        M_mesh_1dReduced=mesh;
#else
        M_mesh_1dReduced = reloadMesh<mesh_1dreduced_type>(smpath,this->worldCommPtr());
#endif
    }
    else
    {
        if (Environment::vm().count(prefixvm(this->prefix(),"1dreduced-geofile")))
        {
            this->log("SolidMechanics","createMesh1dReduced", "use 1dreduced-geofile" );
            std::string geofile=soption(_name="1dreduced-geofile",_prefix=this->prefix() );
            std::string path = this->rootRepository();
            std::string mshfile = path + "/" + prefix1dreduced + ".msh";
            this->setMeshFile(mshfile);

            fs::path curPath=fs::current_path();
            bool hasChangedRep=false;
            if ( curPath != fs::path(this->rootRepository()) )
            {
                this->log("createMeshModel","", "change repository (temporary) for build mesh from geo : "+ this->rootRepository() );
                hasChangedRep=true;
                Environment::changeRepository( _directory=boost::format(this->rootRepository()), _subdir=false );
            }

            gmsh_ptrtype geodesc = geo( _filename=geofile,
                                        _prefix=prefix1dreduced,
                                        _worldcomm=this->worldCommPtr() );
            // allow to have a geo and msh file with a filename equal to prefix
            geodesc->setPrefix(prefix1dreduced);
            M_mesh_1dReduced = createGMSHMesh(_mesh=new mesh_1dreduced_type,_desc=geodesc,
                                              _prefix=prefix1dreduced,_worldcomm=this->worldCommPtr(),
                                              _partitions=this->worldComm().localSize() );

            // go back to previous repository
            if ( hasChangedRep )
                Environment::changeRepository( _directory=boost::format(curPath.string()), _subdir=false );
        }
        else
        {
            this->loadConfigMeshFile1dReduced( prefix1dreduced );
        }
        this->saveMeshFile( smpath );
    }

    this->log("SolidMechanics","initMesh1dReduced", "finish" );

#endif
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    this->log("SolidMechanics1dReduced","initFunctionSpaces", "start" );

    // function space and elements
    M_spaceDispVect = space_displacement_type::New(_mesh=M_mesh );
    M_spaceDisp = M_spaceDispVect->compSpace();
    // scalar field
    M_fieldDisp.reset( new element_displacement_component_type( M_spaceDisp, "structure displacement" ));
    M_fieldVelocity.reset( new element_displacement_component_type( M_spaceDisp, "structure velocity" ));
    M_fieldAcceleration.reset( new element_displacement_component_type( M_spaceDisp, "structure acceleration" ));
    // vectorial field
    M_fieldDisp_vect.reset(new element_displacement_type( M_spaceDispVect, "structure vect 1d displacement" ));
    M_fieldVelocity_vect.reset(new element_displacement_type( M_spaceDispVect, "velocity vect 1d displacement" ));

    // backend : use worldComm of Xh_1dReduced
    // M_backend_1dReduced = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh_1dReduced->worldCommPtr() );

    this->log("SolidMechanics1dReduced","initFunctionSpaces", "finish" );
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    std::string thekey = this->keyword();
    M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( std::move(thekey), "displacement_imposed" );
    thekey = this->keyword();
    M_volumicForcesProperties = this->modelProperties().boundaryConditions().getScalarFields( { { this->keyword(), "VolumicForces" } } );
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    this->log("SolidMechanics1dReduced","initTimeStep", "start" );
    this->timerTool("Constructor").start();

    double ti = this->timeInitial();
    double tf = this->timeFinal();
    double dt = this->timeStep();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    fs::path saveTsDir = fs::path(this->rootRepository())/fs::path( prefixvm(this->prefix(),prefixvm(this->subPrefix(),"ts")) );

    if ( M_timeStepping == "Newmark" )
    {
        M_timeStepNewmark = newmark( _space=M_spaceDisp,
                                     _name="displacement"+suffixName,
                                     _prefix=this->prefix(),
                                     _initial_time=ti, _final_time=tf, _time_step=dt,
                                     _restart=this->doRestart(),_restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                     _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
        M_timeStepNewmark->setfileFormat( myFileFormat );
        M_timeStepNewmark->setPathSave( ( saveTsDir/"displacement-1dreduced" ).string() );
    }
    else
    {
        CHECK( false ) << "invalid timeStepping : " << M_timeStepping;
    }

    if ( !this->doRestart() )
    {
        // up current time
        if ( M_timeStepping == "Newmark" )
            this->updateTime( M_timeStepNewmark->time() );
    }
    else
    {
        if ( M_timeStepping == "Newmark" )
        {
            // restart time step
            double tir = M_timeStepNewmark->restart();
            // load a previous solution as current solution
            *M_fieldDisp = M_timeStepNewmark->previousUnknown();
            this->updateInterfaceDispFrom1dDisp();
            // up initial time
            this->setTimeInitial( tir );
            // up current time
            this->updateTime( tir ) ;
        }
    }

    double tElapsed = this->timerTool("Constructor").stop("initTimeStep");
    this->log("SolidMechanics1dReduced","initTimeStep", (boost::format("finish in %1% s") %tElapsed).str() );
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initInitialConditions()
{}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("SolidMechanics1dReduced","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->setPostProcessExportsAllFieldsAvailable( {"displacement"} );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    this->setPostProcessSaveAllFieldsAvailable( {"displacement" } );
    super_type::initPostProcess();

    if ( !this->postProcessExportsFields().empty() )
    {
        std::string geoExportType="static";//change_coords_only, change, static
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export-1dReduced",
                               _geo=geoExportType,
                               _path=this->exporterPath() );

        // restart exporter
        if ( M_exporter->doExport() && this->doRestart() && this->restartPath().empty() )
            M_exporter->restart(this->timeInitial());
    }

    double tElpased = this->timerTool("Constructor").stop("initPostProcess");
    this->log("SolidMechanics1dReduced","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = 1;//this->nBlockMatrixGraph();

    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    int indexBlock=0;

    myblockGraph(indexBlock,indexBlock) =stencil(_test=this->functionSpace1dReduced(),_trial=this->functionSpace1dReduced(),
                                                 //_pattern_block=this->blockPattern(),
                                                 _diag_is_nonzero=(nBlock==1),
                                                 _close=(nBlock==1) )->graph();
    ++indexBlock;

    myblockGraph.close();
    return myblockGraph;
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    if ( M_timeStepping == "Newmark" )
    {
        // start time step
        if ( !this->doRestart() )
            M_timeStepNewmark->start( *M_fieldDisp );
        // up current time
        this->updateTime( M_timeStepNewmark->time() );
    }

}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    if ( M_timeStepping == "Newmark" )
    {
        // next time step
        M_timeStepNewmark->next( *M_fieldDisp );
        // up current time
        this->updateTime( M_timeStepNewmark->time() );
    }
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );

    this->setParameterValues( paramValues );
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }

    M_bcDirichlet.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );
}


SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::predictorDispl()
{
    //this->fieldDisplacementScal1dReduced().add(M_newmark_displ_1dReduced->timeStep(),M_newmark_displ_1dReduced->currentVelocity() );
    double dt = M_timeStepNewmark->timeStep();
    if (M_timeStepNewmark->iteration() == 1)
    {
        M_fieldDisp->add( dt, M_timeStepNewmark->previousVelocity(0) );
    }
    else
    {
        M_fieldDisp->add( (3./2.)*dt, M_timeStepNewmark->previousVelocity(0) );
        M_fieldDisp->add( (-1./2.)*dt, M_timeStepNewmark->previousVelocity(1) );
    }
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::updateVelocity()
{
    // if ( M_timeStepping != "Newmark" )
    //     return;

    M_timeStepNewmark->updateFromDisp(*M_fieldDisp);
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::updateInterfaceDispFrom1dDisp()
{
    M_fieldDisp_vect->on( _range=elements(this->mesh()),
                          _expr=idv(M_fieldDisp)*oneY() );
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::updateInterfaceVelocityFrom1dVelocity()
{
    M_fieldVelocity_vect->on( _range=elements(this->mesh()),
                              _expr=idv(M_timeStepNewmark->currentVelocity())*oneY() );
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
typename SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::element_displacement_ptrtype
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::extendVelocity1dReducedVectorial( element_displacement_component_type const& vel1d ) const
{
    auto res = M_spaceDispVect->elementPtr( vf::idv(vel1d)*vf::oneY() );
    return res;
}

//---------------------------------------------------------------------------------------------------//


} // namespace FeelModels
} // namespace Feel
