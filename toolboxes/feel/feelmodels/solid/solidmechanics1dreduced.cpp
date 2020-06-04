/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/solid/solidmechanics1dreduced.hpp>

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
{}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::init()
{

    if ( !M_mesh )
        this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    this->initFunctionSpaces();

    this->initBoundaryConditions();

}


SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
#if 0
    // axi-sym
    M_thickness_1dReduced = doption(_name="1dreduced-thickness",_prefix=this->prefix());
    M_radius_1dReduced = doption(_name="1dreduced-radius",_prefix=this->prefix());
#endif
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initMaterialProperties()
{}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initMesh()
{
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
    this->log("SolidMechanics","initFunctionSpaces", "start" );

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

    this->log("SolidMechanics","initFunctionSpaces", "finish" );

}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initTimeStep()
{
#if 0
        this->log("SolidMechanics","initTimeStep", "start" );
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
            M_newmark_displ_1dReduced = newmark( _space=M_Xh_1dReduced,
                                                 _name="displacement"+suffixName,
                                                 _prefix=this->prefix(),
                                                 _initial_time=ti, _final_time=tf, _time_step=dt,
                                                 _restart=this->doRestart(),_restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                                 _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
            M_newmark_displ_1dReduced->setfileFormat( myFileFormat );
            M_newmark_displ_1dReduced->setPathSave( ( saveTsDir/"displacement-1dreduced" ).string() );
        }
        else
        {
            CHECK( false ) << "invalid timeStepping : " << M_timeStepping;
        }

    if ( !this->doRestart() )
        {
            // up current time
            if ( M_timeStepping == "Newmark" )
                this->updateTime( M_newmark_displ_1dReduced->time() );
        }
        else
        {
            if ( M_timeStepping == "Newmark" )
            {
                // restart time step
                double tir = M_newmark_displ_1dReduced->restart();
                // load a previous solution as current solution
                *M_disp_1dReduced = M_newmark_displ_1dReduced->previousUnknown();
                this->updateInterfaceDispFrom1dDisp();
                // up initial time
                this->setTimeInitial( tir );
                // up current time
                this->updateTime( tir ) ;
            }
        }

#endif
}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initInitialConditions()
{}

SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::initPostProcess()
{
#if 0
        this->log("SolidMechanics","createExporters1dReduced", "start" );

    //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
    std::string geoExportType="static";

    if (!M_isHOVisu)
    {
        M_exporter_1dReduced = exporter( _mesh=this->mesh1dReduced(),
                                      //_name=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export-1dReduced")),
                                      _name="Export-1dReduced",
                                      _geo=geoExportType,
                                      _worldcomm=M_Xh_1dReduced->worldComm(),
                                      _path=this->exporterPath() );
    }

        // restart exporter
    if (this->doRestart() && this->restartPath().empty())
        if ( M_exporter_1dReduced && M_exporter_1dReduced->doExport() )
            M_exporter_1dReduced->restart(this->timeInitial());


    this->log("SolidMechanics","createExporters1dReduced", "finish" );

#endif
}


SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_1DREDUCED_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    //if ( M_timeStepping == "Newmark" )
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
    //if ( M_timeStepping == "Newmark" )
    {
        // next time step
        M_timeStepNewmark->next( *M_fieldDisp );
        // up current time
        this->updateTime( M_timeStepNewmark->time() );
    }
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
